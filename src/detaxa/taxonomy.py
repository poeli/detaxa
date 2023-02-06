#!/usr/bin/env python

# Po-E (Paul) Li
# B-11, Los Alamos National Lab
# Date: 05/15/2016

import imp
import sys
import os
import tarfile
import requests
import logging

logger = logging.getLogger(__name__)

# Set default path of `taxonomy_db/` and `major_level_to_abbr.json`:
# The default `taxonomy_db/` path is your the location
lib_path = os.path.dirname(os.path.realpath(__file__))
taxonomy_dir = f"{lib_path}/taxonomy_db"
abbr_json_path = f"{taxonomy_dir}/major_level_to_abbr.json"

if os.path.isdir( "./taxonomy_db" ):
    taxonomy_dir = "./taxonomy_db"
    # if there is a new `major_level_to_abbr.json` in in the taxonomy_dir, use the json file.
    # Otherwise, use the default one comes with this package
    if os.path.isfile( f"{taxonomy_dir}/major_level_to_abbr.json" ):
        abbr_json_path = f"{taxonomy_dir}/major_level_to_abbr.json"
    else:
        pass
elif os.path.isdir( os.getcwd()+"/taxonomy_db" ):
    taxonomy_dir = os.getcwd()+"/taxonomy_db"

# init global dir
taxDepths      = {}
taxParents     = {}
taxRanks       = {}
taxNames       = {}
taxMerged      = {}
taxNumChilds   = {}
accTid         = {}
tidLineage     = {}
tidLineageDict = {}
nameTid        = {}

major_level_to_abbr = {}
abbr_to_major_level = {}


def taxid2rank( taxID, guess_strain=True ):
    taxID = _checkTaxonomy(taxID)
    if taxID == "unknown": return "unknown"

    if taxID == '1':
        return "root"

    if taxRanks[taxID] == "no rank" and guess_strain:
        # a leaf taxonomy is a strain
        if taxidIsLeaf(taxID):
            return "strain"
        # if not
        else:
            nmtid = taxid2nearestMajorTaxid(taxID)
            nmrank = _getTaxRank(nmtid)
            if nmrank == "species":
                return "species - others"
            else:
                return "others"
    
    return taxRanks[taxID]

def taxid2name(taxID):
    taxID = _checkTaxonomy(taxID)
    if taxID == "unknown":
        return "unknown"
    else:
        return _getTaxName(taxID)

def taxid2depth(taxID):
    taxID = _checkTaxonomy(taxID)
    if taxID == "unknown":
        return "unknown"
    else:
        return _getTaxDepth(taxID)

def taxid2type(taxID):
    taxID = _checkTaxonomy(taxID)
    if taxID == "unknown": return "unknown"

    origID = taxID
    lastID = taxID
    taxID = taxParents[taxID]

    while taxID != '1' and taxRanks[taxID] != 'species':
        lastID = taxID
        taxID = taxParents[taxID]

    if taxRanks[taxID] != 'species':
        taxID = 0
    else:
        taxID = lastID
        if taxID == origID: taxID = 0

    return taxID

def taxid2parent(taxID):
    taxID = _checkTaxonomy(taxID)
    if taxID == "unknown": return "unknown"

    taxID = taxParents[taxID]
    while taxID != '1' and taxRanks[taxID] == 'no rank':
        taxID = taxParents[taxID]

    return taxID

def name2taxid(name, rank=None, partial_match=False, method='all', reset=False):
    """
    name2taxid: convert organism name to taxid
    """
    global nameTid

    # reset previous mapping results
    if reset:
        nameTid = {}
    
    if not name in nameTid:
        matched_taxid = []
        for taxid in taxNames:
            if partial_match==True:
                if not name in taxNames[taxid]:
                    continue
            else:
                if name!=taxNames[taxid]:
                    continue
            
            if rank:
                if _getTaxRank(taxid)==rank:
                    matched_taxid.append(taxid)
            else:
                matched_taxid.append(taxid)

            # return when the first match found
            if method=='first' and len(matched_taxid):
                nameTid[name] = matched_taxid
                return nameTid[name]
        
        nameTid[name] = matched_taxid
        return nameTid[name]
    else:
        return nameTid[name]

def taxid2nameOnRank( taxID, target_rank=None ):
    taxID = _checkTaxonomy(taxID)
    if taxID == "unknown": return "unknown"
    if taxID == 1: return "root"
    if target_rank == "root": return "root"

    if taxID:
        rank = _getTaxRank(taxID)
        name = _getTaxName(taxID)

        if target_rank == "strain" and taxidIsLeaf(taxID):
            return name

        while taxID:
            if rank.upper() == target_rank.upper(): return name
            if name == 'root': break
            taxID = _getTaxParent(taxID)
            rank = _getTaxRank(taxID)
            name = _getTaxName(taxID)
    else:
        return ""

def taxid2taxidOnRank( taxID, target_rank=None ):
    taxID = _checkTaxonomy(taxID)
    if taxID == "unknown": return "unknown"

    if taxID:
        rank = _getTaxRank(taxID)
        name = _getTaxName(taxID)

        if target_rank == rank or ( target_rank == 'strain' and rank == 'no rank'): return taxID
        if target_rank == "root": return 1

        while taxID:
            if rank.upper() == target_rank.upper(): return taxID
            if name == 'root': break

            taxID = _getTaxParent(taxID)
            rank = _getTaxRank(taxID)
            name = _getTaxName(taxID)
    else:
        return ""

def taxidIsLeaf(taxID):
    taxID = _checkTaxonomy(taxID)
    if taxID == "unknown": return False
    if not taxID in taxNumChilds:
        return True
    else:
        return False

def _taxid2fullLink(taxID):
    taxID = _checkTaxonomy(taxID)
    if taxID == "unknown": return {}
    link = _autoVivification()

    while taxID != '1':
        name = _getTaxName(taxID)
        if not name: break
        
        parID = _getTaxParent(taxID)
        link[parID] = taxID
        taxID = parID
    
    return link

def taxid2fullLineage( taxID, sep='|', use_rank_abbr=False, space2underscore=True):
    link = _taxid2fullLink(taxID)
    texts = []
    if len(link):
        for p_taxID in link:
            taxID = link[p_taxID]
            rank = _getTaxRank(taxID)
            name = _getTaxName(taxID)

            if use_rank_abbr and (rank in major_level_to_abbr):
                rank =  major_level_to_abbr[rank]

            if sep == ';':
                texts.append(f"{rank}__{name}")
            else:
                texts.append(f"{rank}|{taxID}|{name}")
    
    texts.reverse()

    if space2underscore:
        return sep.join(texts).replace(' ', '_')
    else:
        return sep.join(texts)

def taxid2fullLinkDict(taxID):
    return _taxid2fullLink( taxID)

def taxid2nearestMajorTaxid(taxID):
    taxID = _checkTaxonomy(taxID)
    if taxID == "unknown": return "unknown"
    ptid = _getTaxParent(taxID)
    while ptid != '1':
        tmp = _getTaxRank( ptid )
        if tmp in major_level_to_abbr:
            return ptid
        else:
            ptid = _getTaxParent( ptid )

    return "1"

def taxid2lineage( tid, all_major_rank=True, print_strain=False, space2underscore=False, sep="|"):
    lineage = _taxid2lineage( tid, all_major_rank, print_strain, space2underscore)
    texts = []
    for rank in major_level_to_abbr:
        if rank in lineage:
            if print_strain==False and rank=="strain":
                continue
            if sep == ";":
                texts.append( f"{major_level_to_abbr[rank]}__{lineage[rank]['name']}" )
            else:
                texts.append( f"{rank}|{lineage[rank]['taxid']}|{lineage[rank]['name']}" ) 
    
    if space2underscore:
        return sep.join(texts).replace(' ', '_')
    else:
        return sep.join(texts) 

def taxid2lineageDICT( tid, all_major_rank=True, print_strain=False, space2underscore=False):
    return _taxid2lineage( tid, all_major_rank, print_strain, space2underscore)

def _taxid2lineage(tid, all_major_rank, print_strain, space2underscore):
    taxID = _checkTaxonomy( tid )
    if taxID == "unknown": return {}
    if taxID in tidLineageDict: return tidLineageDict[taxID]

    info = _autoVivification()
    level = {abbr: '' for abbr in abbr_to_major_level}

    rank = _getTaxRank(taxID)
    orig_rank = rank
    name = _getTaxName(taxID)
    str_name = name
    if space2underscore: str_name = str_name.replace(" ", "_")

    while taxID:
        if rank in major_level_to_abbr:
            if space2underscore: name = name.replace(" ", "_")
            level[major_level_to_abbr[rank]] = name

            #for output JSON
            info[rank]["name"] = name
            info[rank]["taxid"] = taxID

        taxID = _getTaxParent(taxID)
        rank = _getTaxRank(taxID)
        name = _getTaxName(taxID)

        if name == 'root': break

    # try to get the closest "no_rank" taxa to "type" representing subtype/group (mainly for virus)
    typeTID = taxid2type(tid)
    if typeTID:
        info["type"]["name"]  = _getTaxName(typeTID)
        info["type"]["taxid"] = typeTID

    last = str_name

    ranks = list(abbr_to_major_level.keys())
    ranks.reverse()
    idx = 0
    
    # input taxid is a major rank
    if orig_rank in major_level_to_abbr:
        idx = ranks.index( major_level_to_abbr[orig_rank] )
    # if not, find the next major rank
    else:
        nmtid = taxid2nearestMajorTaxid( tid )
        nmrank = taxid2rank( nmtid )
        if nmrank == "root":
            idx = 7
        else:
            idx = ranks.index( major_level_to_abbr[nmrank] )

    for lvl in ranks[idx:]:
        if all_major_rank == False:
            if not level[lvl]: continue

        if not level[lvl]:
            level[lvl] = "%s - no_%s_rank"%(last,lvl)
            info[abbr_to_major_level[lvl]]["name"]  = "%s - no_%s_rank"%(last,lvl)
            info[abbr_to_major_level[lvl]]["taxid"] = 0

        last=level[lvl]

    if orig_rank == "strain":
        info["strain"]["name"]  = str_name
        info["strain"]["taxid"] = tid

    tidLineageDict[tid] = info
    return info

def _getTaxDepth(taxID):
    if taxID in taxMerged: taxID = taxMerged[taxID]
    return taxDepths[taxID]

def _getTaxName(taxID):
    if taxID in taxMerged: taxID = taxMerged[taxID]
    return taxNames[taxID]

def _getTaxParent(taxID):
    if taxID in taxMerged: taxID = taxMerged[taxID]
    return taxParents[taxID]

def _getTaxRank(taxID):
    if taxID in taxMerged: taxID = taxMerged[taxID]
    return taxRanks[taxID]

def lca_taxid(taxids):
    """ lca_taxid
    Return lowest common ancestor (LCA) taxid of input taxids
    """
    ranks = ['strain','species','genus','family','order','class','phylum','superkingdom']

    merged_dict = _autoVivification()
    for tid in taxids:
        lineage = taxid2lineageDICT(tid, 1, 1)
        for r in ranks:
            if not r in lineage:
                ttid = "0"
            else:
                ttid = lineage[r]['taxid']

            if ttid in merged_dict[r]:
                merged_dict[r][ttid] += 1
            else:
                merged_dict[r][ttid] = 1

    for r in ranks:
        if len(merged_dict[r]) == 1:
            for ttid in merged_dict[r]:
                # skip if no tid in this rank
                if ttid=="0":
                    continue
                return ttid

    return '1'

def acc2taxid( acc,  accession2taxid_file=None):
    if not accession2taxid_file:
        accession2taxid_file = f"{taxonomy_dir}/accession2taxid.tsv"

    #remove version number#
    acc = acc.split('.')[0]

    if not acc in accTid:
        logger.info( f"acc2taxid from file: {accession2taxid_file}" )
        with open( accession2taxid_file ) as f:
            f.seek(0, 2)
            start = 0
            end = f.tell()
            accCur = ""
            
            while( acc != accCur and start < end ):
                posNew = (end+start)/2
                f.seek( posNew )
                if posNew != start: f.readline()
                line = f.readline()    
                
                logger.debug( "start: %15d, posNew: %15d, end: %15d, line: %s" % (start, posNew, end, line) )
                if line :
                    (accNew, tid) = line.split('\t')
                else:
                    break

                logger.debug( f'[acc, accNew, accCur]=[{acc}, {accNew}, {accCur}]')
                if acc > accNew and accCur != accNew and accNew:
                    if accNew: posNew = f.tell()
                    start = posNew
                    if start >= end: end = start+1
                else:
                    end = posNew
                
                accCur = accNew
            f.close()

            if accCur == acc:
                accTid[acc] = tid.strip()
            else:
                accTid[acc] = ""

    return accTid[acc]

def loadTaxonomy(dbpath=None,
                 cus_taxonomy_file=None, 
                 cus_taxonomy_format="tsv",
                 auto_download=True):
    """ loadTaxonomy
    Method to load taxonomy

    Arguments:
    dbpath [STR] The path to search NCBI taxonomy files (default None)
    cus_taxonomy_file [STR] The path of the custom taxonomy file (default None)
    cus_taxonomy_format ['tsv','mgnify_lineage','gtdb_taxonomy','gtdb_metadata'] The format of the custom taxonomy file. (default 'tsv')
    auto_download [BOOL] Download NCBI taxonomy files when local database not found. (default True)

    Return: None
    """
    global taxonomy_dir, abbr_json_path
    
    if dbpath:
        taxonomy_dir = dbpath

    logger.debug( f"Taxonomy directory: {taxonomy_dir}" )

    # loading major levels to json file
    _loadAbbrJson(abbr_json_path)

    #NCBI ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
    taxdump_tgz_file = taxonomy_dir+"/taxdump.tar.gz"

    #raw taxonomy dmp files from NCBI
    names_dmp_file = taxonomy_dir+"/names.dmp"
    nodes_dmp_file = taxonomy_dir+"/nodes.dmp"
    merged_dmp_file = taxonomy_dir+"/merged.dmp"

    #parsed taxonomy tsv file
    taxonomy_file = taxonomy_dir+"/taxonomy.tsv"
    merged_taxonomy_file = taxonomy_dir+"/taxonomy.merged.tsv"

    #custom taxonomy file
    if not cus_taxonomy_file:
        cus_taxonomy_file = taxonomy_dir+"/taxonomy.custom.tsv"

    # checking if taxonomy files provided
    if not os.path.isfile( taxdump_tgz_file ) \
        and not os.path.isfile( merged_taxonomy_file ) \
        and not os.path.isfile( taxonomy_file ) \
        and not (os.path.isfile( names_dmp_file ) and os.path.isfile( nodes_dmp_file )) \
        and not os.path.isfile( cus_taxonomy_file ):
        
        logger.info( f"No taxonomy files not found." )
        if auto_download:
            NCBITaxonomyDownload(taxonomy_dir)
        else:
            logger.info( f"Auto-download is off." )
            logger.fatal( f"No available taxonomy files." )
            _die( "[ERROR] No available taxonomy files." )

    # try to load taxonomy from taxonomy.tsv
    if os.path.isfile( nodes_dmp_file ) and  os.path.isfile( names_dmp_file ):
        loadNCBITaxonomy(taxdump_tgz_file, names_dmp_file, nodes_dmp_file, merged_dmp_file)
    elif os.path.isfile(taxdump_tgz_file):
        loadNCBITaxonomy(taxdump_tgz_file, names_dmp_file, nodes_dmp_file, merged_dmp_file)

    if os.path.isfile(taxonomy_file):
        logger.info( "Open taxonomy file: %s"% taxonomy_file )
        loadTaxonomyTSV(taxonomy_file)

    # try to load custom taxonomy from taxonomy.custom.tsv
    if os.path.isfile(cus_taxonomy_file) and (cus_taxonomy_format=='tsv'):
        logger.info( "Open custom taxonomy node file (tsv format): %s"% cus_taxonomy_file)
        loadTaxonomyTSV(cus_taxonomy_file)
    # try to load custom taxonomy from lineage file
    elif os.path.isfile(cus_taxonomy_file) and (cus_taxonomy_format=='mgnify_lineage'):
        logger.info( "Open custom taxonomy node file (lineage format): %s"% cus_taxonomy_file)
        loadMgnifyTaxonomy(cus_taxonomy_file)
    # try to load custom taxonomy from GTDB file
    elif os.path.isfile(cus_taxonomy_file) and (cus_taxonomy_format in ['gtdb_taxonomy','gtdb_metadata']):
        loadGTDBTaxonomy(cus_taxonomy_file, cus_taxonomy_format)
    elif os.path.isfile(cus_taxonomy_file):
        logger.fatal( f"invalid cus_taxonomy_format: {cus_taxonomy_format}" )
        _die(f"[ERROR] Invalid cus_taxonomy_format: {cus_taxonomy_format}")


def _loadAbbrJson(abbr_json_path):
    import json
    global major_level_to_abbr, abbr_to_major_level

    # Opening JSON file
    with open(abbr_json_path) as f:    
        major_level_to_abbr = json.load(f)
        f.close

    if len(major_level_to_abbr):
        abbr_to_major_level = {v: k for k, v in major_level_to_abbr.items()}
    else:
        logger.fatal( f"None of the major level to aberration loaded from {abbr_json_path}." )
        _die(f"[ERROR] None of the major level to aberration loaded from {abbr_json_path}.")


def NCBITaxonomyDownload(dir=None):
    global taxonomy_dir

    if not dir:
        dir = taxonomy_dir
        logger.info( f"No destination taxonomy dir input. Use default: {dir}..." )

    url = 'http://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
    logger.info( f"Auto downloading taxanomy from {url}..." )
    # download taxonomy file if auto_download enabled
    r = requests.get(url)
    if not os.path.exists( taxonomy_dir ):
        os.makedirs( taxonomy_dir )
    
    taxdump_tgz_file = f'{taxonomy_dir}/taxdump.tar.gz'

    with open(taxdump_tgz_file, 'wb') as f:
        f.write(r.content)
    if os.path.getsize( taxdump_tgz_file ):
        logger.info( f"Saved to {taxdump_tgz_file}." )
    else:
        logger.fatal( f"Failed to download or save taxonomy files." )
        _die( "[ERROR] Failed to download or save taxonomy files." )    

    # extract
    tax_tar = tarfile.open(taxdump_tgz_file, "r:gz")
    logger.info( f"Extracting nodes.dmp..." )
    tax_tar.extract('nodes.dmp', taxonomy_dir)
    logger.info( f"Extracting names.dmp..." )
    tax_tar.extract('names.dmp', taxonomy_dir)
    logger.info( f"Extracting merged.dmp..." )
    tax_tar.extract('merged.dmp', taxonomy_dir)
    tax_tar.close()
    # delete taxdump_tgz_file
    os.remove(taxdump_tgz_file)

def loadTaxonomyTSV(tsv_taxonomy_file):

    # loading major levels to json file
    _loadAbbrJson(abbr_json_path)

    try:
        with open(tsv_taxonomy_file) as f:
            for line in f:
                line = line.rstrip('\r\n')
                if not line: continue
                tid, depth, parent, rank, name = line.split('\t')
                taxParents[tid] = parent
                taxDepths[tid] = depth
                taxRanks[tid] = rank
                taxNames[tid] = name
                if parent in taxNumChilds:
                    taxNumChilds[parent] += 1
                else:
                    taxNumChilds[parent] = 1
            f.close()
            logger.info( f"Done parsing custom tsv taxonomy file." )
    except IOError:
        _die( "Failed to open custom tsv taxonomy file: %s." % tsv_taxonomy_file )

def loadNCBITaxonomy(taxdump_tgz_file=None, 
                     names_dmp_file=None, 
                     nodes_dmp_file=None, 
                     merged_dmp_file=None):

    # loading major levels to json file
    _loadAbbrJson(abbr_json_path)

    # try to load taxonomy from taxonomy.tsv
    if os.path.isfile( nodes_dmp_file ) and  os.path.isfile( names_dmp_file ):
        try:
            # read name from names.dmp
            logger.info( f"Open taxonomy name file: {names_dmp_file}" )
            with open(names_dmp_file) as f:
                for line in f:
                    tid, name, tmp, nametype = line.rstrip('\r\n').split('\t|\t')
                    if not nametype.startswith("scientific name"):
                        continue
                    taxNames[tid] = name
                f.close()
                logger.info( f"Done parsing taxonomy name file." )    

            # read taxonomy info from nodes.dmp
            logger.info( f"Open taxonomy node file: {nodes_dmp_file}" )
            with open(nodes_dmp_file) as f:
                for line in f:
                    fields = line.rstrip('\r\n').split('\t|\t')
                    tid = fields[0]
                    parent = fields[1]
                    taxParents[tid] = parent
                    taxDepths[tid] = taxDepths[parent]+1 if parent in taxDepths else 0 # could have potiential bug if child node is parsed before parent node.
                    taxRanks[tid] = fields[2]
                    if parent in taxNumChilds:
                        taxNumChilds[parent] += 1
                    else:
                        taxNumChilds[parent] = 1
                f.close()
                logger.info( f"Done parsing taxonomy node file." )
        except IOError:
            _die( "Failed to open taxonomy files (taxonomy.tsv, nodes.dmp and names.dmp)." )
    elif os.path.isfile( taxdump_tgz_file ):
        try:
            logger.info( f"Open taxonomy file: {taxdump_tgz_file}" )
            tar = tarfile.open(taxdump_tgz_file, "r:gz")
            
            # read name from names.dmp
            logger.info( "Extract taxonomy names file: names.dmp" )
            member = tar.getmember("names.dmp")
            f = tar.extractfile(member)
            for line in f.readlines():
                tid, name, tmp, nametype = line.decode('utf8').rstrip('\r\n').split('\t|\t')
                if not nametype.startswith("scientific name"):
                    continue
                taxNames[tid] = name
            f.close()
            
            # read taxonomy info from nodes.dmp
            logger.info( "Extract taxonomy nodes file: nodes.dmp" )
            member = tar.getmember("nodes.dmp")
            f = tar.extractfile(member)
            for line in f.readlines():
                fields = line.decode('utf8').rstrip('\r\n').split('\t|\t')
                tid = fields[0]
                parent = fields[1]
                taxParents[tid] = parent
                taxDepths[tid] = taxDepths[parent]+1 if parent in taxDepths else 0 # could have potiential bug if child node is parsed before parent node.
                taxRanks[tid] = fields[2]
                if parent in taxNumChilds:
                    taxNumChilds[parent] += 1
                else:
                    taxNumChilds[parent] = 1
            f.close()

            # read taxonomy info from merged.dmp
            logger.info( "Extract taxonomy merged file: merged.dmp" )
            member = tar.getmember("merged.dmp")
            f = tar.extractfile(member)
            for line in f.readlines():
                fields = line.decode('utf8').rstrip('\r\n').split('\t|')
                taxMerged[fields[0]] = fields[1].strip('\t')
            f.close()

        except IOError:
            _die( "Failed to load taxonomy from %s"%taxdump_tgz_file )
    
    #try to load merged taxids
    if os.path.isfile( merged_dmp_file ):
        logger.info( "Open merged taxonomy node file: %s"% merged_dmp_file )
        with open(merged_dmp_file) as f:
            for line in f:
                fields = line.rstrip('\r\n').split('\t|')
                taxMerged[fields[0]] = fields[1].strip('\t')
            f.close()
            logger.info( f"Done parsing merged taxonomy file." )
    
def loadMgnifyTaxonomy(mgnify_taxonomy_file=None):
    """
    loadMgnifyTaxonomy()
    """

    # loading major levels to json file
    _loadAbbrJson(abbr_json_path)

    # try to load custom taxonomy from lineage file
    if os.path.isfile(mgnify_taxonomy_file):
        logger.info( "Open custom taxonomy node file (lineage format): %s"% mgnify_taxonomy_file)
        try:
            with open(mgnify_taxonomy_file) as f:
                for line in f:
                    line = line.rstrip('\r\n')
                    if not line: continue
                    if line.startswith('#'): continue
                    if not line.startswith('sk__'):
                        logger.warn( f"A text line of lineage has to start with 'sk__'...skipped: {line}" )
                        continue

                    temp = line.split(';')
                    p_name = ''
                    rank = ''
                    name = ''

                    for i in range(1, len(temp)+1):
                        # this taxa
                        (rank_abbr, name) = temp[-i].split('__')

                        import re
                        re_taxa = re.compile("^([^_]+)__(.*)$")
                        rank_abbr, name = re_taxa.match(temp[-i]).groups()

                        # for na taxon (no_{rank_abbr}_rank)
                        if name=="": name = p_name

                        # paranet taxa
                        try:
                            p_rank_abbr, p_name = re_taxa.match(temp[-(i+1)]).groups()
                            if p_name=="":
                                p_name = f'{name} - no_{rank_abbr}_rank'
                        except:
                            # for the superkingdom rank, assign parant taxid to 1 (root)
                            p_name = '1'
                            if not '1' in taxRanks: taxRanks['1'] = 'root'
                            if not '1' in taxNames: taxNames['1'] = 'root'

                        if rank_abbr in abbr_to_major_level:
                            rank = abbr_to_major_level[rank_abbr]
                        else:
                            rank = rank_abbr
                            
                        tid = name
                        taxParents[tid] = p_name
                        taxRanks[tid] = rank
                        taxNames[tid] = name
                        if p_name in taxNumChilds:
                            taxNumChilds[p_name] += 1
                        else:
                            taxNumChilds[p_name] = 1
                f.close()
                logger.info( f"Done parsing custom taxonomy file." )
        except IOError:
            _die( "Failed to open custom taxonomy file: %s." % mgnify_taxonomy_file )

    logger.info( f"Done parsing taxonomy files (total {len(taxParents)} taxa loaded)" )

def loadGTDBTaxonomy(gtdb_taxonomy_file=None, gtdb_taxonomy_format="gtdb_metadata"):
    """
    loadGTDBTaxonomy()
    """

    # loading major levels to json file
    _loadAbbrJson(abbr_json_path, type='gtdb')

    # try to load custom taxonomy from GTDB file
    if os.path.isfile(gtdb_taxonomy_file) and (gtdb_taxonomy_format in ['gtdb_taxonomy','gtdb_metadata']):
        logger.info( f"Open custom taxonomy node file ({gtdb_taxonomy_format}): %s"% gtdb_taxonomy_file)
        try:
            with open(gtdb_taxonomy_file) as f:
                for line in f:
                    line = line.rstrip('\r\n')
                    if not line: continue
                    if line.startswith('#'): continue
                    if line.startswith('accession'): continue

                    acc, lineage = '', ''

                    if gtdb_taxonomy_format=='gtdb_taxonomy':
                        try:
                            acc, lineage = line.split('\t')
                        except:
                            logger.fatal( f"Incorrect GTDB taxonomy .tsv format: {gtdb_taxonomy_file}" )
                            _die( f"[ERROR] 2 columns are required for GTDB taxonomy .tsv format: {gtdb_taxonomy_file}" )
                        lineage = f'{lineage};x__{acc}'
                    elif gtdb_taxonomy_format=='gtdb_metadata':
                        temp = line.split('\t')
                        if len(temp)!=110:
                            logger.fatal( f"Incorrect GTDB metadata .tsv format: {gtdb_taxonomy_file}" )
                            _die( f"[ERROR] 110 columns are required for GTDB metadata .tsv format: {gtdb_taxonomy_file}" )
                        acc = temp[0]
                        # col_17: 'gtdb_taxonomy'; col_63: 'ncbi_organism_name'
                        lineage = f'{temp[16]};x__{temp[62]}'
                    else:
                        logger.fatal( f"Incorrect format: {gtdb_taxonomy_format}: {gtdb_taxonomy_file}" )
                        _die( f"[ERROR] Incorrect format: {gtdb_taxonomy_format}: {gtdb_taxonomy_file}" )

                    temp = lineage.split(';')
                    p_name = ''
                    rank = ''
                    name = ''

                    for i in range(1, len(temp)+1):
                        # this taxa
                        import re
                        re_taxa = re.compile("^([^_]+)__(.*)$")
                        rank_abbr, name = re_taxa.match(temp[-i]).groups()
                        if i==1 and rank_abbr=='x':
                            name = f'{name} ({acc})'
                        # for na taxon (no_{rank_abbr}_rank)
                        if name=="": name = p_name

                        # paranet taxa
                        try:
                            p_rank_abbr, p_name = re_taxa.match(temp[-(i+1)]).groups()
                            if p_name=="":
                                p_name = f'{name} - no_{p_rank_abbr}_rank'
                        except:
                            # for the *first* taxa in lineage line (usually superkingdom), assign parant taxid to 1 (root)
                            p_name = '1'
                            if not '1' in taxRanks: taxRanks['1'] = 'root'
                            if not '1' in taxNames: taxNames['1'] = 'root'

                        if rank_abbr=='d':
                            rank = 'superkingdom'
                        elif rank_abbr=='x':
                            rank = 'strain'
                        if rank_abbr in abbr_to_major_level:
                            rank = abbr_to_major_level[rank_abbr]
                        else:
                            rank = rank_abbr
                            
                        tid = acc.split('.')[0] if i==1 and rank=='strain' else name
                        taxParents[tid] = p_name
                        taxRanks[tid] = rank
                        taxNames[tid] = name
                        if p_name in taxNumChilds:
                            taxNumChilds[p_name] += 1
                        else:
                            taxNumChilds[p_name] = 1
                f.close()
                logger.info( f"Done parsing custom taxonomy file." )
        except IOError:
            _die( "Failed to open custom taxonomy file: %s." % gtdb_taxonomy_file )

    logger.info( f"Done parsing taxonomy files (total {len(taxParents)} taxa loaded)" )

class _autoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

def _die( msg ):
    sys.exit(msg)

def _checkTaxonomy(taxID):
    if not len(taxParents):
        logger.fatal("Taxonomy not loaded. \"loadTaxonomy()\" must be called first.")
        _die("Taxonomy not loaded. \"loadTaxonomy()\" must be called first.")

    if taxID:
        # taxID must be in string type
        taxID = str(taxID)

        # convert to merged taxID first if needs
        if taxID in taxMerged:
            logger.info( f"Merged taxID found: {taxID} -> {taxMerged[taxID]}." )
            taxID = taxMerged[taxID]

        if (taxID in taxNames) and (taxID in taxParents):
            return taxID
        else:
            return "unknown"
