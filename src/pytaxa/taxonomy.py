#!/usr/bin/env python
__version__="0.5.0"

# Po-E (Paul) Li
# B-11, Los Alamos National Lab
# Date: 05/15/2016

import sys
import os
import tarfile
import requests
import logging

####################
# Global variables #
####################

# Set default path of taxonomy_db/:
# The default `taxonomy_db/` path is your the location
lib_path = os.path.dirname(os.path.realpath(__file__))
taxonomy_dir = lib_path + "/taxonomy_db"

if os.path.isdir( "./taxonomy_db" ):
    taxonomy_dir = "./taxonomy_db"
elif os.path.isdir( os.getcwd()+"/taxonomy_db" ):
    taxonomy_dir = os.getcwd()+"/taxonomy_db"

taxDepths      = {}
taxParents     = {}
taxRanks       = {}
taxNames       = {}
taxMerged      = {}
taxNumChilds   = {}
accTid         = {}
tidLineage     = {}
tidLineageDict = {}

major_level_to_abbr = {
    'superkingdom' : 'sk',
    'phylum'       : 'p',
    'class'        : 'c',
    'order'        : 'o',
    'family'       : 'f',
    'genus'        : 'g',
    'species'      : 's',
    'strain'       : 'n',
}
abbr_to_major_level = {
    'sk'           : 'superkingdom',
    'p'            : 'phylum',
    'c'            : 'class',
    'o'            : 'order',
    'f'            : 'family',
    'g'            : 'genus',
    's'            : 'species',
    'n'            : 'strain',
}

####################
#      Methods     #
####################

def taxidStatus( taxID ):
    if taxID in taxMerged:
        return taxMerged[taxID]

    if taxID in taxNames and taxID in taxNames and taxID in taxRanks:
        if '.' in taxID:
            return "valid custom"
        return "valid"
    else:
        return "invalid"

def acc2taxid( acc ):
    _checkTaxonomy()

    accession2taxid_file = f"{taxonomy_dir}/accession2taxid.tsv"
    #remove version number#
    acc = acc.split('.')[0]

    logging.info( f"acc2taxid from file: {accession2taxid_file}" )

    if not acc in accTid:
        with open( accession2taxid_file ) as f:
            f.seek(0, 2)
            start = 0
            end = f.tell()
            accCur = ""
            
            logging.info( f"acc2taxid from file: {accession2taxid_file}" )
            
            while( acc != accCur and start < end ):
                
                posNew = (end+start)/2
                
                f.seek( posNew )
        
                if posNew != start: f.readline()

                line = f.readline()    
                
                logging.info( "start: %15d, posNew: %15d, end: %15d, line: %s" % (start, posNew, end, line) )
                if line :
                    (accNew, tid) = line.split('\t')
                else:
                    break

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

    tid = _checkTaxonomy(accTid[acc])

    return tid

def taxid2rank( taxID, guess_strain=True ):
    taxID = _checkTaxonomy( taxID )
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

def taxid2name( taxID ):
    taxID = _checkTaxonomy( taxID )
    if taxID == "unknown":
        return "unknown"
    else:
        return _getTaxName(taxID)

def taxid2depth( taxID ):
    taxID = _checkTaxonomy( taxID )
    if taxID == "unknown":
        return "unknown"
    else:
        return _getTaxDepth(taxID)

def taxid2type( taxID ):
    taxID = _checkTaxonomy( taxID )
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

def taxid2parent( taxID ):
    taxID = _checkTaxonomy( taxID )
    if taxID == "unknown": return "unknown"

    taxID = taxParents[taxID]
    while taxID != '1' and taxRanks[taxID] == 'no rank':
        taxID = taxParents[taxID]

    return taxID

def taxid2nameOnRank( taxID, target_rank=None ):
    taxID = _checkTaxonomy( taxID )
    if taxID == "unknown": return "unknown"

    if taxID == 1: return "root"
    if target_rank == "root": return "root"

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

    return ""

def taxid2taxidOnRank( taxID, target_rank=None ):
    taxID = _checkTaxonomy( taxID )
    if taxID == "unknown": return "unknown"

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

    return ""

def taxidIsLeaf( taxID ):
    taxID = _checkTaxonomy( taxID )
    if taxID == "unknown": return False
    if not taxID in taxNumChilds:
        return True
    else:
        return False

def taxid2fullLineage( taxID ):
    taxID = _checkTaxonomy( taxID )
    if taxID == "unknown": return "unknown"
    fullLineage = ""

    while taxID != '1':
        rank = _getTaxRank(taxID)
        name = _getTaxName(taxID)
        if not name: break
        fullLineage += "%s|%s|%s|"%(rank,taxID,name)
        taxID = taxParents[taxID]

    return fullLineage

def taxid2fullLinkDict( taxID ):
    taxID = _checkTaxonomy( taxID )
    if taxID == "unknown": return "unknown"
    link = {}

    while taxID != '1':
        name = _getTaxName(taxID)
        if not name: break

        parID = taxParents[taxID]
        link[parID] = taxID
        taxID = parID

    return link

def taxid2nearestMajorTaxid( taxID ):
    taxID = _checkTaxonomy( taxID )
    if taxID == "unknown": return "unknown"
    ptid = _getTaxParent( taxID )
    while ptid != '1':
        tmp = taxid2rank( ptid )
        if tmp in major_level_to_abbr:
            return ptid
        else:
            ptid = _getTaxParent( ptid )

    return "1"

def taxid2lineage( tid, print_all_rank=True, print_strain=False, replace_space2underscore=True, output_type="auto"):
    return _taxid2lineage( tid, print_all_rank, print_strain, replace_space2underscore, output_type)

def taxid2lineageDICT( tid, print_all_rank=True, print_strain=False, replace_space2underscore=False, output_type="DICT" ):
    return _taxid2lineage( tid, print_all_rank, print_strain, replace_space2underscore, output_type )

def taxid2lineageTEXT( tid, print_all_rank=True, print_strain=False, replace_space2underscore=True, output_type="DICT"):
    lineage = _taxid2lineage( tid, print_all_rank, print_strain, replace_space2underscore, output_type)
    texts = []
    for rank in major_level_to_abbr:
        if rank in lineage:
            texts.append( f"{major_level_to_abbr[rank]}__{lineage[rank]['name']}" ) 
    
    return ";".join(texts).replace(" ","_")

def _taxid2lineage(tid, print_all_rank, print_strain, replace_space2underscore, output_type):
    taxID = _checkTaxonomy( tid )
    if taxID == "unknown": return "unknown"

    if output_type == "DICT":
        if taxID in tidLineageDict: return tidLineageDict[taxID]
    else:
        if taxID in tidLineage: return tidLineage[taxID]

    info = _autoVivification()
    lineage = []

    level = {
        'sk' : '',
        'p' : '',
        'c' : '',
        'o' : '',
        'f' : '',
        'g' : '',
        's' : '',
        'n' : ''
    }

    rank = taxid2rank(taxID)
    orig_rank = rank
    name = _getTaxName(taxID)
    str_name = name
    if replace_space2underscore: str_name.replace(" ", "_")

    while taxID:
        if rank in major_level_to_abbr:
            if replace_space2underscore: name.replace(" ", "_")
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

    ranks = ['n','s','g','f','o','c','p','sk']
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
        if print_all_rank == 0:
            if not level[lvl]: continue

        if not level[lvl]:
            level[lvl] = "%s - no_%s_rank"%(last,lvl)
            info[abbr_to_major_level[lvl]]["name"]  = "%s - no_%s_rank"%(last,lvl)
            info[abbr_to_major_level[lvl]]["taxid"] = 0

        last=level[lvl]
        #lineage.append( "%s__%s"%(lvl, level[lvl]) )
        lineage.append( level[lvl] )

    lineage.reverse()

    if print_strain:
        if orig_rank == "strain":
            #lineage.append( "n__%s"%(str_name) )
            lineage.append( "%s"%(str_name) )
            info["strain"]["name"]  = str_name
            info["strain"]["taxid"] = tid

    if output_type == "DICT":
        tidLineageDict[tid] = info
        return info
    else:
        tidLineage[tid] = "|".join(lineage)
        return "|".join(lineage)

def _getTaxDepth( taxID ):
    return taxDepths[taxID]

def _getTaxName( taxID ):
    return taxNames[taxID]

def _getTaxParent( taxID ):
    return taxParents[taxID]

def _getTaxRank( taxID ):
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

def loadTaxonomy( dbpath=None, 
                  cus_taxonomy_file=None, 
                  cus_taxonomy_format="tsv", 
                  auto_download=True
    ):
    global taxonomy_dir

    if dbpath:
        taxonomy_dir = dbpath

    logging.info( f"Taxonomy directory: {taxonomy_dir}" )

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

    # checking if taxonomy files downloaded
    if not os.path.isfile( taxdump_tgz_file ) \
        and not os.path.isfile( merged_taxonomy_file ) \
        and not os.path.isfile( taxonomy_file ) \
        and not (os.path.isfile( names_dmp_file ) and os.path.isfile( nodes_dmp_file )) \
        and not os.path.isfile( cus_taxonomy_file ):
        
        logging.info( f"Local taxonomy files not found." )
        if auto_download:
            url = 'http://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
            logging.info( f"Auto downloading taxanomy from {url}..." )
            # download taxonomy file if auto_download enabled
            r = requests.get(url)
            if not os.path.exists( taxonomy_dir ):
                os.makedirs( taxonomy_dir )
            
            taxdump_tgz_file = f'{taxonomy_dir}/taxdump.tar.gz'

            with open(taxdump_tgz_file, 'wb') as f:
                f.write(r.content)
            if os.path.getsize( taxdump_tgz_file ):
                logging.info( f"Saved to {taxdump_tgz_file}." )
            else:
                logging.fatal( f"Failed to download or save taxonomy files." )
                _die( "[ERROR] Failed to download or save taxonomy files." )    

            # extract
            tax_tar = tarfile.open(taxdump_tgz_file, "r:gz")
            tax_tar.extract('nodes.dmp', taxonomy_dir)
            tax_tar.extract('names.dmp', taxonomy_dir)
            tax_tar.extract('merged.dmp', taxonomy_dir)
            tax_tar.close()
            # delete taxdump_tgz_file
            os.remove(taxdump_tgz_file)

        else:
            logging.info( f"Auto-download is off." )
            logging.fatal( f"No available taxonomy files." )
            _die( "[ERROR] No available taxonomy files." )

    # try to load taxonomy from taxonomy.tsv
    if os.path.isfile( nodes_dmp_file ) and  os.path.isfile( names_dmp_file ):
        try:
            # read name from names.dmp
            logging.info( f"Open taxonomy name file: {names_dmp_file}" )
            with open(names_dmp_file) as f:
                for line in f:
                    tid, name, tmp, nametype = line.rstrip('\r\n').split('\t|\t')
                    if not nametype.startswith("scientific name"):
                        continue
                    taxNames[tid] = name
                f.close()
                logging.info( f"Done parsing taxonomy name file." )    

            # read taxonomy info from nodes.dmp
            logging.info( f"Open taxonomy node file: {nodes_dmp_file}" )
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
                logging.info( f"Done parsing taxonomy node file." )
        except IOError:
            _die( "Failed to open taxonomy files (taxonomy.tsv, nodes.dmp and names.dmp)." )
    elif os.path.isfile( taxdump_tgz_file ):
        try:
            logging.info( f"Open taxonomy file: {taxdump_tgz_file}" )
            tar = tarfile.open(taxdump_tgz_file, "r:gz")
            
            # read name from names.dmp
            logging.info( "Extract taxonomy names file: names.dmp" )
            member = tar.getmember("names.dmp")
            f = tar.extractfile(member)
            for line in f.readlines():
                tid, name, tmp, nametype = line.decode('utf8').rstrip('\r\n').split('\t|\t')
                if not nametype.startswith("scientific name"):
                    continue
                taxNames[tid] = name
            f.close()
            
            # read taxonomy info from nodes.dmp
            logging.info( "Extract taxonomy nodes file: nodes.dmp" )
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
        except IOError:
            _die( "Failed to load taxonomy from %s"%taxdump_tgz_file )
    elif os.path.isfile( taxonomy_file ):
        logging.info( "Open taxonomy file: %s"% taxonomy_file )
        try:
            with open(taxonomy_file) as f:
                for line in f:
                    tid, depth, parent, rank, name = line.rstrip('\r\n').split('\t')
                    taxParents[tid] = parent
                    taxDepths[tid] = depth
                    taxRanks[tid] = rank
                    taxNames[tid] = name
                    if parent in taxNumChilds:
                        taxNumChilds[parent] += 1
                    else:
                        taxNumChilds[parent] = 1
                f.close()
                logging.info( f"Done parsing taxonomy file." )
        except IOError:
            _die( "Failed to open taxonomy file: %s." % taxonomy_file )

    #try to load merged taxids
    if os.path.isfile( merged_taxonomy_file ):
        logging.info( "Open merged taxonomy node file: %s"% merged_taxonomy_file )
        with open(merged_taxonomy_file) as f:
            for line in f:
                line = line.rstrip('\r\n')
                if not line: continue
                mtid, tid = line.split('\t')
                taxMerged[mtid] = tid
            f.close()
            logging.info( f"Done parsing merged taxonomy file." )

    # try to load custom taxonomy from taxonomy.custom.tsv
    if os.path.isfile(cus_taxonomy_file) and (cus_taxonomy_format=='tsv'):
        logging.info( "Open custom taxonomy node file (tsv format): %s"% cus_taxonomy_file)
        try:
            with open(cus_taxonomy_file) as f:
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
                logging.info( f"Done parsing custom taxonomy file." )
        except IOError:
            _die( "Failed to open custom taxonomy file: %s." % cus_taxonomy_file )

    # try to load custom taxonomy from lineage file
    if os.path.isfile(cus_taxonomy_file) and (cus_taxonomy_format=='lineage'):
        logging.info( "Open custom taxonomy node file (lineage format): %s"% cus_taxonomy_file)
        try:
            with open(cus_taxonomy_file) as f:
                for line in f:
                    line = line.rstrip('\r\n')
                    if not line: continue
                    if line.startswith('#'): continue
                    if not line.startswith('sk__'):
                        logging.warn( f"A text line of lineage has to start with 'sk__'...skipped: {line}" )
                        continue

                    temp = line.split(';')
                    p_name = ''
                    rank = ''
                    name = ''

                    for i in range(1, len(temp)+1):
                        # this taxa
                        (rank_abbr, name) = temp[-i].split('__')
                        # for na taxon (no_{rank_abbr}_rank)
                        if name=="":
                            name = p_name

                        # paranet taxa
                        try:
                            p_name = temp[-(i+1)].split('__')[1]
                            if p_name=="":
                                p_name = f'{name} - no_{rank_abbr}_rank'
                        except:
                            # for the superkingdom rank, assign parant taxid to 1 (root)
                            p_name = '1'

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
                logging.info( f"Done parsing custom taxonomy file." )
        except IOError:
            _die( "Failed to open custom taxonomy file: %s." % cus_taxonomy_file )


    logging.info( "Done parsing taxonomy files (%d taxons loaded)" % len(taxParents) )

##########################
##  Internal functions  ##
##########################

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

def _checkTaxonomy(taxID="", acc=""):
    if not len(taxParents):
        _die("Taxonomy not loaded. \"loadTaxonomy()\" must be called first.")

    taxID = str(taxID)

    if taxID:
        if taxID in taxMerged:
            taxID = taxMerged[taxID]

    if (taxID in taxNames) and (taxID in taxParents):
        return taxID
    else:
        return "unknown"