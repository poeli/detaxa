#!/usr/bin/env python3

# Po-E (Paul) Li
# B-11, Los Alamos National Lab
# Date: 01/07/2020

import sys
import io
import os.path
import json
import gzip
import tarfile
import logging
import csv

class Taxonomy(object):
    def __init__(
            self,
            dbpath="taxonomy_db",
            cus_taxonomy_file=None,
            auto_download=True,
            force_download=False,
            major_lvl_abbr={
                'superkingdom':'sk',
                'kingdom':'k',
                'phylum':'p',
                'class':'c',
                'order':'o',
                'family':'f',
                'genus':'g',
                'species':'s'
            }):
        self.dbpath = dbpath
        self.cus_taxonomy_file = cus_taxonomy_file
        self.auto_download = auto_download
        self.major_lvl_abbr = major_lvl_abbr
        self.abbr_major_lvl = dict([[v,k] for k,v in self.major_lvl_abbr.items()])

        self.taxDepths      = {}
        self.taxParents     = {}
        self.taxRanks       = {}
        self.taxNames       = {}
        self.taxMerged      = {}
        self.taxNumChilds   = {}
        self.accTid         = {}
        self.tidLineage     = {}
        self.tidLineageDict = {}

    def loadTaxonomy(self):
        logging.debug( f"Taxonomy directory: {self.dbpath}" )

        #NCBI ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
        taxdump_tgz_file = self.dbpath+"/taxdump.tar.gz"

        #raw taxonomy dmp files from NCBI
        names_dmp_file = self.dbpath+"/names.dmp"
        nodes_dmp_file = self.dbpath+"/nodes.dmp"
        merged_dmp_file = self.dbpath+"/merged.dmp"

        #parsed taxonomy tsv file
        taxonomy_file = self.dbpath+"/taxonomy.tsv"
        merged_taxonomy_file = self.dbpath+"/taxonomy.merged.tsv"

        #custom taxonomy file
        if not self.cus_taxonomy_file:
            self.cus_taxonomy_file = self.dbpath+"/taxonomy.custom.tsv"

        # checking if taxonomy files downloaded
        if self.force_download or (not os.path.isfile( taxdump_tgz_file ) \
          and not os.path.isfile( merged_taxonomy_file ) \
          and not (os.path.isfile( names_dmp_file ) and os.path.isfile( nodes_dmp_file )) \
          and not os.path.isfile( self.cus_taxonomy_file)):
            import requests
            logging.debug( "Local taxonomy files not found." )
            if self.auto_download or self.force_download:
                url = 'http://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
                logging.debug( f"Auto downloading taxanomy from {url}..." )
                # download taxonomy file if auto_download enabled
                r = requests.get(url)
                if not os.path.exists( self.dbpath):
                    os.makedirs( self.dbpath )
                with open(f'{self.dbpath}/taxdump.tar.gz', 'wb') as f:
                    f.write(r.content)
                if os.path.getsize( f'{self.dbpath}/taxdump.tar.gz'):
                    logging.debug( f'Saved to {self.dbpath}/taxdump.tar.gz.' )
                else:
                    self._die( "[ERROR] Failed to download taxonomy files.")
            else:
                logging.debug( "Auto-download is off.")
                self._die( "[ERROR] No available taxonomy files.")

        # unpacking taxdump.tar.gz
        if os.path.isfile(taxdump_tgz_file):
            try:
                logging.info( f"Unpacking taxonomy file: {taxdump_tgz_file}")
                tar = tarfile.open(taxdump_tgz_file, "r:gz")
                tar.extractall(self.dbpath)
                tar.close()
            except IOError:
                self._die( "Failed to load taxonomy from %s\n"%taxdump_tgz_file )
            # delete taxdump.tar.gz
            os.remove(taxdump_tgz_file)
        
        if os.path.isfile(taxonomy_file):
            logging.debug( f"Open taxonomy file: %s\n"% taxonomy_file )
            try:
                with open(taxonomy_file) as f:
                    for line in f:
                        tid, depth, parent, rank, name = line.rstrip('\r\n').split('\t')
                        self.taxParents[tid] = parent
                        self.taxDepths[tid] = depth
                        self.taxRanks[tid] = rank
                        self.taxNames[tid] = name
                        if parent in self.taxNumChilds:
                            self.taxNumChilds[parent] += 1
                        else:
                            self.taxNumChilds[parent] = 1
                    f.close()
                    logging.debug( f"Done parsing taxonomy file.")

                #try to load merged taxids
                if os.path.isfile( merged_taxonomy_file):
                    logging.debug( f"Open merged taxonomy node file: %s\n"% merged_taxonomy_file )
                    with open(merged_taxonomy_file) as f:
                        for line in f:
                            line = line.rstrip('\r\n')
                            if not line: continue
                            mtid, tid = line.split('\t')
                            self.taxMerged[mtid] = tid
                        f.close()
                        logging.debug( f"Done parsing merged taxonomy file.")
            except IOError:
                self._die( "Failed to open taxonomy file: %s.\n" % taxonomy_file )
        # loading from *.dmp files
        elif os.path.isfile( names_dmp_file ):
            try:
                # read name from names.dmp
                logging.debug( f"Open taxonomy name file: {names_dmp_file}\n"  )
                with open(names_dmp_file) as f:
                    for line in f:
                        tid, name, tmp, nametype = line.rstrip('\r\n').split('\t|\t')
                        if not nametype.startswith("scientific name"):
                            continue
                        self.taxNames[tid] = name
                    f.close()
                    logging.debug( f"Done parsing taxonomy name file.")	

                # read taxonomy info from nodes.dmp
                logging.debug( f"Open taxonomy node file: {nodes_dmp_file}")
                with open(nodes_dmp_file) as f:
                    for line in f:
                        fields = line.rstrip('\r\n').split('\t|\t')
                        tid = fields[0]
                        parent = fields[1]
                        self.taxParents[tid] = parent
                        self.taxDepths[tid] = self.taxDepths[parent]+1 if parent in self.taxDepths else 0 # could have potiential bug if child node is parsed before parent node.
                        self.taxRanks[tid] = fields[2]
                        if parent in self.taxNumChilds:
                            self.taxNumChilds[parent] += 1
                        else:
                            self.taxNumChilds[parent] = 1
                    f.close()
                    logging.debug( f"Done parsing taxonomy node file.")

                if os.path.isfile( merged_dmp_file):
                    logging.debug( f"Open merged taxonomy node file: {merged_dmp_file}")
                    with open(merged_dmp_file) as f:
                        for line in f:
                            line = line.rstrip('\r\n')
                            if not line: continue
                            fields = line.replace("\t","").split('|')
                            mtid = fields[0]
                            tid = fields[1]
                            self.taxMerged[mtid] = tid
                        f.close()
                        logging.debug( f"Done parsing merged taxonomy file.")
            except IOError:
                self._die( "Failed to open taxonomy files (taxonomy.tsv, nodes.dmp and names.dmp).")

        # try to load custom taxonomy from taxonomy.custom.tsv
        if os.path.isfile( cus_taxonomy_file):
            logging.debug( f"Open custom taxonomy node file: %s\n"% cus_taxonomy_file)
            try:
                with open(cus_taxonomy_file) as f:
                    for line in f:
                        line = line.rstrip('\r\n')
                        if not line: continue
                        tid, depth, parent, rank, name = line.split('\t')
                        self.taxParents[tid] = parent
                        self.taxDepths[tid] = depth
                        self.taxRanks[tid] = rank
                        self.taxNames[tid] = name
                        if parent in self.taxNumChilds:
                            self.taxNumChilds[parent] += 1
                        else:
                            self.taxNumChilds[parent] = 1
                    f.close()
                    logging.debug( f"Done parsing custom taxonomy file.")
            except IOError:
                self._die( "Failed to open custom taxonomy file: %s.\n" % cus_taxonomy_file )

        logging.debug( f"Done parsing taxonomy files (%d taxons loaded)\n" % len(self.taxParents) )

    def _die(self, msg):
        sys.exit(msg)

    def _checkTaxonomy(self, taxID="", acc=""):
        if not len(self.taxParents):
            self._die("Taxonomy not loaded. \"loadTaxonomy()\" must be called first.\n")

        if taxID:
            if taxID in self.taxMerged:
                taxID = self.taxMerged[taxID]

        if taxID in self.taxNames and taxID in self.taxParents and taxID in self.taxParents:
            return taxID
        else:
            return "unknown"

    def to_json_tree(self):
        import json

        links = [(self.taxParents[c], c) for c in self.taxParents]

        parents, children = zip(*links)
        root_nodes = {x for x in parents if x not in children}

        for node in root_nodes:
            links.append(('Root', node))

        def get_nodes(node):
            d = {}
            d['value'] = self.taxNames[node].lower()
            d['label'] = self.taxNames[node]
            d['id'] = node
            children = get_children(node)
            if children:
                d['children'] = [getself._nodes(child) for child in children]
            return d

        def get_children(node):
            return [x[1] for x in links if x[0] == node]

        tree = getself._nodes('Root')
        print(json.dumps(tree, indent=4))

    def taxidStatus(self, taxID):
        if taxID in self.taxMerged:
            return self.taxMerged[taxID]

        if taxID in self.taxNames and taxID in self.taxNames and taxID in self.taxRanks:
            if '.' in taxID:
                return "valid custom"
            return "valid"
        else:
            return "invalid"

    def acc2taxid(self, acc):
        self._checkTaxonomy()

        accession2taxid_file = f"{self.dbpath}/accession2taxid.tsv"
        #remove version number#
        acc = acc.split('.')[0]

        logging.debug( f"acc2taxid from file: {accession2taxid_file}")

        if not acc in self.accTid:
            with open( accession2taxid_file ) as f:
                f.seek(0, 2)
                start = 0
                end = f.tell()
                accCur = ""
                
                logging.debug( f"acc2taxid from file: {accession2taxid_file}")
                
                while(acc != accCur and start < end):
                    
                    posNew = (end+start)/2
                    
                    f.seek( posNew )
            
                    if posNew != start: f.readline()

                    line = f.readline()	
                    
                    logging.debug( f"start: %15d, posNew: %15d, end: %15d, line: %s" % (start, posNew, end, line) )
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
                    self.accTid[acc] = tid.strip()
                else:
                    self.accTid[acc] = ""

        tid = self._checkTaxonomy(self.accTid[acc])

        return tid

    def taxid2rank(self, taxID, guess_strain=True):
        taxID = self._checkTaxonomy( taxID )
        if taxID == "unknown": return "unknown"

        if taxID == '1':
            return "root"

        if self.taxRanks[taxID] == "no rank" and guess_strain:
            # a leaf taxonomy is a strain
            if taxidIsLeaf(taxID):
                return "strain"
            # if not
            else:
                nmtid = taxid2nearestMajorTaxid(taxID)
                nmrank = self._getTaxRank(nmtid)
                if nmrank == "species":
                    return "species - others"
                else:
                    return "others"
        
        return self.taxRanks[taxID]

    def taxid2name(self, taxID):
        taxID = self._checkTaxonomy( taxID )
        if taxID == "unknown":
            return "unknown"
        else:
            return self._getTaxName(taxID)

    def taxid2depth(self, taxID):
        taxID = self._checkTaxonomy( taxID )
        if taxID == "unknown":
            return "unknown"
        else:
            return self._getTaxDepth(taxID)

    def taxid2type(self, taxID):
        taxID = self._checkTaxonomy( taxID )
        if taxID == "unknown": return "unknown"

        origID = taxID
        lastID = taxID
        taxID = self.taxParents[taxID]

        while taxID != '1' and self.taxRanks[taxID] != 'species':
            lastID = taxID
            taxID = self.taxParents[taxID]

        if self.taxRanks[taxID] != 'species':
            taxID = 0
        else:
            taxID = lastID
            if taxID == origID: taxID = 0

        return taxID

    def taxid2parent(self, taxID):
        taxID = self._checkTaxonomy( taxID )
        if taxID == "unknown": return "unknown"

        taxID = self.taxParents[taxID]
        while taxID != '1' and self.taxRanks[taxID] == 'no rank':
            taxID = self.taxParents[taxID]

        return taxID

    def taxid2nameOnRank(self, taxID, target_rank=None):
        taxID = self._checkTaxonomy( taxID )
        if taxID == "unknown": return "unknown"

        if taxID == 1: return "root"
        if target_rank == "root": return "root"

        rank = self._getTaxRank(taxID)
        name = self._getTaxName(taxID)

        if target_rank == "strain" and taxidIsLeaf(taxID):
            return name

        while taxID:
            if rank.upper() == target_rank.upper(): return name
            if name == 'root': break
            taxID = self._getTaxParent(taxID)
            rank = self._getTaxRank(taxID)
            name = self._getTaxName(taxID)

        return ""

    def taxid2taxidOnRank(self, taxID, target_rank=None):
        taxID = self._checkTaxonomy( taxID )
        if taxID == "unknown": return "unknown"

        rank = self._getTaxRank(taxID)
        name = self._getTaxName(taxID)

        if target_rank == rank or ( target_rank == 'strain' and rank == 'no rank'): return taxID
        if target_rank == "root": return 1

        while taxID:
            if rank.upper() == target_rank.upper(): return taxID
            if name == 'root': break

            taxID = self._getTaxParent(taxID)
            rank = self._getTaxRank(taxID)
            name = self._getTaxName(taxID)

        return ""

    def taxidIsLeaf(self, taxID):
        taxID = self._checkTaxonomy( taxID )
        if taxID == "unknown": return False
        if not taxID in self.taxNumChilds:
            return True
        else:
            return False

    def taxid2fullLineage(self, taxID):
        taxID = self._checkTaxonomy( taxID )
        if taxID == "unknown": return "unknown"
        fullLineage = ""

        while taxID != '1':
            rank = self._getTaxRank(taxID)
            name = self._getTaxName(taxID)
            if not name: break
            fullLineage += "%s|%s|%s|"%(rank,taxID,name)
            taxID = self.taxParents[taxID]

        return fullLineage

    def taxid2fullLinkDict(self, taxID):
        taxID = self._checkTaxonomy( taxID )
        if taxID == "unknown": return "unknown"
        link = {}

        while taxID != '1':
            name = self._getTaxName(taxID)
            if not name: break

            parID = self.taxParents[taxID]
            link[parID] = taxID
            taxID = parID

        return link

    def taxid2nearestMajorTaxid(self, taxID):
        taxID = self._checkTaxonomy( taxID )
        if taxID == "unknown": return "unknown"
        ptid = self._getTaxParent( taxID )
        while ptid != '1':
            tmp = taxid2rank( ptid )
            if tmp in self.major_lvl_abbr:
                return ptid
            else:
                ptid = self._getTaxParent( ptid )

        return "1"

    def taxid2lineage(self, tid, print_all_rank=True, print_strain=False, replace_space2underscore=True, output_type="auto"):
        return self._taxid2lineage( tid, print_all_rank, print_strain, replace_space2underscore, output_type)

    def taxid2lineageDICT(self, tid, print_all_rank=True, print_strain=False, replace_space2underscore=False, output_type="DICT"):
        return self._taxid2lineage( tid, print_all_rank, print_strain, replace_space2underscore, output_type )

    def taxid2lineageTEXT(self, tid, print_all_rank=True, print_strain=False, replace_space2underscore=True, output_type="DICT"):
        lineage = self._taxid2lineage( tid, print_all_rank, print_strain, replace_space2underscore, output_type)
        texts = []
        for rank in self.major_lvl_abbr:
            if rank in lineage:
                texts.append( f"{self.major_lvl_abbr[rank]}__{lineage[rank]['name']}" ) 
        
        return ";".join(texts).replace(" ","_")

    def _taxid2lineage(self, tid, print_all_rank, print_strain, replace_space2underscore, output_type):
        taxID = self._checkTaxonomy( tid )
        if taxID == "unknown": return "unknown"

        if output_type == "DICT":
            if taxID in self.tidLineageDict: return self.tidLineageDict[taxID]
        else:
            if taxID in self.tidLineage: return self.tidLineage[taxID]

        info = self._autoVivification()
        lineage = []

        level = {
            'k' : '',
            'p' : '',
            'c' : '',
            'o' : '',
            'f' : '',
            'g' : '',
            's' : ''
        }

        rank = taxid2rank(taxID)
        orig_rank = rank
        name = self._getTaxName(taxID)
        str_name = name
        if replace_space2underscore: str_name.replace(" ", "_")

        while taxID:
            if rank in self.major_lvl_abbr:
                if replace_space2underscore: name.replace(" ", "_")
                level[self.major_lvl_abbr[rank]] = name

                #for output JSON
                info[rank]["name"] = name
                info[rank]["taxid"] = taxID

            taxID = self._getTaxParent(taxID)
            rank = self._getTaxRank(taxID)
            name = self._getTaxName(taxID)

            if name == 'root': break

        # try to get the closest "no_rank" taxa to "type" representing subtype/group (mainly for virus)
        typeTID = taxid2type(tid)
        if typeTID:
            info["type"]["name"]  = self._getTaxName(typeTID)
            info["type"]["taxid"] = typeTID

        last = str_name

        ranks = ['s','g','f','o','c','p','k']
        idx = 0
        
        # input taxid is a major rank
        if orig_rank in self.major_lvl_abbr:
            idx = ranks.index( self.major_lvl_abbr[orig_rank] )
        # if not, find the next major rank
        else:
            nmtid = taxid2nearestMajorTaxid( tid )
            nmrank = taxid2rank( nmtid )
            if nmrank == "root":
                idx = 7
            else:
                idx = ranks.index( self.major_lvl_abbr[nmrank] )

        for lvl in ranks[idx:]:
            if print_all_rank == 0:
                if not level[lvl]: continue

            if not level[lvl]:
                level[lvl] = "%s - no_%s_rank"%(last,lvl)
                info[self.abbr_major_lvl[lvl]]["name"]  = "%s - no_%s_rank"%(last,lvl)
                info[self.abbr_major_lvl[lvl]]["taxid"] = 0

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
            self.tidLineageDict[tid] = info
            return info
        else:
            self.tidLineage[tid] = "|".join(lineage)
            return "|".join(lineage)

    def _getTaxDepth(self, taxID):
        return self.taxDepths[taxID]

    def _getTaxName(self, taxID):
        return self.taxNames[taxID]

    def _getTaxParent(self, taxID):
        return self.taxParents[taxID]

    def _getTaxRank(self, taxID):
        return self.taxRanks[taxID]

    def lcaself._taxid(self, taxids):
        """ lca_taxid
        Return lowest common ancestor (LCA) taxid of input taxids
        """
        ranks = ['strain','species','genus','family','order','class','phylum','superkingdom']

        merged_dict = self._autoVivification()
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

class _autoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.self.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value
