#!/usr/bin/env python

#
# This script is use to convert EDGE taxonomy list files to OTU formats
# 

import logging
import detaxa.taxonomy as t
import sys

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(module)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M',
)

def cli():
    t.loadTaxonomy()

    infile = sys.argv[1]
    tol_class_reads = 0
    lineage_dict = {}

    with open(infile, 'r') as fh:
        for line in fh:
            line = line.strip('\n')
            (lvl, taxa, rollup, assigned, taxid) = line.split('\t')

            # skip headers
            if lvl=='LEVEL':
                continue
            elif lvl=='unclassified':
                logging.info(f'Total number of unclassified reads: {rollup}')
                continue
            elif lvl=='root':
                tol_class_reads = int(rollup)
                logging.info(f'Total number of classified reads: {tol_class_reads}')
                continue

            # To match OTU talbes, we will only process taxa that have reads 
            # being assigned to them directly. Not all taxa are at the major 
            # ranks. We will count these reads to their least major taxa.
            if int(assigned)>0:
                lineage = t.taxid2lineage(taxid, sep=';', space2underscore=True)

                # skipping taxa not found in the taxonomy tree
                if lineage=="":
                    logging.info(f'unmapped taxa: {line}')
                    continue

                if lineage in lineage_dict:
                    lineage_dict[lineage] += int(assigned)
                else:
                    lineage_dict[lineage] = int(assigned)
            else:
                continue
        
        for lineage in sorted(lineage_dict):
            print(f'{lineage}\t{lineage_dict[lineage]}\t{lineage_dict[lineage]/tol_class_reads}')

        fh.close()

if __name__ == '__main__':
    cli()
