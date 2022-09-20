#!/usr/bin/env perl



ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/dead_nucl.accession2taxid.gz
ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz
ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz

zcat ../taxonomy/accession2taxid/*nucl*.gz | cut -f 2,3 | sed 's/\.[[:digit:]]*//' | LC_ALL=C sort > accession2taxid.nucl.tsv
ln -s accession2taxid.nucl.tsv accession2taxid.tsv

