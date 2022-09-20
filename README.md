# Taxonomy

The goal of the project is to provide a Python-based tool that is able to retrieve taxonomic information / lineage given a taxonomy ID or NCBI accession number. The tool supports NCBI taxonomy dump files and user defined taxonomies.

## Overview of functions

Loading taxonomy information:

- loadTaxonomy( dbpath=taxonomyDir )
- loadRefSeqCatelog( refseq_catelog_file, seq_type="nc" )

Convert accession number:

- acc2taxid( acc )

Convert taxonomy id:

- taxid2rank( taxID, guess_strain=True )
- taxid2name( taxID )
- taxid2depth( taxID )
- taxid2type( taxID )
- taxid2parent( taxID )
- taxid2nameOnRank( taxID, r )
- taxid2taxidOnRank( taxID, r )
- taxidIsLeaf( taxID )
- taxid2fullLineage( taxID )
- taxid2fullLinkDict( taxID )
- taxid2nearestMajorTaxid( taxID )
- taxid2lineage( tid, print_all_rank=1, print_strain=0)
- taxid2lineageDICT( tid, print_all_rank=1, print_strain=0)

## Taxonomy information files

The `taxonomy.py` will take 4 types of NCBI taxonomy information:

1. NCBI taxonomy dump file (compressed):
ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

2. NCBI taxonomy dump files (decompressed)
decompressed from taxdump.tar.gz
- nodes.dmp
- names.dmp
- merged.dmp (optional; for merged taxonomic nodes)

3. Pre-processed taxonomy information tsv files:
- taxonomy.tsv
- taxonomy.custom.tsv (optional; for custom taxonomies)
- taxonomy.merged.tsv (optional; for merged taxonomic nodes)

4. (optional) Mapping table of accession number to taxid in TSV format
- accession2taxid.tsv (optional; for function `acc2taxid()`)

## Usage

As a standalone converter:

```sh
# Pre-processed taxonomy tsv files in ./database directory
./taxonomy.py

# Decompressed NCBI taxonomy dump files in ./taxonomy directory
./taxonomy.py taxonomy
```

or can be used as a python module:

```python
#import taxonomy as module
import taxonomy as t

#load taxonomy info
t.loadTaxonomy( dbPath )

#convert taxid to name
name = taxid2name(tid)
```

## Pre-processing taxonomy information

Retrieve taxonomy information can be downloaded from NCBI FTP site. In fact, `taxdump.tar.gz` is all we need if not converting accession# to taxid.

```sh
# create database directory
mkdir -p taxonomy_db
cd taxonomy_db
rsync -auvzh --delete rsync://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz taxonomy/
rsync -auvzh --delete rsync://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid taxonomy/

tar -xzf taxonomy/taxdump.tar.gz -C taxonomy/
cp taxonomy/names.dmp .
cp taxonomy/nodes.dmp .
cp taxonomy/merged.dmp .

cd taxonomy/accession2taxid/
gzip -dc *nucl*.gz dead_wgs.accession2taxid.gz | cut -f 1,3 | LC_ALL=C sort -T . > ../../accession2taxid.nucl.tsv
ln -s accession2taxid.nucl.tsv accession2taxid.tsv
```

#### Generating taxonomy.tsv

You can use `extractTaxonomy.pl` to generate taoxnomy information file `taxonomy.tsv`. The script will compile a merged tsv file from NCBI taxonomy `nodes.dmp` and `names.dmp`.

```sh
tar -xzf ../taxonomy/taxdump.tar.gz -C ../taxonomy
../scripts/extractTaxonomy.pl ../taxonomy > taxonomy.tsv
```
The format of taxonomy.tsv (and taxonomy.custom.tsv) :

| Column | Description     | 
|--------|-----------------| 
| 1      | Taxid           | 
| 2      | depth           | 
| 3      | Parent taxid    | 
| 4      | Rank            | 
| 5      | Scientific name | 

#### Generating taxonomy.custom.tsv

`taxonomy.custom.tsv` uses the same format as `taxonomy.tsv` to expand taxonomic nodes. It's optional but recommended to add some pesudo-nodes for viral taxonomy as below:

```sh
bind "set disable-completion on"
/bin/cat <<EOM >taxonomy.custom.tsv
131567	1	1	subroot	cellular organisms
35237	2	10239	phylum	dsDNA viruses, no RNA stage
35325	2	10239	phylum	dsRNA viruses
12333	2	10239	phylum	unclassified phages
12429	2	10239	phylum	unclassified viruses
12877	2	10239	phylum	Satellites
29258	2	10239	phylum	ssDNA viruses
35268	2	10239	phylum	Retro-transcribing viruses
186616	2	10239	phylum	environmental samples
439488	2	10239	phylum	ssRNA viruses
451344	2	10239	phylum	unclassified archaeal viruses
552364	2	10239	phylum	unclassified virophages
686617	2	10239	phylum	unassigned viruses
1425366	2	10239	phylum	Virus-associated RNAs
1714266	2	10239	phylum	Virus families not assigned to an order
EOM
bind "set disable-completion off"
```

#### Generating taxonomy.merged.tsv

```sh
perl -pe 's/\s+\|//g' ../taxonomy/merged.dmp > taxonomy.merged.tsv
```

#### Generating `accession2taxid.tsv`

```sh
zcat ../taxonomy/accession2taxid/*nucl*.gz | cut -f 2,3 | sed 's/\.[[:digit:]]*//' | LC_ALL=C sort > accession2taxid.nucl.tsv
ln -s accession2taxid.nucl.tsv accession2taxid.tsv
```
