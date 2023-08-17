# DeTaxa

The goal of this python package is to offer a flexible, adaptable, and well-defined solution for taxonomy and accession number lookup. DeTaxa empowers users to seamlessly integrate their self-defined taxonomies into existing taxonomy systems. This library readily accommodates a range of widely used taxonomic systems, including NCBI taxonomy, EBI MGnify lineage, and GTDB taxonomy. By utilizing DeTaxa, users have the ability to create and import their customized taxonomies complete with lineages into their preferred taxonomic system. DeTaxa maintains compatibility with taxonomy files produced by [Krona](https://github.com/marbl/Krona) (Ondov et. al., 2011).

## Installation

Use python setup-tool or pip to install this package:
```
python setup.py install
```
or
```
pip install .
```

(Optional) You can run `detaxa update` to download current taxanomy file from NCBI.

## Usage

Use as a python module:

```python
#import taxonomy as module
import detaxa.taxonomy as t

#load taxonomy info
t.loadTaxonomy()

#convert taxid to name
name = t.taxid2name(tid)
```

or, run as a standalone converter:

```sh
$ detaxa query -i 2697049
```

## Acknowledgement
Part of the codes are inspired and ported from Krona taxonomy tool written in Perl.