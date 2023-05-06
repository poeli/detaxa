#! python
import os, sys, logging
import argparse as ap
from . import taxonomy as t
from . import __version__
import click

@click.group(help=f"""deTaxa taxonomy utility v{__version__}.""")
def cli():
    pass

@cli.command()
@click.argument('taxid', required=True, type=str)
@click.option('-d', '--database',
              help='path of taxonomy_db/',
              required=False,
              default=None,
              type=str)
@click.option('-c', '--custom-taxa',
              help='path of custom taxonomy file',
              required=False,
              default=None,
              type=str)
@click.option('-f', '--custom-fmt',
              help="custom taxonomy format 'tsv', 'lineage', 'gtdb_taxonomy' and 'gtdb_metadata'",
              required=False,
              default='tsv',
              type=click.Choice(['tsv', 'lineage', 'gtdb_taxonomy', 'gtdb_metadata'], case_sensitive=False)
              )
@click.option('--debug',
              help='debug mode',
              is_flag=True,
              default=False)

def taxid(taxid, database, custom_taxa, custom_fmt, debug):
    if debug:
        logging.basicConfig(
            level=logging.DEBUG,
            format='%(asctime)s [%(levelname)s] %(module)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M',
        )

    if custom_fmt.startswith('gtdb'):
        t.loadGTDBTaxonomy(cus_taxonomy_file=custom_taxa, cus_taxonomy_format=custom_fmt)
    else:
        t.loadTaxonomy( database, cus_taxonomy_file=custom_taxa, cus_taxonomy_format=custom_fmt)

    if taxid:
        print( "taxid2name( %s )                 => %s" % (taxid, t.taxid2name(taxid)) )
        print( "taxid2rank( %s )                 => %s" % (taxid, t.taxid2rank(taxid)) )
        print( "taxid2type( %s )                 => %s" % (taxid, t.taxid2type(taxid)) )
        print( "taxid2parent( %s )               => %s" % (taxid, t.taxid2parent(taxid)) )
        print( "taxidIsLeaf( %s )                => %s" % (taxid, t.taxidIsLeaf(taxid)) )
        print( "taxid2nearestMajorTaxid( %s )    => %s" % (taxid, t.taxid2nearestMajorTaxid(taxid)) )
        print( "taxid2nameOnRank( %s, 'genus')   => %s" % (taxid, t.taxid2nameOnRank(taxid, "genus")) )
        print( "taxid2taxidOnRank( %s, 'genus')  => %s" % (taxid, t.taxid2taxidOnRank(taxid, "genus")) )
        print( "taxid2nameOnRank( %s, 'phylum')  => %s" % (taxid, t.taxid2nameOnRank(taxid, "phylum")) )
        print( "taxid2taxidOnRank( %s, 'phylum') => %s" % (taxid, t.taxid2taxidOnRank(taxid, "phylum")) )
        print( "taxid2lineage( %s, sep=';' )     => %s" % (taxid, t.taxid2lineage(taxid, sep=';')) )
        print( "taxid2lineageDICT( %s )          => %s" % (taxid, t.taxid2lineageDICT(taxid)) )
        print( "taxid2fullLineage( %s )          => %s" % (taxid, t.taxid2fullLineage(taxid)) )
        print( "taxid2fullLineage( %s, sep=';' ) => %s" % (taxid, t.taxid2fullLineage(taxid, sep=';')) )
        print( "taxid2fullLinkDict( %s )         => %s" % (taxid, t.taxid2fullLinkDict(taxid)) )
    else:
        print( "No taxid found." )


@cli.command()
@click.argument('name', required=True, type=str)
@click.option('-d', '--database',
              help='path of taxonomy_db/',
              required=False,
              default=None,
              type=str)
@click.option('-c', '--custom-taxa',
              help='path of custom taxonomy file',
              required=False,
              default=None,
              type=str)
@click.option('-f', '--custom-fmt',
              help="custom taxonomy format 'tsv', 'lineage', 'gtdb_taxonomy' and 'gtdb_metadata'",
              required=False,
              default='tsv',
              type=click.Choice(['tsv', 'lineage', 'gtdb_taxonomy', 'gtdb_metadata'], case_sensitive=False)
              )
@click.option('-r', '--rank',
              help="rank of the input taxa",
              required=False,
              default=None,
              type=str
              )
@click.option('-p', '--partial',
              help="partial match of input taxa, default is False",
              is_flag=True,
              default=False
              )
@click.option('--debug',
              help='debug mode',
              is_flag=True,
              default=False)

def name2tid(name, database, custom_taxa, custom_fmt, rank, partial, debug):
    if debug:
        logging.basicConfig(
            level=logging.DEBUG,
            format='%(asctime)s [%(levelname)s] %(module)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M',
        )
    
    if custom_fmt.startswith('gtdb'):
        t.loadGTDBTaxonomy(cus_taxonomy_file=custom_taxa, cus_taxonomy_format=custom_fmt)
    else:
        t.loadTaxonomy( database, cus_taxonomy_file=custom_taxa, cus_taxonomy_format=custom_fmt)
 
    print(t.name2taxid(name, rank, partial))

@cli.command()
@click.argument('accession', required=True, type=str)
@click.option('-m', '--mapping',
              help='path of mapping table',
              required=False,
              default=None,
              type=str)
@click.option('--debug',
              help='debug mode',
              is_flag=True,
              default=False)
def acc2taxid(accession, mapping_tsv, debug):
    if debug:
        logging.basicConfig(
            level=logging.DEBUG,
            format='%(asctime)s [%(levelname)s] %(module)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M',
        )

    print(t.acc2taxid(accession, mapping_tsv))

@cli.command()
@click.option('-d', '--database',
              help='path of taxonomy_db/',
              required=False,
              default=None,
              type=str)
@click.option('--acc',
              help='update accession2taxid data',
              is_flag=True,
              default=False)
@click.option('--debug',
              help='debug mode',
              is_flag=True,
              default=False)

def update(database, acc, debug):
    if debug:
        logging.basicConfig(
            level=logging.DEBUG,
            format='%(asctime)s [%(levelname)s] %(module)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M',
        )
    else:
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s [%(levelname)s] %(module)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M',
        )
    t.NCBITaxonomyDownload(database, accession2taxid=acc)


if __name__ == '__main__':
    cli()