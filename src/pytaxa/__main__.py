#! python
import os, sys, logging
import argparse as ap
from . import taxonomy as t
from . import __version__
import click

@click.group(help=f"""PyTaxa taxonomy utility v{t.__version__}.""")
def pytaxacli():
    pass

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(module)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M',
)

@pytaxacli.command('query')

@click.option('-i', '--id',
              help='a taxonomy id or accession# (if available)',
              required=True,
              type=str)
@click.option('-d', '--database',
              help='path of taxonomy_db/',
              required=False,
              default=None,
              type=str)

def demo(database, id):
    t.loadTaxonomy( database )
    inid = id

    if inid[0] in "1234567890":
        taxid = inid
    else:
        taxid = t.acc2taxid( inid )
        print( "acc2taxid( %s ) => %s"   % (inid, taxid) )

    if taxid:
        print( "taxid2name( %s )                 => %s" % (taxid, t.taxid2name(taxid)) )
        print( "taxid2rank( %s )                 => %s" % (taxid, t.taxid2rank(taxid)) )
        print( "taxid2type( %s )                 => %s" % (taxid, t.taxid2type(taxid)) )
        print( "taxid2depth( %s )                => %s" % (taxid, t.taxid2depth(taxid)) )
        print( "taxid2parent( %s )               => %s" % (taxid, t.taxid2parent(taxid)) )
        print( "taxidIsLeaf( %s )                => %s" % (taxid, t.taxidIsLeaf(taxid)) )
        print( "taxid2nearestMajorTaxid( %s )    => %s" % (taxid, t.taxid2nearestMajorTaxid(taxid)) )
        print( "taxid2nameOnRank( %s, 'genus')   => %s" % (taxid, t.taxid2nameOnRank(taxid, "genus")) )
        print( "taxid2taxidOnRank( %s, 'genus')  => %s" % (taxid, t.taxid2taxidOnRank(taxid, "genus")) )
        print( "taxid2nameOnRank( %s, 'phylum')  => %s" % (taxid, t.taxid2nameOnRank(taxid, "phylum")) )
        print( "taxid2taxidOnRank( %s, 'phylum') => %s" % (taxid, t.taxid2taxidOnRank(taxid, "phylum")) )
        print( "taxid2lineage( %s )              => %s" % (taxid, t.taxid2lineage(taxid)) )
        print( "taxid2lineageDICT( %s, 1, 1 )    => %s" % (taxid, t.taxid2lineageDICT(taxid,1,1)) )
        print( "taxid2lineageTEXT( %s, 1, 1 )    => %s" % (taxid, t.taxid2lineageTEXT(taxid,1,1)) )
        print( "taxid2fullLineage( %s )          => %s" % (taxid, t.taxid2fullLineage(taxid)) )
        print( "taxid2fullLinkDict( %s )         => %s" % (taxid, t.taxid2fullLinkDict(taxid)) )
    else:
        print( "No taxid found." )

if __name__ == '__main__':
    pytaxacli()