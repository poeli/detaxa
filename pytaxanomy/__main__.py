#! python
import os, sys, logging
import argparse as ap

__version__ = '0.1.0'

if __name__ == '__main__':
    from .taxonomy import Taxonomy

    t = Taxonomy(
        dbpath=sys.argv[1] if len(sys.argv)>1 else None,
        cus_taxonomy_file=sys.argv[2] if len(sys.argv)>2 else None
    )

    #loading taxonomy
    t.loadTaxonomy()

    inid = 0
    try:
        inid = input("\nEnter acc/taxid: ")
    except:
        inid = 0

    while(inid):
        if inid[0] in "1234567890":
            taxid = inid
        else:
            taxid = acc2taxid( inid )
            print( "acc2taxid( %s ) => %s"   % (inid, taxid) )

        if taxid:
            print( "taxid2name( %s )                 => %s" % (taxid, taxid2name(taxid)) )
            print( "taxid2rank( %s )                 => %s" % (taxid, taxid2rank(taxid)) )
            print( "taxid2type( %s )                 => %s" % (taxid, taxid2type(taxid)) )
            print( "taxid2parent( %s )               => %s" % (taxid, taxid2parent(taxid)) )
            print( "taxidIsLeaf( %s )                => %s" % (taxid, taxidIsLeaf(taxid)) )
            print( "taxid2nearestMajorTaxid( %s )    => %s" % (taxid, taxid2nearestMajorTaxid(taxid)) )
            print( "taxid2nameOnRank( %s, 'genus')   => %s" % (taxid, taxid2nameOnRank(taxid, "genus")) )
            print( "taxid2taxidOnRank( %s, 'genus')  => %s" % (taxid, taxid2taxidOnRank(taxid, "genus")) )
            print( "taxid2nameOnRank( %s, 'phylum')  => %s" % (taxid, taxid2nameOnRank(taxid, "phylum")) )
            print( "taxid2taxidOnRank( %s, 'phylum') => %s" % (taxid, taxid2taxidOnRank(taxid, "phylum")) )
            print( "taxid2lineage( %s )              => %s" % (taxid, taxid2lineage(taxid)) )
            print( "taxid2lineageDICT( %s, 1, 1 )    => %s" % (taxid, taxid2lineageDICT(taxid,1,1)) )
            print( "taxid2lineageTEXT( %s, 1, 1 )    => %s" % (taxid, taxid2lineageTEXT(taxid,1,1)) )
            print( "taxid2fullLineage( %s )          => %s" % (taxid, taxid2fullLineage(taxid)) )
            print( "taxid2fullLinkDict( %s )         => %s" % (taxid, taxid2fullLinkDict(taxid)) )
        else:
            print( "No taxid found.\n" )

        try:
            inid = input("\nEnter acc/taxid: ")
        except:
            inid = 0
