"""
Basic tests for detaxa package
"""
import pytest
import detaxa
from detaxa.taxonomy import taxid2name, taxid2rank


def test_version():
    """Test that version is accessible"""
    assert hasattr(detaxa, '__version__')
    assert isinstance(detaxa.__version__, str)


def test_taxid2name_basic():
    """Test basic taxid2name functionality"""
    # Test with root taxid
    name = taxid2name("1")
    assert name == "root"
    
    # Test with unknown taxid
    name = taxid2name("999999999")
    assert name == "unknown"


def test_taxid2rank_basic():
    """Test basic taxid2rank functionality"""
    # Test with root taxid
    rank = taxid2rank("1")
    assert rank == "root"
    
    # Test with unknown taxid
    rank = taxid2rank("999999999")
    assert rank == "unknown"


def test_imports():
    """Test that main functions can be imported"""
    from detaxa.taxonomy import (
        taxid2name, taxid2rank, taxid2parent, 
        taxid2lineage, name2taxid
    )
    
    # Basic import test
    assert callable(taxid2name)
    assert callable(taxid2rank)
    assert callable(taxid2parent)
    assert callable(taxid2lineage)
    assert callable(name2taxid)
