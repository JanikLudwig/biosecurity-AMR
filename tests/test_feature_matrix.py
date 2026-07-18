import pytest
from genome_firewall.build_feature_matrix import determine_feature, extract_mutation

def test_extract_mutation():
    assert extract_mutation("S84L") == "S84L"
    assert extract_mutation("mutation in S84L") == "S84L"
    assert extract_mutation("H481Y") == "H481Y"
    assert extract_mutation("Some text without mut") == ""

def test_determine_feature_gene():
    row = {
        'Gene symbol': 'mecA',
        'Element name': '',
        'Closest reference accession': '',
        'Element type': 'AMR'
    }
    feat_id, f_type, gene_sym, mut = determine_feature(row)
    assert feat_id == 'gene::mecA'
    assert f_type == 'gene'
    assert gene_sym == 'mecA'
    
def test_determine_feature_mutation():
    row = {
        'Gene symbol': 'gyrA',
        'Element name': 'S84L',
        'Closest reference accession': '',
        'Element type': 'POINT'
    }
    feat_id, f_type, gene_sym, mut = determine_feature(row)
    assert feat_id == 'mutation::gyrA::S84L'
    assert f_type == 'mutation'
    assert gene_sym == 'gyrA'
    assert mut == 'S84L'

def test_determine_feature_fallback():
    row = {
        'Gene symbol': '',
        'Element name': 'blaZ_family',
        'Closest reference accession': '',
        'Element type': 'AMR'
    }
    feat_id, f_type, gene_sym, mut = determine_feature(row)
    assert feat_id == 'gene::blaZ_family'
    
def test_determine_feature_synonyms():
    # Test using synonyms
    row = {
        'Gene': 'blaZ',
        'Sequence name': 'blaZ',
        'Closest reference name': 'blaZ',
        'Element type': 'AMR'
    }
    feat_id, f_type, gene_sym, mut = determine_feature(row)
    assert feat_id == 'gene::blaZ'
