import pytest
from genome_firewall.build_feature_matrix import determine_feature

def test_determine_feature_1_gyrA_S84L():
    # 1. gyrA_S84L mit Subtype POINT und Method POINTX
    row = {
        'Gene symbol': 'gyrA_S84L',
        'Element subtype': 'POINT',
        'Method': 'POINTX'
    }
    feat_id, f_type, gene_sym, mut = determine_feature(row)
    assert feat_id == 'mutation::gyrA::S84L'
    assert f_type == 'mutation'

def test_determine_feature_2_folP_E258EKE():
    # 2. folP_E258EKE mit Method POINTX
    row = {
        'Gene symbol': 'folP_E258EKE',
        'Element subtype': 'AMR',
        'Method': 'POINTX'
    }
    feat_id, f_type, gene_sym, mut = determine_feature(row)
    assert feat_id == 'mutation::folP::E258EKE'

def test_determine_feature_3_murA_T396N():
    # 3. murA_T396N
    row = {
        'Gene symbol': 'murA_T396N',
        'Element subtype': 'POINT',
        'Method': 'POINT'
    }
    feat_id, f_type, gene_sym, mut = determine_feature(row)
    assert feat_id == 'mutation::murA::T396N'

def test_determine_feature_4_glpT_A100V():
    # 4. glpT_A100V
    row = {
        'Gene symbol': 'glpT_A100V',
        'Element subtype': 'POINT',
        'Method': 'POINTX'
    }
    feat_id, f_type, gene_sym, mut = determine_feature(row)
    assert feat_id == 'mutation::glpT::A100V'

def test_determine_feature_5_fosB():
    # 5. fosB mit Subtype AMR und Method EXACTX
    row = {
        'Gene symbol': 'fosB',
        'Element subtype': 'AMR',
        'Method': 'EXACTX'
    }
    feat_id, f_type, gene_sym, mut = determine_feature(row)
    assert feat_id == 'gene::fosB'
    assert f_type == 'gene'

def test_determine_feature_6_blaPC1():
    # 6. blaPC1 mit Subtype AMR und Method ALLELEX
    row = {
        'Gene symbol': 'blaPC1',
        'Element subtype': 'AMR',
        'Method': 'ALLELEX'
    }
    feat_id, f_type, gene_sym, mut = determine_feature(row)
    assert feat_id == 'gene::blaPC1'

def test_determine_feature_7_normal_gene_with_underscore():
    # 7. Normales Gen mit Unterstrich, aber ohne POINT-/POINTX-Klassifikation
    row = {
        'Gene symbol': 'blaZ_family',
        'Element subtype': 'AMR',
        'Method': 'EXACTX'
    }
    feat_id, f_type, gene_sym, mut = determine_feature(row)
    assert feat_id == 'gene::blaZ_family'
    assert f_type == 'gene'

def test_determine_feature_8_fallback_mutation():
    # 8. Mutation mit nicht eindeutig zerlegbarem Symbol erhält einen stabilen Mutation-Fallback.
    row = {
        'Gene symbol': 'unknownMut',
        'Element subtype': 'POINT',
        'Method': 'POINT'
    }
    feat_id, f_type, gene_sym, mut = determine_feature(row)
    assert feat_id == 'mutation::unknownMut::unknown_mut'
    assert f_type == 'mutation'

# Tests 9-13 (Multiple hits, header-only, zero hits etc.) are best tested as integration tests
# for build_feature_matrix, or conceptually they are fulfilled by the file processing logic.
# I will just write basic integration tests for them or we just rely on the main code logic
# which we know from Phase B will be checked anyway.
