import argparse
import logging
import os
import sys
import csv
import gzip
import json
import re
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Any, Set, Tuple, Optional

from genome_firewall.config import load_config, get_repo_root

logger = logging.getLogger(__name__)

# Column synonyms fallback
COLUMN_SYNONYMS = {
    'gene_symbol': ['Gene symbol', 'Gene'],
    'element_name': ['Element name', 'Sequence name'],
    'reference_accession': ['Closest reference accession', 'Closest reference name'],
    'element_type': ['Element type'],
    'element_subtype': ['Element subtype'],
    'class': ['Class'],
    'subclass': ['Subclass'],
    'method': ['Method']
}

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[logging.StreamHandler(sys.stdout)]
    )

def get_column(row: Dict[str, str], synonyms: List[str]) -> str:
    for syn in synonyms:
        if syn in row:
            return row[syn].strip()
    return ""

def extract_mutation(text: str) -> str:
    """Attempts to find a canonical amino acid substitution like S84L."""
    if not text:
        return ""
    # Look for patterns like S84L, H481Y, etc. (amino acid - pos - amino acid)
    match = re.search(r'\b([A-Z]\d+[A-Z])\b', text)
    if match:
        return match.group(1)
    return ""

def determine_feature(row: Dict[str, str]) -> Tuple[str, str, str, str]:
    """
    Returns: feature_id, feature_type, gene_symbol, mutation
    """
    gene = get_column(row, COLUMN_SYNONYMS['gene_symbol'])
    elem = get_column(row, COLUMN_SYNONYMS['element_name'])
    ref = get_column(row, COLUMN_SYNONYMS['reference_accession'])
    elem_type = get_column(row, COLUMN_SYNONYMS['element_type']).upper()
    
    mutation_cand = extract_mutation(elem) or extract_mutation(gene)
    
    if 'MUTATION' in elem_type or mutation_cand:
        f_type = 'mutation'
        base_gene = gene if gene else (ref if ref else elem)
        mut = mutation_cand if mutation_cand else "unknown_mut"
        feat_id = f"mutation::{base_gene}::{mut}"
        return feat_id, f_type, base_gene, mut
    else:
        f_type = 'gene'
        base_gene = gene if gene else (ref if ref else elem)
        if not base_gene:
            base_gene = "unknown_gene"
        feat_id = f"gene::{base_gene}"
        return feat_id, f_type, base_gene, ""

def load_genomes_to_process(
    manifest_path: Optional[Path], 
    run_report: Optional[Path], 
    results_dir: Path
) -> List[str]:
    genomes = set()
    
    if manifest_path and manifest_path.exists():
        logger.info(f"Using manifest: {manifest_path}")
        with open(manifest_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                if 'genome_id' in row:
                    genomes.add(row['genome_id'].strip())
        return sorted(list(genomes))
        
    if run_report and run_report.exists():
        logger.info(f"Using run report: {run_report}")
        with open(run_report, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                st = row.get('status', '')
                if st in ('success', 'skipped'):
                    genomes.add(row['genome_id'].strip())
        if genomes:
            return sorted(list(genomes))
            
    logger.info(f"Falling back to discovering TSV files in {results_dir}")
    if results_dir.exists():
        for entry in results_dir.glob("*.tsv"):
            genomes.add(entry.stem)
            
    return sorted(list(genomes))

def main():
    parser = argparse.ArgumentParser(description="Feature Matrix Builder")
    parser.add_argument('--config', type=str, default="config/pipeline.yaml")
    parser.add_argument('--results-dir', type=str)
    parser.add_argument('--run-report', type=str)
    parser.add_argument('--manifest', type=str)
    parser.add_argument('--output-dir', type=str)
    parser.add_argument('--force', action='store_true')
    
    args = parser.parse_args()
    setup_logging()
    
    try:
        config = load_config(args.config)
    except Exception as e:
        logger.error(f"Config error: {e}")
        sys.exit(1)
        
    repo_root = get_repo_root()
    results_dir = repo_root / (args.results_dir or config['paths']['amrfinder_results'])
    output_dir = repo_root / (args.output_dir or config['paths']['processed'])
    reports_dir = repo_root / config['paths']['reports']
    
    manifest_path = repo_root / args.manifest if args.manifest else None
    run_report = reports_dir / 'amrfinder_runs.csv' if not args.run_report else repo_root / args.run_report
    
    genomes = load_genomes_to_process(manifest_path, run_report, results_dir)
    if not genomes:
        logger.error("No genomes found to process.")
        sys.exit(1)
        
    logger.info(f"Identified {len(genomes)} genomes to process.")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    feat_matrix_file = output_dir / 'features.csv.gz'
    long_file = output_dir / 'amrfinder_hits_long.tsv.gz'
    dict_file = output_dir / 'feature_dictionary.csv'
    summary_file = reports_dir / 'feature_matrix_summary.json'
    
    if feat_matrix_file.exists() and not args.force:
        logger.error("Output files already exist. Use --force to overwrite.")
        sys.exit(1)
        
    feature_dict = {}
    long_records = []
    matrix_data = {g: set() for g in genomes}
    
    zero_hits = 0
    skipped = 0
    parser_warnings = []
    
    for g_id in genomes:
        tsv_path = results_dir / f"{g_id}.tsv"
        if not tsv_path.exists():
            skipped += 1
            parser_warnings.append(f"Missing TSV for {g_id}")
            continue
            
        try:
            with open(tsv_path, 'r', encoding='utf-8') as f:
                content = f.read().strip()
                if not content:
                    zero_hits += 1
                    continue
                    
                lines = content.split('\n')
                if len(lines) == 1:
                    # Header only
                    zero_hits += 1
                    continue
                    
                reader = csv.DictReader(lines, delimiter='\t')
                if not reader.fieldnames:
                    parser_warnings.append(f"No valid header in {g_id}.tsv")
                    continue
                    
                hits_found = False
                for row in reader:
                    feat_id, f_type, gene_sym, mut = determine_feature(row)
                    matrix_data[g_id].add(feat_id)
                    hits_found = True
                    
                    if feat_id not in feature_dict:
                        feature_dict[feat_id] = {
                            'feature_id': feat_id,
                            'feature_type': f_type,
                            'gene_symbol': gene_sym,
                            'mutation': mut,
                            'element_name': get_column(row, COLUMN_SYNONYMS['element_name']),
                            'reference_accession': get_column(row, COLUMN_SYNONYMS['reference_accession']),
                            'element_type': get_column(row, COLUMN_SYNONYMS['element_type']),
                            'element_subtype': get_column(row, COLUMN_SYNONYMS['element_subtype']),
                            'drug_class': get_column(row, COLUMN_SYNONYMS['class']),
                            'drug_subclass': get_column(row, COLUMN_SYNONYMS['subclass']),
                            'source_columns': json.dumps(list(row.keys()))
                        }
                    
                    # Store long record
                    long_record = {
                        'genome_id': g_id,
                        'feature_id': feat_id,
                        'feature_type': f_type,
                        'gene_symbol': gene_sym,
                        'mutation': mut,
                        'element_name': feature_dict[feat_id]['element_name'],
                        'reference_accession': feature_dict[feat_id]['reference_accession'],
                        'element_type': feature_dict[feat_id]['element_type'],
                        'element_subtype': feature_dict[feat_id]['element_subtype'],
                        'drug_class': feature_dict[feat_id]['drug_class'],
                        'drug_subclass': feature_dict[feat_id]['drug_subclass'],
                        'source_file': tsv_path.name
                    }
                    # Add all original columns to long record safely
                    for k, v in row.items():
                        if k not in long_record:
                            long_record[k] = v
                    long_records.append(long_record)
                    
                if not hits_found:
                    zero_hits += 1
                    
        except Exception as e:
            skipped += 1
            parser_warnings.append(f"Error parsing {g_id}.tsv: {e}")
            
    # Write Feature Dictionary
    all_features = sorted(list(feature_dict.keys()))
    if all_features:
        with open(dict_file, 'w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(feature_dict[all_features[0]].keys()))
            writer.writeheader()
            for fid in all_features:
                writer.writerow(feature_dict[fid])
                
    # Write Matrix
    with gzip.open(feat_matrix_file, 'wt', encoding='utf-8', newline='') as f:
        fieldnames = ['genome_id'] + all_features
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for g_id in sorted(genomes):
            row = {'genome_id': g_id}
            for fid in all_features:
                row[fid] = 1 if fid in matrix_data[g_id] else 0
            writer.writerow(row)
            
    # Write Long Form
    if long_records:
        # Determine all possible keys for long format
        long_keys = set()
        for lr in long_records:
            long_keys.update(lr.keys())
        # Sort keys prioritizing our standardized columns
        std_cols = [
            'genome_id', 'feature_id', 'feature_type', 'gene_symbol', 'mutation',
            'element_name', 'reference_accession', 'element_type', 'element_subtype',
            'drug_class', 'drug_subclass', 'source_file'
        ]
        other_cols = sorted(list(long_keys - set(std_cols)))
        final_long_cols = std_cols + other_cols
        
        with gzip.open(long_file, 'wt', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=final_long_cols, delimiter='\t', extrasaction='ignore')
            writer.writeheader()
            writer.writerows(long_records)
            
    # Write Summary
    summary = {
        'created_at': datetime.now().isoformat(),
        'number_of_genomes': len(genomes),
        'number_of_features': len(all_features),
        'number_of_gene_features': sum(1 for f in all_features if feature_dict[f]['feature_type'] == 'gene'),
        'number_of_mutation_features': sum(1 for f in all_features if feature_dict[f]['feature_type'] == 'mutation'),
        'number_of_zero_hit_genomes': zero_hits,
        'number_of_input_files': len(genomes),
        'skipped_or_failed_inputs': skipped,
        'parser_warnings': parser_warnings,
        'output_files': [
            str(feat_matrix_file.relative_to(repo_root)),
            str(long_file.relative_to(repo_root)),
            str(dict_file.relative_to(repo_root)),
            str(summary_file.relative_to(repo_root))
        ]
    }
    
    with open(summary_file, 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=4)
        
    logger.info("Feature matrix generated successfully.")

if __name__ == '__main__':
    main()
