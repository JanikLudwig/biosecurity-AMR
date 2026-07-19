import argparse
import json
import pandas as pd
import numpy as np
from pathlib import Path

def safe_display_path(path_str):
    if not path_str:
        return ""
    
    path = Path(path_str).resolve()
    repo_root = Path(__file__).resolve().parent.parent
    
    try:
        rel = path.relative_to(repo_root)
        return rel.as_posix()
    except ValueError:
        return path.name

def main(args=None):
    parser = argparse.ArgumentParser(description="Generate Pilot Report")
    parser.add_argument('--labels', type=str, default="data/processed/aureus_pilot_100_labels.csv.gz", help="Path to processed labels CSV")
    parser.add_argument('--manifest', type=str, default="data/manifests/aureus_pilot_100.csv", help="Path to manifest CSV")
    parser.add_argument('--run-report', type=str, default="reports/amrfinder_runs.csv", help="Path to AMRFinder run report")
    parser.add_argument('--feature-summary', type=str, default="reports/feature_matrix_summary.json", help="Path to feature matrix summary JSON")
    parser.add_argument('--features', type=str, default="data/processed/features.csv.gz", help="Path to feature matrix CSV")
    parser.add_argument('--output-json', type=str, default="reports/pilot_100/pilot_summary.json", help="Output JSON path")
    parser.add_argument('--output-markdown', type=str, default="reports/pilot_100/pilot_summary.md", help="Output Markdown path")
    
    parsed_args = parser.parse_args(args)
    
    # 1. Load data
    man_df = pd.read_csv(parsed_args.manifest, dtype=str)
    selected_gids = man_df["genome_id"].tolist()
    pilot_genome_count = len(selected_gids)
    
    labels_df = pd.read_csv(parsed_args.labels, dtype=str)
    # Filter labels to only include those in the manifest
    labels_df = labels_df[labels_df["genome_id"].isin(selected_gids)]
    unique_label_genomes = labels_df["genome_id"].nunique()
    unique_genome_antibiotic_pairs = len(labels_df)
    
    # 2. Run times
    run_df = pd.read_csv(parsed_args.run_report)
    run_df = run_df[run_df["original_genome_id"].isin(selected_gids)]
    
    success_df = run_df[run_df["status"] == "success"]
    successes = len(success_df)
    fails = len(run_df[run_df["status"] == "failed"])
    
    if successes > 0:
        runtimes = success_df["runtime_seconds"].dropna()
        min_time = runtimes.min() if len(runtimes) > 0 else 0
        max_time = runtimes.max() if len(runtimes) > 0 else 0
        mean_time = runtimes.mean() if len(runtimes) > 0 else 0
        median_time = runtimes.median() if len(runtimes) > 0 else 0
        total_time = runtimes.sum() if len(runtimes) > 0 else 0
    else:
        min_time = max_time = mean_time = median_time = total_time = 0
        
    # 3. Feature Matrix
    with open(parsed_args.feature_summary, "r") as f:
        mat_sum = json.load(f)
        
    # 4. Label distribution
    dist = {}
    for col in ["cefoxitin", "ciprofloxacin", "erythromycin"]:
        ab_df = labels_df[labels_df["antibiotic"] == col]
        counts = ab_df["label"].value_counts().to_dict()
        dist[col] = counts
        
    # 5. Top mutations
    matrix_df = pd.read_csv(parsed_args.features)
    matrix_df = matrix_df[matrix_df["genome_id"].isin(selected_gids)]
    mutation_cols = [c for c in matrix_df.columns if c.startswith("mutation::")]
    if mutation_cols and len(matrix_df) > 0:
        mut_freq = matrix_df[mutation_cols].sum().sort_values(ascending=False)
        top_mutations = mut_freq[mut_freq > 0].head(5).to_dict()
        genomes_with_mutations = int((matrix_df[mutation_cols].sum(axis=1) > 0).sum())
    else:
        top_mutations = {}
        genomes_with_mutations = 0
    
    # 6. Prepare report dict
    report_dict = {
        "data_sources": {
            "labels": safe_display_path(parsed_args.labels),
            "manifest": safe_display_path(parsed_args.manifest),
            "run_report": safe_display_path(parsed_args.run_report)
        },
        "number_of_unique_label_genomes": unique_label_genomes,
        "number_of_unique_genome_antibiotic_pairs": unique_genome_antibiotic_pairs,
        "pilot_genome_count": pilot_genome_count,
        "amrfinder_success": successes,
        "amrfinder_failed": fails,
        "runtimes_s": {
            "min": round(min_time, 2),
            "max": round(max_time, 2),
            "mean": round(mean_time, 2),
            "median": round(median_time, 2),
            "total": round(total_time, 2)
        },
        "gene_features": mat_sum.get("number_of_gene_features", 0),
        "mutation_features": mat_sum.get("number_of_mutation_features", 0),
        "feature_columns_without_id": mat_sum.get("number_of_features", 0),
        "dataframe_shape": f"{mat_sum.get('number_of_genomes', 0)}x{mat_sum.get('number_of_features', 0) + 1}",
        "zero_hit_genomes": mat_sum.get("number_of_zero_hit_genomes", 0),
        "genomes_with_mutations": genomes_with_mutations,
        "top_mutations": top_mutations,
        "label_distribution": dist
    }
    
    Path(parsed_args.output_json).parent.mkdir(parents=True, exist_ok=True)
    with open(parsed_args.output_json, "w", encoding="utf-8") as f:
        json.dump(report_dict, f, indent=4)
        
    md_text = f"""# AMRFinderPlus Pilot Report ({pilot_genome_count} Genome)

## 1. Datenquellen
- **Labels:** `{report_dict['data_sources']['labels']}`
- **Manifest:** `{report_dict['data_sources']['manifest']}`
- **Run Report:** `{report_dict['data_sources']['run_report']}`

## 2. Datenbestand (Pilot)
- **Pilotgenome:** {report_dict['pilot_genome_count']}
- **Eindeutige Genome mit Labels:** {report_dict['number_of_unique_label_genomes']}
- **Label-Paare gesamt:** {report_dict['number_of_unique_genome_antibiotic_pairs']}

## 3. AMRFinderPlus Pilotlauf
- **Erfolgsquote:** {report_dict['amrfinder_success']} erfolgreich / {report_dict['amrfinder_failed']} fehlgeschlagen
- **Laufzeiten (Sekunden pro Genom):** 
  - Min: {report_dict['runtimes_s']['min']}s | Max: {report_dict['runtimes_s']['max']}s
  - Mean: {report_dict['runtimes_s']['mean']}s | Median: {report_dict['runtimes_s']['median']}s
  - Total: {report_dict['runtimes_s']['total']}s

## 4. Feature Matrix
- **Zahl Genfeatures:** {report_dict['gene_features']}
- **Zahl Mutationsfeatures:** {report_dict['mutation_features']}
- **Zahl Feature-Spalten ohne genome_id:** {report_dict['feature_columns_without_id']}
- **DataFrame-Form einschließlich genome_id:** {report_dict['dataframe_shape']}
- **Nullzeilen (Genome ohne Treffer):** {report_dict['zero_hit_genomes']}
- **Zahl Genome mit mindestens einer Mutation:** {report_dict['genomes_with_mutations']}
- **Häufigste Mutationsfeatures:**
"""
    for m, c in top_mutations.items():
        md_text += f"  - {m}: {c}\n"

    md_text += f"\n## 5. Kohorten-Verteilung (R/S)\n"
    for ab, counts in dist.items():
        r_c = counts.get('R', 0)
        s_c = counts.get('S', 0)
        md_text += f"- **{ab.capitalize()}:** R={r_c}, S={s_c}\n"
        
    Path(parsed_args.output_markdown).parent.mkdir(parents=True, exist_ok=True)
    with open(parsed_args.output_markdown, "w", encoding="utf-8") as f:
        f.write(md_text)

if __name__ == "__main__":
    main()
