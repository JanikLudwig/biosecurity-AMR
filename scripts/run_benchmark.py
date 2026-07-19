import os
import pandas as pd
import subprocess
import time
import json
import hashlib
from pathlib import Path

def normalize_tsv(filepath):
    df = pd.read_csv(filepath, sep='\t')
    df = df.sort_values(by=list(df.columns))
    return df.to_csv(index=False, sep='\t')

def hash_tsv(normalized_content):
    return hashlib.sha256(normalized_content.encode('utf-8')).hexdigest()

def run_benchmark():
    # 1. Create subset of 20 genomes deterministically
    pilot_df = pd.read_csv("data/manifests/aureus_pilot_100.csv")
    # Take first 20 as deterministic subset
    bench_df = pilot_df.head(20)
    bench_df.to_csv("data/manifests/aureus_benchmark_20.csv", index=False)
    
    configs = [
        {"name": "A", "workers": 1, "threads": 4},
        {"name": "B", "workers": 2, "threads": 3},
        {"name": "C", "workers": 3, "threads": 2}
    ]
    
    results_summary = []
    
    for cfg in configs:
        name = cfg["name"]
        w = cfg["workers"]
        t = cfg["threads"]
        
        out_dir = f"data/interim/benchmark/{name.lower()}"
        log_dir = f"logs/benchmark/{name.lower()}"
        rep_file = f"reports/benchmark/{name.lower()}_runs.csv"
        
        cmd = [
            "python", "-m", "genome_firewall.run_amrfinder",
            "--manifest", "data/manifests/aureus_benchmark_20.csv",
            "--output-dir", out_dir,
            "--log-dir", log_dir,
            "--report", rep_file,
            "--workers", str(w),
            "--threads", str(t),
            "--force"
        ]
        
        start_time = time.time()
        print(f"Running config {name}: workers={w}, threads={t}")
        proc = subprocess.run(cmd, capture_output=True, text=True)
        end_time = time.time()
        
        total_time = end_time - start_time
        
        # Read the report
        rep_df = pd.read_csv(rep_file)
        successes = len(rep_df[rep_df['status'] == 'success'])
        fails = len(rep_df[rep_df['status'] == 'failed'])
        mean_time = rep_df['runtime_seconds'].mean()
        median_time = rep_df['runtime_seconds'].median()
        genomes_per_min = (successes / total_time) * 60 if total_time > 0 else 0
        
        results_summary.append({
            "config": name,
            "workers": w,
            "threads": t,
            "workers_x_threads": w * t,
            "total_runtime_s": round(total_time, 2),
            "successes": successes,
            "fails": fails,
            "mean_runtime_s": round(mean_time, 2),
            "median_runtime_s": round(median_time, 2),
            "genomes_per_minute": round(genomes_per_min, 2)
        })
        
    Path("reports/benchmark").mkdir(parents=True, exist_ok=True)
    pd.DataFrame(results_summary).to_csv("reports/benchmark/benchmark_results.csv", index=False)
    with open("reports/benchmark/benchmark_summary.json", "w") as f:
        json.dump(results_summary, f, indent=4)
        
    # Scientific identity check
    print("Checking scientific identity...")
    gids = bench_df["genome_id"].tolist()
    identity_results = []
    identical_count = 0
    diff_count = 0
    
    for gid in gids:
        # replace unsafe chars
        norm_id = str(gid).replace(":", "_").replace("/", "_")
        
        hashes = {}
        contents = {}
        for cfg in configs:
            name = cfg["name"]
            tsv_path = f"data/interim/benchmark/{name.lower()}/{norm_id}.tsv"
            if os.path.exists(tsv_path):
                norm_content = normalize_tsv(tsv_path)
                hashes[name] = hash_tsv(norm_content)
                contents[name] = norm_content
            else:
                hashes[name] = "MISSING"
        
        unique_hashes = set(hashes.values())
        if len(unique_hashes) == 1 and "MISSING" not in unique_hashes:
            identical_count += 1
        else:
            diff_count += 1
            identity_results.append({
                "genome_id": gid,
                "hashes": hashes
            })
            
    md_text = "# Benchmark Summary\n\n"
    for r in results_summary:
        md_text += f"## Config {r['config']}\n"
        md_text += f"- Workers: {r['workers']}\n- Threads: {r['threads']}\n"
        md_text += f"- Total Runtime: {r['total_runtime_s']}s\n"
        md_text += f"- Successes: {r['successes']}\n"
        md_text += f"- Mean Runtime: {r['mean_runtime_s']}s\n"
        md_text += f"- Genomes/Min: {r['genomes_per_minute']}\n\n"
        
    md_text += f"## Scientific Identity\n- Identical: {identical_count}\n- Different: {diff_count}\n"
    if diff_count > 0:
        md_text += "- Details:\n"
        for diff in identity_results:
            md_text += f"  - {diff['genome_id']}: {diff['hashes']}\n"
            
    with open("reports/benchmark/benchmark_summary.md", "w") as f:
        f.write(md_text)
        
    print(f"Identical: {identical_count}, Diff: {diff_count}")
    
if __name__ == "__main__":
    run_benchmark()
