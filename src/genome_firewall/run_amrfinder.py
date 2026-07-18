import argparse
import logging
import os
import shutil
import subprocess
import sys
import csv
import json
import re
import platform
from pathlib import Path
from datetime import datetime
from typing import List, Optional, Tuple, Dict, Any

from genome_firewall.config import load_config, get_repo_root

logger = logging.getLogger(__name__)

def setup_logging(log_dir: Path):
    log_dir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )

def setup_genome_logger(genome_id: Path, log_dir: Path) -> logging.Logger:
    g_logger = logging.getLogger(f"genome.{genome_id}")
    g_logger.setLevel(logging.INFO)
    # Clear existing handlers
    if g_logger.hasHandlers():
        g_logger.handlers.clear()
    
    log_file = log_dir / f"{genome_id}.log"
    fh = logging.FileHandler(log_file, encoding='utf-8')
    fh.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s'))
    g_logger.addHandler(fh)
    return g_logger

def normalize_genome_id(name: str) -> str:
    """Replaces unsafe characters with underscores."""
    return re.sub(r'[^a-zA-Z0-9._-]', '_', name)

def select_backend(backend: str) -> str:
    if backend == 'auto':
        if platform.system() == 'Windows':
            return 'docker'
        else:
            if shutil.which('amrfinder'):
                return 'native'
            else:
                return 'docker'
    
    if backend == 'native':
        if not shutil.which('amrfinder'):
            raise RuntimeError("Backend 'native' requested, but 'amrfinder' is not in PATH.")
        return 'native'
    
    if backend == 'docker':
        if not shutil.which('docker'):
            raise RuntimeError("Backend 'docker' requested, but 'docker' is not in PATH.")
        
        # Test if docker is running
        try:
            subprocess.run(['docker', 'info'], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            raise RuntimeError("Backend 'docker' requested, but Docker daemon is not running.")
        return 'docker'
    
    raise ValueError(f"Unknown backend: {backend}")

def load_manifest(manifest_path: Path) -> List[Tuple[str, Path]]:
    if not manifest_path.exists():
        raise FileNotFoundError(f"Manifest not found: {manifest_path}")
    
    genomes = []
    seen_ids = set()
    manifest_dir = manifest_path.parent
    
    with open(manifest_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if 'genome_id' not in row or 'fasta_path' not in row:
                raise ValueError("Manifest must contain 'genome_id' and 'fasta_path' columns.")
            
            raw_id = row['genome_id'].strip()
            if not raw_id:
                raise ValueError("Empty genome_id found in manifest.")
            
            fasta_path = Path(row['fasta_path'])
            if not fasta_path.is_absolute():
                fasta_path = manifest_dir / fasta_path
            
            norm_id = normalize_genome_id(raw_id)
            if norm_id in seen_ids:
                raise ValueError(f"Duplicate normalized genome_id found in manifest: {norm_id}")
            seen_ids.add(norm_id)
            
            if not fasta_path.exists():
                raise FileNotFoundError(f"FASTA file not found for {raw_id}: {fasta_path}")
            
            genomes.append((raw_id, fasta_path))
    
    return sorted(genomes, key=lambda x: x[0])

def discover_fastas(input_dir: Path) -> List[Tuple[str, Path]]:
    if not input_dir.exists():
        return []
    
    extensions = {'.fa', '.fna', '.fasta'}
    genomes = []
    seen_ids = set()
    
    # Check for unsupported compressed files
    for entry in input_dir.iterdir():
        if entry.is_file():
            if entry.name.endswith('.gz') and any(entry.name.endswith(ext + '.gz') for ext in extensions):
                logger.error(f"Compressed FASTA not supported: {entry.name}. Please decompress it first.")
                sys.exit(1)
            
            if entry.suffix in extensions:
                raw_id = entry.stem
                norm_id = normalize_genome_id(raw_id)
                if norm_id in seen_ids:
                    logger.error(f"Identifier collision after normalization for {raw_id} -> {norm_id}")
                    sys.exit(1)
                seen_ids.add(norm_id)
                genomes.append((raw_id, entry))
    
    return sorted(genomes, key=lambda x: x[0])

def build_docker_command(
    fasta_path: Path, 
    output_temp: Path, 
    image: str, 
    organism: str, 
    genome_id: str, 
    threads: int, 
    plus: bool
) -> List[str]:
    input_dir = fasta_path.parent.resolve()
    output_dir = output_temp.parent.resolve()
    
    fasta_name = fasta_path.name
    out_name = output_temp.name
    
    cmd = [
        "docker", "run", "--rm",
        "--entrypoint", "amrfinder",
        "--mount", f"type=bind,source={input_dir},target=/input,readonly",
        "--mount", f"type=bind,source={output_dir},target=/output",
        image,
        "-n", f"/input/{fasta_name}",
        "-O", organism,
        "--name", genome_id,
        "--threads", str(threads),
        "-o", f"/output/{out_name}"
    ]
    if plus:
        cmd.append("--plus")
    return cmd

def build_native_command(
    fasta_path: Path, 
    output_temp: Path, 
    organism: str, 
    genome_id: str, 
    threads: int, 
    plus: bool
) -> List[str]:
    cmd = [
        "amrfinder",
        "-n", str(fasta_path.resolve()),
        "-O", organism,
        "--name", genome_id,
        "--threads", str(threads),
        "-o", str(output_temp.resolve())
    ]
    if plus:
        cmd.append("--plus")
    return cmd

def validate_existing_output(output_path: Path) -> bool:
    if not output_path.exists():
        return False
    try:
        with open(output_path, 'r', encoding='utf-8') as f:
            first_line = f.readline()
            if not first_line:
                return False
            return True # Header only or with data is valid
    except Exception:
        return False

def run_single_genome(
    raw_id: str,
    fasta_path: Path,
    output_dir: Path,
    log_dir: Path,
    backend: str,
    image: str,
    organism: str,
    threads: int,
    plus: bool,
    force: bool
) -> Dict[str, Any]:
    norm_id = normalize_genome_id(raw_id)
    final_output = output_dir / f"{norm_id}.tsv"
    temp_output = output_dir / f"{norm_id}.tsv.part"
    
    result = {
        'original_genome_id': raw_id,
        'genome_id': norm_id,
        'input_fasta': str(fasta_path),
        'backend': backend,
        'docker_image': image if backend == 'docker' else '',
        'organism': organism,
        'plus_enabled': str(plus).lower(),
        'command': '',
        'started_at': '',
        'finished_at': '',
        'runtime_seconds': '',
        'exit_code': '',
        'status': '',
        'output_tsv': str(final_output),
        'log_file': str(log_dir / f"{norm_id}.log"),
        'error_message': '',
        'amrfinder_version': '',
        'database_version': ''
    }
    
    g_logger = setup_genome_logger(norm_id, log_dir)
    g_logger.info(f"Processing genome: {raw_id} (normalized: {norm_id})")
    
    if not force and validate_existing_output(final_output):
        g_logger.info("Valid output already exists. Skipping.")
        result['status'] = 'skipped'
        return result
    
    if backend == 'docker':
        cmd = build_docker_command(fasta_path, temp_output, image, organism, norm_id, threads, plus)
    else:
        cmd = build_native_command(fasta_path, temp_output, organism, norm_id, threads, plus)
    
    result['command'] = json.dumps(cmd)
    g_logger.info(f"Backend: {backend}")
    g_logger.info(f"Command: {' '.join(cmd)}")
    
    start_time = datetime.now()
    result['started_at'] = start_time.isoformat()
    
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        end_time = datetime.now()
        result['finished_at'] = end_time.isoformat()
        result['runtime_seconds'] = str(round((end_time - start_time).total_seconds(), 2))
        result['exit_code'] = str(proc.returncode)
        
        g_logger.info(f"Exit code: {proc.returncode}")
        g_logger.info(f"Runtime: {result['runtime_seconds']}s")
        if proc.stdout:
            g_logger.info(f"STDOUT:\n{proc.stdout.strip()}")
        if proc.stderr:
            g_logger.info(f"STDERR:\n{proc.stderr.strip()}")
            
        if proc.returncode == 0 and validate_existing_output(temp_output):
            # Parse version if possible from stdout or stderr
            ver_match = re.search(r'amrfinder\s+([0-9.]+)', proc.stdout + proc.stderr)
            if ver_match:
                result['amrfinder_version'] = ver_match.group(1)
            
            db_match = re.search(r'Database:\s+([0-9-.]+)', proc.stdout + proc.stderr)
            if db_match:
                result['database_version'] = db_match.group(1)
            
            temp_output.replace(final_output)
            result['status'] = 'success'
            g_logger.info("Successfully completed.")
        else:
            result['status'] = 'failed'
            err_msg = "Process failed or invalid output."
            result['error_message'] = err_msg
            g_logger.error(err_msg)
            if temp_output.exists():
                temp_output.unlink()
    
    except Exception as e:
        end_time = datetime.now()
        result['finished_at'] = end_time.isoformat()
        result['runtime_seconds'] = str(round((end_time - start_time).total_seconds(), 2))
        result['status'] = 'failed'
        result['error_message'] = str(e)
        g_logger.error(f"Execution exception: {e}")
        if temp_output.exists():
            temp_output.unlink()
    
    return result

def write_run_report(reports_dir: Path, results: List[Dict[str, Any]]):
    reports_dir.mkdir(parents=True, exist_ok=True)
    report_file = reports_dir / 'amrfinder_runs.csv'
    
    if not results:
        return
        
    fieldnames = list(results[0].keys())
    write_header = not report_file.exists()
    
    with open(report_file, 'a', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        writer.writerows(results)

def main():
    parser = argparse.ArgumentParser(description="Plattformübergreifender AMRFinderPlus-Runner")
    parser.add_argument('--backend', choices=['auto', 'docker', 'native'], help="Backend selection")
    parser.add_argument('--config', type=str, default="config/pipeline.yaml", help="Path to config file")
    parser.add_argument('--image', type=str, help="Docker image")
    parser.add_argument('--organism', type=str, help="Organism for AMRFinderPlus")
    parser.add_argument('--threads', type=int, help="Number of threads")
    parser.add_argument('--input-dir', type=str, help="Input directory for FASTA files")
    parser.add_argument('--output-dir', type=str, help="Output directory")
    parser.add_argument('--log-dir', type=str, help="Log directory")
    parser.add_argument('--report', type=str, help="Path for report CSV (overrides default reports dir)")
    parser.add_argument('--manifest', type=str, help="Path to manifest CSV")
    parser.add_argument('--force', action='store_true', help="Overwrite existing valid outputs")
    parser.add_argument('--limit', type=int, help="Limit number of processed genomes")
    parser.add_argument('--fail-fast', action='store_true', help="Stop on first failure")
    parser.add_argument('--plus', action='store_true', help="Enable --plus flag for AMRFinderPlus")
    
    args = parser.parse_args()
    
    try:
        config = load_config(args.config)
    except Exception as e:
        logger.error(f"Config error: {e}")
        sys.exit(1)
        
    repo_root = get_repo_root()
    
    # Resolve parameters (CLI overrides config)
    backend_pref = args.backend or config['amrfinder']['backend']
    image = args.image or config['amrfinder']['docker_image']
    organism = args.organism or config['amrfinder']['organism']
    threads = args.threads if args.threads is not None else config['amrfinder']['threads']
    plus = args.plus or config['amrfinder'].get('plus', False)
    
    input_dir = repo_root / (args.input_dir or config['paths']['genomes'])
    output_dir = repo_root / (args.output_dir or config['paths']['amrfinder_results'])
    log_dir = repo_root / (args.log_dir or config['paths']['amrfinder_logs'])
    reports_dir = repo_root / (args.report if args.report else config['paths']['reports'])
    
    if threads <= 0:
        logger.error("Threads must be a positive integer.")
        sys.exit(1)
        
    setup_logging(log_dir)
    
    try:
        backend = select_backend(backend_pref)
    except RuntimeError as e:
        logger.error(e)
        sys.exit(1)
        
    logger.info(f"Resolved backend: {backend}")
    
    if args.manifest:
        manifest_path = repo_root / args.manifest
        try:
            genomes = load_manifest(manifest_path)
        except Exception as e:
            logger.error(f"Manifest error: {e}")
            sys.exit(1)
    else:
        genomes = discover_fastas(input_dir)
        
    if args.limit and args.limit > 0:
        genomes = genomes[:args.limit]
        
    if not genomes:
        logger.info("No genomes found to process.")
        sys.exit(0)
        
    output_dir.mkdir(parents=True, exist_ok=True)
    
    results = []
    has_failure = False
    
    for raw_id, fasta_path in genomes:
        res = run_single_genome(
            raw_id=raw_id,
            fasta_path=fasta_path,
            output_dir=output_dir,
            log_dir=log_dir,
            backend=backend,
            image=image,
            organism=organism,
            threads=threads,
            plus=plus,
            force=args.force
        )
        results.append(res)
        
        if res['status'] == 'failed':
            has_failure = True
            if args.fail_fast:
                logger.error(f"Fail-fast triggered by genome {raw_id}")
                break
                
    write_run_report(reports_dir, results)
    
    if has_failure:
        logger.error("One or more processing steps failed.")
        sys.exit(1)
    else:
        logger.info("All processed successfully.")
        sys.exit(0)

if __name__ == '__main__':
    main()
