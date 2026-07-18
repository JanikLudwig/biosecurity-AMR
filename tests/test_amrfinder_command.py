import pytest
from pathlib import Path
from genome_firewall.run_amrfinder import build_docker_command, build_native_command, normalize_genome_id

def test_build_docker_command():
    fasta_path = Path("C:/data/genomes/test.fa")
    output_temp = Path("C:/data/out/test.tsv.part")
    image = "ncbi/amr:4.2.7-2026-05-15.1"
    organism = "Staphylococcus_aureus"
    genome_id = "test"
    threads = 4
    
    cmd = build_docker_command(
        fasta_path=fasta_path,
        output_temp=output_temp,
        image=image,
        organism=organism,
        genome_id=genome_id,
        threads=threads,
        plus=False
    )
    
    assert cmd[0] == "docker"
    assert cmd[1] == "run"
    assert "--rm" in cmd
    assert "--entrypoint" in cmd
    assert cmd[cmd.index("--entrypoint") + 1] == "amrfinder"
    assert image in cmd
    
    # Mount checks
    mounts = [c for c in cmd if c.startswith("type=bind,source=")]
    assert len(mounts) == 2
    
    input_mount = [m for m in mounts if "target=/input" in m][0]
    assert "readonly" in input_mount
    assert str(fasta_path.parent.resolve()) in input_mount
    
    output_mount = [m for m in mounts if "target=/output" in m][0]
    assert "readonly" not in output_mount
    assert str(output_temp.parent.resolve()) in output_mount
    
    assert "-n" in cmd
    assert cmd[cmd.index("-n") + 1] == "/input/test.fa"
    
    assert "-o" in cmd
    assert cmd[cmd.index("-o") + 1] == "/output/test.tsv.part"
    
    assert "-O" in cmd
    assert cmd[cmd.index("-O") + 1] == organism
    
    assert "--name" in cmd
    assert cmd[cmd.index("--name") + 1] == genome_id
    
    assert "--threads" in cmd
    assert cmd[cmd.index("--threads") + 1] == str(threads)
    
    assert "--plus" not in cmd

def test_build_docker_command_plus():
    cmd = build_docker_command(
        fasta_path=Path("dummy"),
        output_temp=Path("dummy"),
        image="ncbi/amr:4.2.7",
        organism="dummy",
        genome_id="dummy",
        threads=1,
        plus=True
    )
    assert "--plus" in cmd

def test_build_native_command():
    fasta_path = Path("/data/genomes/test.fa")
    output_temp = Path("/data/out/test.tsv.part")
    organism = "Staphylococcus_aureus"
    genome_id = "test"
    threads = 4
    
    cmd = build_native_command(
        fasta_path=fasta_path,
        output_temp=output_temp,
        organism=organism,
        genome_id=genome_id,
        threads=threads,
        plus=False
    )
    
    assert cmd[0] == "amrfinder"
    assert "-n" in cmd
    assert cmd[cmd.index("-n") + 1] == str(fasta_path.resolve())
    assert "-o" in cmd
    assert cmd[cmd.index("-o") + 1] == str(output_temp.resolve())

def test_normalize_genome_id():
    assert normalize_genome_id("NCTC8325") == "NCTC8325"
    assert normalize_genome_id("Strain A-1_b.2") == "Strain_A-1_b.2"
    assert normalize_genome_id("Unsafe$Char!") == "Unsafe_Char_"
