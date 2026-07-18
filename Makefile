.PHONY: setup-amrfinder amrfinder features test

setup-amrfinder:
	powershell -ExecutionPolicy Bypass -File scripts/setup_amrfinder_docker.ps1

amrfinder:
	powershell -ExecutionPolicy Bypass -File scripts/run_amrfinder_docker.ps1

features:
	python -m genome_firewall.build_feature_matrix --config config/pipeline.yaml

test:
	pytest -v


download-bvbrc-release:
	bash scripts/download_bvbrc_release_files.sh

audit:
	python -m genome_firewall.audit_bvbrc

prepare-cohort:
	python -m genome_firewall.prepare_cohort

pilot-download:
	python -m genome_firewall.download_fastas --mode pilot
