# Genome Firewall

Ein defensiver Hackathon-Prototyp zur Analyse von bakteriellen FASTA-Dateien auf AMR-Marker (Antimicrobial Resistance) mittels AMRFinderPlus.

> **Wichtiger Hinweis:** Die Abwesenheit eines Resistenzmarkers beweist nicht, dass ein Organismus gegenüber einem bestimmten Antibiotikum sensibel ist. Dies ist ein Forschungsprototyp und kein klinisches Entscheidungssystem!

## Architektur & Konzept

- **Python-Pipeline:** Die Datenverarbeitung erfolgt primär über Python (Nativ unter Windows/Linux).
- **AMRFinderPlus:** Wird unter Windows vollständig isoliert im gepinnten Docker-Container (`ncbi/amr:4.2.7-2026-05-15.1`) ausgeführt, um Abhängigkeitsprobleme (Bioconda unter win-64) zu vermeiden. Die Outputs werden per Mount auf dem Host-System abgelegt und dann von Python (Pandas) gelesen.
- **Feature-Matrix:** Das System extrahiert bekannte Gene und spezifische Mutationen (z.B. `S84L` in `gyrA`) für `Staphylococcus_aureus` und erzeugt eine binäre Presence/Absence Matrix.

## Windows setup with Conda and Docker

Folge diesen Schritten strikt in dieser Reihenfolge, um die Umgebung unter Windows 11 aufzusetzen.

### Voraussetzungen:
- Git
- Anaconda oder Miniconda
- Docker Desktop (läuft und verwendet Linux-Container-Modus)

### Schritte:

1. **In das Repository wechseln:**
   ```powershell
   cd path\to\biosecurity-AMR
   ```

2. **Conda-Umgebung in der Anaconda Prompt erstellen:**
   ```powershell
   conda env create -f environment.yml
   ```
   *(Hinweis: Diese Umgebung enthält bewusst **kein** `ncbi-amrfinderplus`, da dies nicht für `win-64` verfügbar ist. Wir verwenden stattdessen Docker!)*

3. **Umgebung aktivieren:**
   ```powershell
   conda activate genome-firewall
   ```

4. **Projekt als editierbares Paket installieren:**
   ```powershell
   python -m pip install -e .
   ```

5. **Docker-AMRFinderPlus einrichten:**
   ```powershell
   powershell -ExecutionPolicy Bypass -File .\scripts\setup_amrfinder_docker.ps1
   ```

6. **FASTA-Dateien bereitstellen:**
   Kopiere `.fa` / `.fna` / `.fasta` Dateien (nicht komprimiert) in folgenden Ordner:
   ```text
   data\raw\genomes
   ```

7. **AMRFinderPlus ausführen:**
   ```powershell
   powershell -ExecutionPolicy Bypass -File .\scripts\run_amrfinder_docker.ps1 -Limit 1 -FailFast
   ```
   *Oder direkt über Python:*
   ```powershell
   python -m genome_firewall.run_amrfinder --backend docker --config config/pipeline.yaml --limit 1 --fail-fast
   ```

8. **Feature-Matrix erzeugen:**
   ```powershell
   python -m genome_firewall.build_feature_matrix --config config/pipeline.yaml
   ```

9. **Tests ausführen:**
   ```powershell
   pytest -v
   ```

## Linux/WSL/Codespaces Setup

Für Linux-basierte Systeme kann AMRFinderPlus nativ installiert werden:

```bash
conda env create -f environment-linux.yml
conda activate genome-firewall-linux
python -m pip install -e .
python -m genome_firewall.run_amrfinder --backend native --config config/pipeline.yaml
```

*Hinweis: Verwende `environment-linux.yml` **nicht** unter Windows.*

## Troubleshooting

- **`conda` wird in PowerShell nicht gefunden:** Öffne stattdessen die "Anaconda Prompt" aus dem Startmenü. Verändere nicht den globalen Windows PATH.
- **Docker-Daemon läuft nicht:** Starte Docker Desktop manuell. Stelle sicher, dass die WSL2 Integration in den Docker-Einstellungen aktiv ist.
- **Keine FASTA-Dateien gefunden:** Die Pipeline sucht nach dekomprimierten Dateien. `.fa.gz` wird aktuell nicht unterstützt. Bitte vorher entpacken.
- **Python-Paket nicht installiert:** Stelle sicher, dass `python -m pip install -e .` im aktivierten Environment aufgerufen wurde, da sonst `genome_firewall` nicht gefunden wird.
- **AMRFinder Datenbank-Update (`amrfinder -u`):** Wird im Rahmen des Docker-Workflows **nicht** verwendet, um Reproduzierbarkeit (Pinned Docker Image) sicherzustellen. Ein Update erfolgt nur durch explizite Änderung des Image-Tags.
## BV-BRC data pipeline

- **Hinweis:** NCTC 8325 ist nur das Smoke-Test-Referenzgenom.
- Die finale Trainingspopulation umfasst viele unterschiedliche Staphylococcus-aureus-Isolate.
- Die Labels stammen aus laborbasierten BV-BRC-AMR-Phänotypen.
- *Computational Prediction* darf nicht als Laborlabel verwendet werden!
- Features (AMRFinderPlus) und Labels (BV-BRC) bleiben stets getrennt.
- Zunächst wird ein Audit ausgeführt: `python -m genome_firewall.audit_bvbrc` oder `make audit`
- Danach werden geeignete Antibiotika in der Config ausgewählt.
- Erst anschließend werden die zugehörigen FASTAs heruntergeladen: `python -m genome_firewall.prepare_cohort` & `python -m genome_firewall.download_fastas --mode pilot`

## Single-genome prediction

Beispiel:

`ash
python scripts/predict_genome.py patient.fna
`

- **Spezies:** Unterstützt ausschließlich *Staphylococcus aureus*.
- **Backend:** Benötigt Docker als Standard für AMRFinderPlus.
- **Vorhersagen:** Erstellt Resistenzwahrscheinlichkeiten (probability_resistant) für 3 Antibiotika (Cefoxitin, Ciprofloxacin, Erythromycin).
- **Output:** Erzeugt JSON und CSV Vorhersagen im Ausgabeverzeichnis.
- **Warnung:** Dies sind development-only Modelle und **kein klinisches Entscheidungssystem**!
