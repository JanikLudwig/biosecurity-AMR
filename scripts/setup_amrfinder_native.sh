#!/usr/bin/env bash
set -e

mkdir -p reports/

echo "Prüfe AMRFinderPlus Installation..."

if ! command -v amrfinder &> /dev/null; then
    echo "Fehler: amrfinder wurde nicht gefunden. Bitte Conda-Umgebung aktivieren (conda activate amr_pipeline)."
    exit 1
fi

echo "Aktualisiere AMRFinderPlus Datenbank (-u)..."
amrfinder -u

echo "Speichere Version in reports/amrfinder_version.txt..."
amrfinder -V > reports/amrfinder_version.txt

echo "Speichere Liste der verfügbaren Organismen..."
amrfinder -l > reports/amrfinder_organisms.txt

echo "Prüfe auf Staphylococcus_aureus Organismus-Unterstützung..."
if grep -q "Staphylococcus_aureus" reports/amrfinder_organisms.txt; then
    echo "Staphylococcus_aureus ist als Option verfügbar."
else
    echo "Fehler: Staphylococcus_aureus ist nicht in der Liste der unterstützten Organismen."
    exit 1
fi

echo "AMRFinderPlus Setup erfolgreich abgeschlossen."
