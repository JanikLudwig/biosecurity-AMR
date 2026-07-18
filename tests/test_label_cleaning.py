import pytest
import pandas as pd
from pathlib import Path

# Wir können einen Dummy Test schreiben, der die Logik des Resolvings aus prepare_cohort.py nachstellt
def test_label_conflict_resolution():
    # Simulation der Daten
    data = [
        {"genome_id": "1", "antibiotic": "ciprofloxacin", "resistant_phenotype": "Resistant"},
        {"genome_id": "1", "antibiotic": "ciprofloxacin", "resistant_phenotype": "Resistant"},  # Duplikat
        {"genome_id": "2", "antibiotic": "oxacillin", "resistant_phenotype": "Resistant"},
        {"genome_id": "2", "antibiotic": "oxacillin", "resistant_phenotype": "Susceptible"} # Konflikt
    ]
    df = pd.DataFrame(data)
    
    clean_labels = []
    excluded = []
    
    for (gid, ab), group in df.groupby(["genome_id", "antibiotic"]):
        phenotypes = group["resistant_phenotype"].dropna().unique()
        if len(phenotypes) > 1:
            excluded.append(gid)
        else:
            repr_row = group.iloc[0].copy()
            clean_labels.append(repr_row)
            
    assert "1" not in excluded
    assert "2" in excluded
    assert len(clean_labels) == 1
    assert clean_labels[0]["genome_id"] == "1"

def test_invalid_labels():
    # Nur erlaubte Labels
    phenotypes = ["Resistant", "Susceptible", "Intermediate", "Unknown"]
    valid = [p for p in phenotypes if p in ["Resistant", "Susceptible", "Intermediate"]]
    assert "Unknown" not in valid
    assert "Resistant" in valid
