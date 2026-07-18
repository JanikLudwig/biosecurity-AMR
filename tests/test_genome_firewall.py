"""End-to-end and unit tests for Genome Firewall v0 (stdlib unittest — no deps)."""

import os
import tempfile
import unittest

from genome_firewall.annotate import AmrHit, AnnotationResult, annotate, parse_amr_table
from genome_firewall.predict import (NO_CALL, RESISTANT, SUSCEPTIBLE,
                                     EV_KNOWN, EV_NONE, predict_sample, _noisy_or)
from genome_firewall import knowledge as kb

_EX = os.path.join(os.path.dirname(__file__), "..", "genome_firewall", "examples")


def _calls(sample):
    return {p.drug_id: p.call for p in sample.predictions}


def _predict_from_tsv(name, species="Escherichia coli", source="amrfinderplus"):
    ann = annotate(backend="tsv", tsv_path=os.path.join(_EX, name), tsv_source=source)
    return predict_sample(ann, species=species)


class TestResistantProfile(unittest.TestCase):
    def test_all_drugs_flagged_resistant(self):
        sample = _predict_from_tsv("ecoli_resistant_amrfinder.tsv")
        calls = _calls(sample)
        for drug in ("ampicillin", "ceftriaxone", "ciprofloxacin",
                     "gentamicin", "cotrimoxazole"):
            self.assertEqual(calls[drug], RESISTANT, f"{drug} should be resistant")

    def test_resistant_evidence_and_markers(self):
        sample = _predict_from_tsv("ecoli_resistant_amrfinder.tsv")
        cef = next(p for p in sample.predictions if p.drug_id == "ceftriaxone")
        self.assertEqual(cef.evidence_category, EV_KNOWN)
        self.assertTrue(any("CTX-M" in m["gene"] for m in cef.supporting_markers))
        self.assertGreaterEqual(cef.confidence, 0.6)


class TestSusceptibleProfile(unittest.TestCase):
    def test_no_markers_gives_susceptible(self):
        sample = _predict_from_tsv("ecoli_susceptible_amrfinder.tsv")
        for p in sample.predictions:
            self.assertEqual(p.call, SUSCEPTIBLE)
            self.assertEqual(p.evidence_category, EV_NONE)
            self.assertEqual(p.target_status, "present")
            # Absence of evidence must be bounded, never high-confidence.
            self.assertLessEqual(p.confidence, 0.8)


class TestWeakProfile(unittest.TestCase):
    def test_weak_markers_yield_no_call(self):
        sample = _predict_from_tsv("ecoli_weak_amrfinder.tsv")
        calls = _calls(sample)
        self.assertEqual(calls["ciprofloxacin"], NO_CALL)   # oqxA only, low-level
        self.assertEqual(calls["cotrimoxazole"], NO_CALL)    # sul only, no dfr
        self.assertEqual(calls["ampicillin"], SUSCEPTIBLE)


class TestEsblDistinction(unittest.TestCase):
    """A plain penicillinase must NOT be reported as ceftriaxone resistance."""

    def test_tem1_hits_ampicillin_not_ceftriaxone(self):
        hit = AmrHit(gene="blaTEM-1", drug_class="BETA-LACTAM",
                     subclass="BETA-LACTAM", identity=100.0, coverage=100.0,
                     source="amrfinderplus")
        ann = AnnotationResult([hit], backend="tsv:amrfinderplus",
                               screening_completeness=0.65)
        sample = predict_sample(ann, species="Escherichia coli")
        calls = _calls(sample)
        self.assertEqual(calls["ampicillin"], RESISTANT)
        self.assertEqual(calls["ceftriaxone"], SUSCEPTIBLE)


class TestSampleGates(unittest.TestCase):
    def test_out_of_scope_species_all_no_call(self):
        ann = annotate(backend="tsv",
                       tsv_path=os.path.join(_EX, "ecoli_resistant_amrfinder.tsv"),
                       tsv_source="amrfinderplus")
        sample = predict_sample(ann, species="Staphylococcus aureus")
        self.assertFalse(sample.species_supported)
        self.assertTrue(all(p.call == NO_CALL for p in sample.predictions))

    def test_failed_qc_all_no_call(self):
        ann = annotate(backend="tsv",
                       tsv_path=os.path.join(_EX, "ecoli_resistant_amrfinder.tsv"),
                       tsv_source="amrfinderplus")
        bad_qc = {"passed": False, "flags": ["assembly_too_short_for_bacterial_genome"]}
        sample = predict_sample(ann, species="Escherichia coli", qc=bad_qc)
        self.assertTrue(all(p.call == NO_CALL for p in sample.predictions))


class TestScreeningCompleteness(unittest.TestCase):
    """cAMRah's broader screen justifies higher 'likely to work' confidence."""

    def test_camrah_beats_amrfinder_on_absence_confidence(self):
        s_amr = _predict_from_tsv("ecoli_susceptible_amrfinder.tsv", source="amrfinderplus")
        s_cam = _predict_from_tsv("ecoli_susceptible_amrfinder.tsv", source="camrah")
        amp_amr = next(p for p in s_amr.predictions if p.drug_id == "ampicillin")
        amp_cam = next(p for p in s_cam.predictions if p.drug_id == "ampicillin")
        self.assertGreater(amp_cam.confidence, amp_amr.confidence)


class TestUnits(unittest.TestCase):
    def test_noisy_or(self):
        self.assertAlmostEqual(_noisy_or([]), 0.0)
        self.assertAlmostEqual(_noisy_or([0.5]), 0.5)
        self.assertAlmostEqual(_noisy_or([0.5, 0.5]), 0.75)
        self.assertGreater(_noisy_or([0.9, 0.78]), 0.9)

    def test_parse_table_keeps_only_amr_rows(self):
        hits = parse_amr_table(os.path.join(_EX, "ecoli_resistant_amrfinder.tsv"))
        self.assertEqual(len(hits), 7)
        self.assertTrue(any(h.is_point_mutation for h in hits))  # gyrA_S83L

    def test_species_support(self):
        self.assertTrue(kb.species_supported("Escherichia coli"))
        self.assertTrue(kb.species_supported(None))
        self.assertFalse(kb.species_supported("Klebsiella pneumoniae"))


class TestDedup(unittest.TestCase):
    def test_identical_genomes_cluster_together(self):
        from genome_firewall.dedup import cluster_genomes
        seq = ">c1\n" + "ACGTACGTACGTACGTGGGGCCCCAAAATTTT" * 40 + "\n"
        other = ">c1\n" + "TTTTGGGGCCCCAAAA" * 90 + "\n"
        with tempfile.TemporaryDirectory() as d:
            a, b, c = (os.path.join(d, n) for n in ("a.fasta", "b.fasta", "c.fasta"))
            for path, data in ((a, seq), (b, seq), (c, other)):  # a == b, c differs
                with open(path, "w") as fh:
                    fh.write(data)
            clusters = cluster_genomes([a, b, c], k=11, num_hashes=100, threshold=0.9)
        # a and b should collapse to one cluster; c stays separate -> 2 clusters.
        self.assertEqual(len(clusters), 2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
