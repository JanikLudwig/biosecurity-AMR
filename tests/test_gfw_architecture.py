"""Focused invariants for the M1/M3 and M2/M4 branch boundary."""

import unittest

import numpy as np
import pandas as pd

from gfw.decide import NO_CALL, WORK, decide_drug
from gfw.engine import Engine
from gfw.panel import DrugEntry
from gfw.predict import DrugModel


class FixedProbabilityModel:
    def __init__(self, probability):
        self.probability = probability
        self.vectors = []

    def predict_proba(self, vector):
        self.vectors.append(vector.copy())
        return np.array([[1.0 - self.probability, self.probability]])


def entry():
    return DrugEntry(
        drug="ciprofloxacin", drug_display="Ciprofloxacin", tier="A",
        n_resistant=100, n_susceptible=100, n_total=200, n_groups=20,
        drug_class="fluoroquinolone", mechanism="", target_kind="protein",
        target_genes=["gyrA"], modelable=True)


class TestBranchBoundary(unittest.TestCase):
    def setUp(self):
        self.entry = entry()

    def test_m3_probability_is_calculated_without_m2(self):
        estimator = FixedProbabilityModel(0.12)
        model = DrugModel("ciprofloxacin", "A", ["mecA", "gyrA_S84L"], estimator,
                          {}, 10, 5, "none")
        engine = Engine(models={"ciprofloxacin": model}, panel=[self.entry],
                        references={}, features=pd.DataFrame(
                            {"mecA": [1], "gyrA_S84L": [0]}, index=["sample"]))
        result = engine.predict_m1("sample")
        self.assertEqual(result.p_resistant, {"ciprofloxacin": 0.12})
        np.testing.assert_array_equal(estimator.vectors[0], [[1, 0]])

    def test_low_probability_with_absent_target_is_no_call(self):
        decision = decide_drug(self.entry, 0.12, "absent", [], [])
        self.assertEqual(decision.call, NO_CALL)
        self.assertEqual(decision.no_call_reason, "drug_target_absent")

    def test_low_probability_with_present_target_can_work(self):
        decision = decide_drug(self.entry, 0.12, "present", [], [{"gene": "gyrA"}])
        self.assertEqual(decision.call, WORK)

    def test_m2_cannot_alter_m3_feature_vector_or_probability(self):
        estimator = FixedProbabilityModel(0.12)
        model = DrugModel("ciprofloxacin", "A", ["mecA", "gyrA_S84L"], estimator,
                          {}, 10, 5, "none")
        m1_features = {"mecA": 1, "unmodelled_symbol": 1}
        probability_before = model.predict_proba(m1_features)
        # This is representative M2 evidence: it is deliberately not accepted by M3.
        m2_evidence = {"target_status": "absent", "detected": [], "missing": ["gyrA"]}
        decision = decide_drug(self.entry, probability_before, m2_evidence["target_status"],
                               list(m1_features), m2_evidence["detected"])
        probability_after = model.predict_proba(m1_features)
        self.assertEqual(decision.call, NO_CALL)
        self.assertEqual(probability_before, probability_after)
        np.testing.assert_array_equal(estimator.vectors[0], estimator.vectors[1])
        np.testing.assert_array_equal(estimator.vectors[0], [[1, 0]])


if __name__ == "__main__":
    unittest.main()
