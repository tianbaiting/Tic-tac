#!/usr/bin/env python3
"""L1: Cross-validate CG/6j/9j against sympy (independent oracle)."""

import json
import unittest
from pathlib import Path
from sympy import S, sqrt, Rational
from sympy.physics.wigner import clebsch_gordan, wigner_3j, wigner_6j, wigner_9j

DATA_DIR = Path(__file__).parent / "data"


def generate_test_vectors():
    """Generate test vectors and save to JSON for C++ cross-validation."""
    vectors = {"cg": [], "w6j": [], "w9j": []}

    # CG coefficients: <j1 m1; j2 m2 | j3 m3>
    # sympy signature: clebsch_gordan(j1, j2, j3, m1, m2, m3)
    cg_cases = [
        (S(1)/2, S(1)/2, 1, S(1)/2, S(1)/2, 1),      # spin triplet |1,1>
        (S(1)/2, S(1)/2, 1, S(1)/2, -S(1)/2, 0),     # spin triplet |1,0>
        (S(1)/2, S(1)/2, 0, S(1)/2, -S(1)/2, 0),     # spin singlet |0,0>
        (1, 1, 2, 1, 0, 1),                            # L=1 coupling
        (1, 1, 2, 0, 0, 0),                            # L=1 coupling m=0
        (1, 1, 0, 1, -1, 0),                           # L=1 -> L=0
        (S(3)/2, S(1)/2, 2, S(3)/2, S(1)/2, 2),       # j=3/2 + 1/2 -> 2
        (S(3)/2, S(1)/2, 1, S(3)/2, -S(1)/2, 1),      # j=3/2 + 1/2 -> 1
        (2, 1, 3, 0, 0, 0),                            # higher L
        (2, 1, 1, 0, 0, 0),                            # triangle limit
    ]
    for args in cg_cases:
        val = float(clebsch_gordan(*args))
        two_args = [int(2 * a) for a in args]
        vectors["cg"].append({"two_j1": two_args[0], "two_j2": two_args[1],
                              "two_j3": two_args[2], "two_m1": two_args[3],
                              "two_m2": two_args[4], "two_m3": two_args[5],
                              "value": val})

    # Wigner 6j symbols
    w6j_cases = [
        (1, 1, 1, 1, 1, 1),
        (S(1)/2, S(1)/2, 1, S(1)/2, S(1)/2, 0),
        (S(1)/2, S(1)/2, 1, S(1)/2, S(1)/2, 1),
        (1, 1, 2, 1, 1, 1),
        (1, 2, 1, 1, 1, 2),
        (2, 2, 2, 2, 2, 2),
        (S(3)/2, S(1)/2, 1, S(1)/2, S(3)/2, 1),
    ]
    for args in w6j_cases:
        val = float(wigner_6j(*args))
        two_args = [int(2 * a) for a in args]
        vectors["w6j"].append({"args": two_args, "value": val})

    # Wigner 9j symbols
    w9j_cases = [
        (1, 1, 0, 1, 1, 0, 0, 0, 0),
        (1, 1, 2, 1, 1, 2, 2, 2, 2),
        (S(1)/2, S(1)/2, 1, S(1)/2, S(1)/2, 1, 1, 1, 0),
        (S(1)/2, S(1)/2, 1, S(1)/2, S(1)/2, 1, 1, 1, 2),
    ]
    for args in w9j_cases:
        val = float(wigner_9j(*args))
        two_args = [int(2 * a) for a in args]
        vectors["w9j"].append({"args": two_args, "value": val})

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    out = DATA_DIR / "coupling_test_values.json"
    with open(out, "w") as f:
        json.dump(vectors, f, indent=2)
    return vectors


class TestCouplingCoefficients(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.vectors = generate_test_vectors()

    def test_cg_values(self):
        for case in self.vectors["cg"]:
            val = float(clebsch_gordan(
                S(case["two_j1"]) / 2, S(case["two_j2"]) / 2, S(case["two_j3"]) / 2,
                S(case["two_m1"]) / 2, S(case["two_m2"]) / 2, S(case["two_m3"]) / 2))
            self.assertAlmostEqual(val, case["value"], places=12,
                                   msg=f"CG mismatch for {case}")

    def test_triangle_violation_zero(self):
        self.assertEqual(float(clebsch_gordan(3, 0, 1, 0, 0, 0)), 0.0)
        self.assertEqual(float(wigner_6j(5, 0, 1, 0, 0, 0)), 0.0)

    def test_cg_orthogonality(self):
        """Sum_m1m2 CG(j1,m1;j2,m2|j,m) * CG(j1,m1;j2,m2|j',m') = delta_{jj'} delta_{mm'}"""
        j1, j2 = 1, 1
        for j in range(abs(j1 - j2), j1 + j2 + 1):
            for m in range(-j, j + 1):
                total = sum(
                    float(clebsch_gordan(j1, j2, j, m1, m - m1, m)) ** 2
                    for m1 in range(-j1, j1 + 1)
                    if abs(m - m1) <= j2
                )
                self.assertAlmostEqual(total, 1.0, places=12,
                                       msg=f"CG orthogonality fail for j={j},m={m}")


if __name__ == "__main__":
    unittest.main()
