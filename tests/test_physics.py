"""Physics invariants and independent manufactured-solution checks."""
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pn_junction import EPS_SI, Q, Junction, analytical, ideal_iv
from poisson_numeric import compare, depletion_charge, solve_poisson


class PhysicsTests(unittest.TestCase):
    def test_default_and_charge_neutrality(self):
        j = Junction()
        self.assertAlmostEqual(j.vbi, 0.77384358, places=7)
        self.assertAlmostEqual(j.na * j.xp / (j.nd * j.xn), 1.0)
        self.assertAlmostEqual(j.xp / j.xn, j.nd / j.na)

    def test_boundary_values_and_equilibrium_bands(self):
        j = Junction()
        x = np.array([-2*j.width, -j.xp, -1e-20, 0, 1e-20, j.xn, 2*j.width])
        d = analytical(j, x)
        np.testing.assert_allclose(d["psi"][[0, 1, -2, -1]], [0, 0, j.vbi, j.vbi], atol=1e-14)
        np.testing.assert_allclose(d["field"][[0, 1, -2, -1]], 0, atol=1e-10)
        self.assertLess(np.max(np.abs(np.diff(d["psi"][2:5]))), 1e-12)
        self.assertTrue(np.all(d["field"] <= 0))
        np.testing.assert_allclose(d["Ef"], 0)
        self.assertAlmostEqual(d["Ei"][0], j.vt * np.log(j.na / j.ni))
        self.assertAlmostEqual(d["Ei"][-1], -j.vt * np.log(j.nd / j.ni))

    def test_poisson_quadratic_manufactured_solution(self):
        x = np.linspace(0, 1, 101)
        # psi = 2 + 3*x + x*(1-x); independent of junction formulas.
        rho = np.full_like(x, 2 * EPS_SI)
        actual = solve_poisson(x, rho, 2, 5)
        np.testing.assert_allclose(actual, 2 + 3*x + x*(1-x), atol=1e-11)

    def test_second_order_smooth_manufactured_solution(self):
        errors = []
        for n in (51, 101, 201):
            x = np.linspace(0, 1, n)
            exact = np.sin(np.pi*x) + 2*x
            rho = EPS_SI * np.pi**2 * np.sin(np.pi*x)
            errors.append(np.max(np.abs(solve_poisson(x, rho, 0, 2) - exact)))
        self.assertTrue(3.9 < errors[0] / errors[1] < 4.1)
        self.assertTrue(3.9 < errors[1] / errors[2] < 4.1)

    def test_asymmetric_poisson_comparisons(self):
        for na, nd in ((1e16, 1e17), (1e17, 1e16), (1e16, 1e16)):
            j = Junction(na, nd)
            coarse, fine = compare(j, 201), compare(j, 1601)
            m = fine["metrics"]
            self.assertLess(m["max_potential_error_V"], 2e-5)
            self.assertLess(m["max_potential_error_V"], coarse["metrics"]["max_potential_error_V"])
            self.assertLess(m["scaled_discrete_residual"], 1e-7)
            self.assertLess(abs(m["net_sheet_charge_C_m2"]), 1e-14)

    def test_ideal_diode(self):
        j = Junction()
        current, isat = ideal_iv(j, np.array([-.5, 0, .1]))
        self.assertAlmostEqual(current[0] / isat, -1.0, places=7)
        self.assertEqual(current[1], 0)
        self.assertGreater(current[2], 0)
        doubled, _ = ideal_iv(j, np.array([.1]), area_m2=2e-8)
        self.assertAlmostEqual(doubled[0] / current[2], 2)

    def test_parameter_trends_and_validation(self):
        self.assertLess(Junction(temperature=400).vbi, Junction(temperature=250).vbi)
        self.assertGreater(Junction(na_cm3=1e17).vbi, Junction().vbi)
        self.assertLess(Junction(na_cm3=1e17).width, Junction().width)
        for kwargs in ({"na_cm3": -1}, {"nd_cm3": float("nan")}, {"temperature": 0}):
            with self.assertRaises(ValueError):
                Junction(**kwargs)
        with self.assertRaises(ValueError):
            solve_poisson(np.array([0, 1, 2, 4, 5.]), np.zeros(5), 0, 1)


if __name__ == "__main__":
    unittest.main()
