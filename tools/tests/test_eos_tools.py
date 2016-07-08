"""Test reduced equation of state tools."""
import unittest
import numpy as np
from simphony.tools import reos
from simphony.tools.maxwell_equal_area_rule import (
    liquid_vapour_density)


class MaxwellTestCase(unittest.TestCase):

    """Test application of the Maxwell's equal area rule."""

    def test_maxwell(self):
        """Density computation with several reduced equations of state."""
        ig = reos.ReducedIdealGas()
        vdw = reos.ReducedVanDerWaals()
        rk = reos.ReducedRedlichKwong()
        so = reos.ReducedSoave()
        pr = reos.ReducedPengRobinson()
        cs = reos.ReducedCarnahanStarling()

        # Note, sorted
        Trs = [0.35, 0.55, 0.75, 0.95]
        Trs_len = len(Trs)

        vdw_results = liquid_vapour_density(Trs, vdw)
        rk_results = liquid_vapour_density(Trs, rk)
        so_results = liquid_vapour_density(Trs, so)
        pr_results = liquid_vapour_density(Trs, pr)
        cs_results = liquid_vapour_density(Trs, cs)

        self.assertEqual(Trs_len, len(vdw_results))
        self.assertEqual(Trs_len, len(rk_results))
        self.assertEqual(Trs_len, len(so_results))
        self.assertEqual(Trs_len, len(pr_results))
        self.assertEqual(Trs_len, len(cs_results))

        for r in np.arange(0, Trs_len-1):
            self.assertEqual(Trs[r], vdw_results[r][0])
            self.assertEqual(Trs[r], rk_results[r][0])
            self.assertEqual(Trs[r], so_results[r][0])
            self.assertEqual(Trs[r], pr_results[r][0])
            self.assertEqual(Trs[r], cs_results[r][0])

            # Reduced molar volume must be
            # greater for the vapour phase than for the liquiq phase
            self.assertGreater(vdw_results[r][2], vdw_results[r][3])
            self.assertGreater(rk_results[r][2], rk_results[r][3])
            self.assertGreater(so_results[r][2], so_results[r][3])
            self.assertGreater(pr_results[r][2], pr_results[r][3])
            self.assertGreater(cs_results[r][2], cs_results[r][3])

            # Reduced density must be
            # less for the vapour phase than for the liquiq phase
            self.assertLess(vdw_results[r][4], vdw_results[r][5])
            self.assertLess(rk_results[r][4], rk_results[r][5])
            self.assertLess(so_results[r][4], so_results[r][5])
            self.assertLess(pr_results[r][4], pr_results[r][5])
            self.assertLess(cs_results[r][4], cs_results[r][5])

            # Reduced equilibrium pressure must be smaller
            # for lower reduced temperature (note that Trs must be sorted)
            self.assertLess(vdw_results[r][1], vdw_results[r+1][1])
            self.assertLess(rk_results[r][1], rk_results[r+1][1])
            self.assertLess(so_results[r][1], so_results[r+1][1])
            self.assertLess(pr_results[r][1], pr_results[r+1][1])
            self.assertLess(cs_results[r][1], cs_results[r+1][1])

            # Reduced molar volume for the vapour phase must be greater
            # for lower reduced temperature (note that Trs must be sorted)
            self.assertGreater(vdw_results[r][2], vdw_results[r+1][2])
            self.assertGreater(rk_results[r][2], rk_results[r+1][2])
            self.assertGreater(so_results[r][2], so_results[r+1][2])
            self.assertGreater(pr_results[r][2], pr_results[r+1][2])
            self.assertGreater(cs_results[r][2], cs_results[r+1][2])

            # Reduced molar volume for the liquid phase must be smaller
            # for lower reduced temperature (note that Trs must be sorted)
            self.assertLess(vdw_results[r][3], vdw_results[r+1][3])
            self.assertLess(rk_results[r][3], rk_results[r+1][3])
            self.assertLess(so_results[r][3], so_results[r+1][3])
            self.assertLess(pr_results[r][3], pr_results[r+1][3])
            self.assertLess(cs_results[r][3], cs_results[r+1][3])

            # Reduced density for the vapour phase must be smaller
            # for lower reduced temperature (note that Trs must be sorted)
            self.assertLess(vdw_results[r][4], vdw_results[r+1][4])
            self.assertLess(rk_results[r][4], rk_results[r+1][4])
            self.assertLess(so_results[r][4], so_results[r+1][4])
            self.assertLess(pr_results[r][4], pr_results[r+1][4])
            self.assertLess(cs_results[r][4], cs_results[r+1][4])

            # Reduced density for the liquiq phase must be greater
            # for lower reduced temperature (note that Trs must be sorted)
            self.assertGreater(vdw_results[r][5], vdw_results[r+1][5])
            self.assertGreater(rk_results[r][5], rk_results[r+1][5])
            self.assertGreater(so_results[r][5], so_results[r+1][5])
            self.assertGreater(pr_results[r][5], pr_results[r+1][5])
            self.assertGreater(cs_results[r][5], cs_results[r+1][5])

        try:
            liquid_vapour_density(Trs, ig)
        except RuntimeError:
            pass
        else:
            msg = 'Computation of densities must fail with Ideal-gas reos.'
            raise AssertionError(msg)

        try:
            liquid_vapour_density([1.01], vdw)
        except ValueError:
            pass
        else:
            raise AssertionError('Tr > 1 must fail with VdW reos.')

        try:
            liquid_vapour_density([1.01], rk)
        except ValueError:
            pass
        else:
            raise AssertionError('Tr > 1 must fail with R-K reos.')

        try:
            liquid_vapour_density([1.01], so)
        except ValueError:
            pass
        else:
            raise AssertionError('Tr > 1 must fail with Soave reos.')

        try:
            liquid_vapour_density([1.01], pr)
        except ValueError:
            pass
        else:
            raise AssertionError('Tr > 1 must fail with P-R reos.')

        try:
            liquid_vapour_density([1.01], cs)
        except ValueError:
            pass
        else:
            raise AssertionError('Tr > 1 must fail with C-S reos.')


if __name__ == '__main__':
    unittest.main()
