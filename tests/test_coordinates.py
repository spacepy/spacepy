import spacepy_testing
import glob
import os
import unittest
from numpy import array
import numpy as np
from spacepy import ctrans
import spacepy.coordinates as spc
from spacepy.time import Ticktock
try:
    import spacepy.irbempy as ib
except ImportError:
    pass #tests will fail, but won't bring down the entire suite
try:
    import astropy.coordinates
    HAVE_ASTROPY = True
except ImportError:
    HAVE_ASTROPY = False
import spacepy.toolbox as tb

__all__ = ['coordsTest', 'coordsTestIrbem', 'QuaternionFunctionTests']


class coordsTest(unittest.TestCase):
    def setUp(self):
        self.cvals = spc.Coords([[1, 2, 4], [1, 2, 2]], 'GEO', 'car', use_irbem=False)

    def tearDown(self):
        pass

    def test_coords(self):
        """Coords should create and do simple conversions"""
        np.testing.assert_equal([1,1], self.cvals.x)
        np.testing.assert_equal([2,2], self.cvals.y)
        np.testing.assert_equal([4,2], self.cvals.z)
        self.cvals.ticks = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO') # add ticktock
        newcoord = self.cvals.convert('GSM', 'sph')

    def test_array_input_1D(self):
        """Coords should build correctly and convert with array input"""
        cvals = spc.Coords(np.array([1, 2, 4]), 'GEO', 'car', use_irbem=False)
        out = cvals.convert('GEO', 'sph')

    def test_array_input_2D(self):
        """Coords should build correctly and convert with array input"""
        cvals = spc.Coords(np.array([[1, 2, 4], [1, 2, 2]]), 'GEO', 'car', use_irbem=False)
        out = cvals.convert('GEO', 'sph')

    def test_bad_position_fails_1D(self):
        """Positions nnot supplied as valid 3-vectors should fail"""
        self.assertRaises(ValueError, spc.Coords, [[1, 2],[1, 2]], 'GEO', 'car', use_irbem=False)

    def test_bad_position_fails_2D(self):
        """Positions nnot supplied as valid 3-vectors should fail"""
        self.assertRaises(ValueError, spc.Coords, np.array([1, 2]), 'GEO', 'car', use_irbem=False)

    def test_bad_position_fails_scalar(self):
        """Positions nnot supplied as valid 3-vectors should fail"""
        self.assertRaises(ValueError, spc.Coords, 1, 'GEO', 'car', use_irbem=False)

    def test_append(self):
        """Test append functionality"""
        c2 = spc.Coords([[6, 7, 8], [9, 10, 11]], 'GEO', 'car', use_irbem=False)
        actual = self.cvals.append(c2)
        expected = [[1, 2, 4], [1, 2, 2], [6, 7, 8], [9, 10, 11]]
        np.testing.assert_equal(expected, actual.data.tolist())

    def test_slice(self):
        """Test slice functionality"""
        expected = spc.Coords([1, 2, 4], 'GEO', 'car', use_irbem=False)
        np.testing.assert_equal(expected.data, self.cvals[0].data)

    def test_slice_with_ticks(self):
        """Test slice functionality with ticks attribute"""
        self.cvals.ticks = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO')
        expected = spc.Coords([1, 2, 4], 'GEO', 'car', use_irbem=False)
        np.testing.assert_equal(expected.data, self.cvals[0].data)

    def test_roundtrip_GEO_ECIMOD(self):
        """Roundtrip should yield input as answer"""
        self.cvals.ticks = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO')
        expected = spc.Coords(self.cvals.data, 'GEO', 'car', use_irbem=False)
        stage1 = self.cvals.convert('ECIMOD', 'car')
        got = stage1.convert('GEO', 'car')
        np.testing.assert_allclose(got.data, expected.data)
        np.testing.assert_equal(got.dtype, 'GEO')

    def test_IRBEMname1(self):
        """CTrans-based conversion should respect IRBEM-style names on input"""
        self.cvals.ticks = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO')
        expected = spc.Coords(self.cvals.data, 'GEO', 'car', use_irbem=False)
        stage1 = self.cvals.convert('GEI', 'car')
        got = stage1.convert('GEO', 'car')
        np.testing.assert_allclose(got.data, expected.data)
        np.testing.assert_equal(got.dtype, 'GEO')

    def test_IRBEMname2(self):
        """CTrans-based conversion should respect IRBEM-style names on input"""
        cvals = spc.Coords(self.cvals.data, 'GEI', 'car', use_irbem=False)
        cvals.ticks = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO')
        stage1 = cvals.convert('GEO', 'car')
        got = stage1.convert('GEI', 'car')
        np.testing.assert_allclose(got.data, self.cvals.data)
        np.testing.assert_equal(got.dtype, spc.SYS_EQUIV['GEI'])

    def test_Coords_same_as_CTrans_1D(self):
        """Make sure 1D position array gives correct answer"""
        tt = Ticktock('2002-02-02T12:00:00', 'ISO')
        pos = [1, 2, 3]
        ctobj = ctrans.CTrans(tt)
        ctobj.calcCoreTransforms()
        expected = np.atleast_2d(ctobj.convert(pos, 'GEO', 'GSE'))
        ccobj = spc.Coords(pos, 'GEO', 'car', ticks=tt, use_irbem=False)
        got = ccobj.convert('GSE', 'car')
        np.testing.assert_allclose(got.data, expected)

    def test_Coords_same_as_CTrans_2D(self):
        """Make sure 2D position array gives correct answer"""
        tval = '2002-02-02T12:00:00'
        tt = Ticktock([tval]*3, 'ISO')
        pos = np.array([[1, 2, 3], [1, 2, 3], [1, 2, 3]])
        ctobj = ctrans.CTrans(tt[0])
        ctobj.calcCoreTransforms()
        expected = ctobj.convert(pos, 'GEO', 'ECI2000')
        ccobj = spc.Coords(pos, 'GEO', 'car', ticks=tt, use_irbem=False)
        got = ccobj.convert('ECI2000', 'car')
        np.testing.assert_allclose(got.data, expected)

    def test_spherical_return_GEO(self):
        """GEO should return correct spherical when requested (no conversion)"""
        tt = Ticktock(['2002-02-02T12:00:00'], 'ISO')
        pos = [3, 0, 3]
        expected = [4.242640687119285, 45, 0]
        ccobj = spc.Coords(pos, 'GEO', 'car', ticks=tt, use_irbem=False)
        got = ccobj.convert('GEO', 'sph')
        np.testing.assert_allclose(got.data, np.atleast_2d(expected))

    def test_spherical_return_GEO_from_ECIMOD(self):
        """GEO should return correct spherical when requested (converted)"""
        tt = Ticktock(2459213.5, 'JD')
        pos = [-0.46409501,  2.96391738,  2.99996826]
        expected = [4.242640687119285, 45, 0]
        ccobj = spc.Coords(pos, 'ECIMOD', 'car', ticks=tt, use_irbem=False)
        got = ccobj.convert('GEO', 'sph')
        np.testing.assert_allclose(got.data, np.atleast_2d(expected), rtol=0, atol=1e-7)

    def test_spherical_input_from_GEO(self):
        """Coords should return correct answer when given spherical (converted)"""
        tt = Ticktock(2459213.5, 'JD')
        pos = [4.242640687119285, 45, 0]
        expected = [-0.46409501,  2.96391738,  2.99996826]
        ccobj = spc.Coords(pos, 'GEO', 'sph', ticks=tt, use_irbem=False)
        got = ccobj.convert('ECIMOD', 'car')
        np.testing.assert_allclose(got.data, np.atleast_2d(expected), rtol=0, atol=1e-7)

    def test_geodetic_from_GEO_spherical(self):
        """Coords should return correct geodetic (converted)"""
        tt = Ticktock(2459213.5, 'JD')
        pos = [4.242640687119285, 45, 0]  # [3, 0, 3] Cartesian
        expected = np.array([20692.69835948, 45.04527886, 0])
        ccobj = spc.Coords(pos, 'GEO', 'sph', ticks=tt, use_irbem=False)
        got = ccobj.convert('GDZ', 'sph')
        np.testing.assert_allclose(got.data, np.atleast_2d(expected), rtol=0, atol=1e-7)

    def test_GDZ_in_kilometers(self):
        """Explicitly set units should be respected"""
        tt = Ticktock(2459218.5, 'JD')
        pos = np.array([2367.83158, -5981.75882, -4263.24591])
        cc_km = spc.Coords(pos, 'GEI', 'car', ticks=tt, units=['km', 'km', 'km'],
                           use_irbem=False)
        got = cc_km.convert('GDZ', 'sph')
        cc_Re = spc.Coords(pos/ctrans.WGS84['A'], 'GEI', 'car', ticks=tt,
                           units=['Re', 'Re', 'Re'], use_irbem=False)
        expected = cc_Re.convert('GDZ', 'sph')
        # same valued output in km regardless of units of input
        np.testing.assert_approx_equal(expected.lati, got.lati, significant=6)
        np.testing.assert_approx_equal(expected.long, got.long, significant=6)
        np.testing.assert_approx_equal(expected.radi, got.radi, significant=6)
        self.assertEqual(got.units[0], 'km')
        self.assertEqual(expected.units[0], 'km')

    def test_GDZ_to_other_returns_in_Re(self):
        """Docs say everything is Re except GDZ, so make sure"""
        tt = Ticktock(2459218.5, 'JD')
        test_alt = 100  # km
        pos = np.array([test_alt, 0, 90])  # km, deg, deg
        cc_km = spc.Coords(pos, 'GDZ', 'sph', ticks=tt, units=['km', 'km', 'km'],
                           use_irbem=False)
        got = cc_km.convert('GEO', 'sph')
        expected_rad = ctrans.WGS84['A'] + test_alt
        # same valued output in Re regardless of units of input
        np.testing.assert_approx_equal(expected_rad, got.radi, significant=6)
        self.assertEqual(got.units[0], 'Re')

    def test_GDZ_cartesian_raises_convert(self):
        """Geodetic coordinates shouldn't be expressed in Cartesian"""
        self.cvals.ticks = Ticktock([2459218.5]*2, 'JD')
        self.assertRaises(ValueError, self.cvals.convert, 'GDZ', 'car')

    def test_GDZ_cartesian_raises_creation(self):
        """Geodetic coordinates shouldn't be expressed in Cartesian"""
        tt = Ticktock(2459218.5, 'JD')
        pos = np.array([2367.83158, -5981.75882, -4263.24591])
        self.assertRaises(ValueError, spc.Coords, pos, 'GDZ', 'car', ticks=tt, use_irbem=False)


class coordsTestIrbem(unittest.TestCase):
    def setUp(self):
        #super(tFunctionTests, self).setUp()
        try:
            self.cvals = spc.Coords([[1, 2, 4],[1, 2, 2]], 'GEO', 'car', use_irbem=True)
        except ImportError:
            pass #tests will fail, but won't bring down the entire suite

    def tearDown(self):
        #super(tFunctionTests, self).tearDown()
        pass

    def test_coords(self):
        """Coords should create and do simple conversions"""
        np.testing.assert_equal([1, 1], self.cvals.x)
        np.testing.assert_equal([2, 2], self.cvals.y)
        np.testing.assert_equal([4, 2], self.cvals.z)
        self.cvals.ticks = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO') # add ticktock
        newcoord = self.cvals.convert('GSM', 'sph')

    def test_append(self):
        """Test append functionality"""
        c2 = spc.Coords([[6, 7, 8], [9, 10, 11]], 'GEO', 'car', use_irbem=True)
        actual = self.cvals.append(c2)
        expected = [[1, 2, 4], [1, 2, 2], [6, 7, 8], [9, 10, 11]]
        np.testing.assert_equal(expected, actual.data.tolist())

    def test_slice(self):
        """Test slice functionality"""
        expected = spc.Coords([1, 2, 4], 'GEO', 'car', use_irbem=True)
        np.testing.assert_equal(expected.data, self.cvals[0].data)

    def test_slice_with_ticks(self):
        """Test slice functionality with ticks attribute"""
        self.cvals.ticks = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO')
        expected = spc.Coords([1, 2, 4], 'GEO', 'car', use_irbem=True)
        np.testing.assert_equal(expected.data, self.cvals[0].data)

    @unittest.skipUnless(HAVE_ASTROPY, 'requires Astropy')
    def test_to_skycoord_without_ticks(self):
        with self.assertRaises(ValueError) as cm:
            sc = self.cvals.to_skycoord()

    @unittest.skipUnless(HAVE_ASTROPY, 'requires Astropy')
    def test_to_skycoord_with_ticks(self):
        self.cvals.ticks = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO') # add ticktock

        sc = self.cvals.to_skycoord()

        assert isinstance(sc, astropy.coordinates.SkyCoord)
        assert sc.frame.name == 'itrs'

        # Check that the data was loaded correctly
        sc_data = sc.cartesian.xyz.to('m').value.T
        np.testing.assert_allclose(sc_data, self.cvals.data * self.cvals.Re, rtol=1e-10)

        # Check that the time was loaded correctly
        np.testing.assert_allclose((sc.obstime - self.cvals.ticks.APT).to('s').value, 0)

    @unittest.skipUnless(HAVE_ASTROPY, 'requires Astropy')
    def test_to_skycoord_with_ticks_and_conversion(self):
        self.cvals.ticks = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO') # add ticktock

        # Convert to a frame other than GEO before calling to_skycoord()
        non_geo = self.cvals.convert('GSE', 'sph')

        sc = non_geo.to_skycoord()

        assert isinstance(sc, astropy.coordinates.SkyCoord)
        assert sc.frame.name == 'itrs'

        # Check that the data converts back correctly
        sc_data = sc.cartesian.xyz.to('m').value.T
        np.testing.assert_allclose(sc_data, self.cvals.data * self.cvals.Re, rtol=1e-10)

        # Check that the time was loaded correctly
        np.testing.assert_allclose((sc.obstime - self.cvals.ticks.APT).to('s').value, 0)

    @unittest.skipUnless(HAVE_ASTROPY, 'requires Astropy')
    def test_from_skycoord(self):
        sc = astropy.coordinates.SkyCoord(x=[1, 2], y=[4, 5], z=[7, 8],
                                          unit='Mm', frame='itrs',
                                          obstime=['2001-02-03T04:00:00', '2005-06-07T08:00:00'])

        coords = spc.Coords.from_skycoord(sc)

        # Check that the data was loaded correctly
        sc_data = sc.cartesian.xyz.to('m').value.T
        np.testing.assert_allclose(coords.data * coords.Re, sc_data, rtol=1e-10)

        # Check that the time was loaded correctly
        np.testing.assert_allclose((sc.obstime - coords.ticks.APT).to('s').value, 0)

    @unittest.skipUnless(HAVE_ASTROPY, 'requires Astropy')
    def test_from_skycoord_with_conversion(self):
        sc = astropy.coordinates.SkyCoord(x=[1, 2], y=[4, 5], z=[7, 8],
                                          unit='Mm', frame='itrs',
                                          obstime=['2001-02-03T04:00:00', '2005-06-07T08:00:00'])

        # Convert to a frame other than GEO before calling from_skycoord()
        coords = spc.Coords.from_skycoord(sc.gcrs)

        # Check that the data was loaded correctly
        sc_data = sc.cartesian.xyz.to('m').value.T
        np.testing.assert_allclose(coords.data * coords.Re, sc_data, rtol=1e-10)

        # Check that the time was loaded correctly
        np.testing.assert_allclose((sc.obstime - coords.ticks.APT).to('s').value, 0)

    def test_roundtrip_GEO_GEI(self):
        """Roundtrip should yield input as answer"""
        self.cvals.ticks = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO')
        expected = spc.Coords(self.cvals.data, 'GEO', 'car', use_irbem=True)
        stage1 = self.cvals.convert('GEI', 'car')
        got = stage1.convert('GEO', 'car')
        np.testing.assert_allclose(got.data, expected.data)
        np.testing.assert_equal(got.dtype, 'GEO')

    def test_GEO_GSE(self):
        """Regression test for IRBEM call"""
        test_cc = spc.Coords([1, 2, 3], 'GEO', 'car', ticks=Ticktock([2459213.5], 'JD'))
        expected = [[-2.118751061419987, -1.8297009536234092, 2.4825253744601286]]
        got = test_cc.convert('GSE', 'car')
        np.testing.assert_allclose(got.data, expected)
        self.assertEqual(got.dtype, 'GSE')

    def test_GDZ_in_kilometers(self):
        """Explicitly set units should be respected"""
        tt = Ticktock(2459218.5, 'JD')
        pos = np.array([2367.83158, -5981.75882, -4263.24591])
        cc_km = spc.Coords(pos, 'GEI', 'car', ticks=tt, units=['km', 'km', 'km'])
        got = cc_km.convert('GDZ', 'sph')
        cc_Re = spc.Coords(pos, 'GEI', 'car', ticks=tt, units=['Re', 'Re', 'Re'])
        expected = cc_Re.convert('GDZ', 'sph')
        np.testing.assert_approx_equal(expected.lati, got.lati, significant=6)
        np.testing.assert_approx_equal(expected.long, got.long, significant=6)
        np.testing.assert_approx_equal(expected.radi, got.radi, significant=4)
        self.assertEqual(got.units[0], 'km')
        self.assertEqual(expected.units[0], 'km')

    def test_GDZ_cartesian_raises_convert(self):
        """Geodetic coordinates shouldn't be expressed in Cartesian"""
        self.cvals.ticks = Ticktock([2459218.5]*2, 'JD')
        self.assertRaises(ValueError, self.cvals.convert, 'GDZ', 'car')

    def test_GDZ_cartesian_raises_creation(self):
        """Geodetic coordinates shouldn't be expressed in Cartesian"""
        tt = Ticktock(2459218.5, 'JD')
        pos = np.array([2367.83158, -5981.75882, -4263.24591])
        self.assertRaises(ValueError, spc.Coords, pos, 'GDZ', 'car', ticks=tt)

    def test_GEI_is_TOD(self):
        """IRBEM inertial isn't labelled, show it's TOD, i.e. GEO Z is TOD Z"""
        self.cvals.ticks = Ticktock([2459218.5]*2, 'JD')
        inertial = self.cvals.convert('GEI', 'car')
        np.testing.assert_allclose(inertial.z, self.cvals.z)

    def test_GSE_to_SMsph_vs_spacepy(self):
        """Compare outputs on GSE Cartesian to SM spherical conversion"""
        tt = Ticktock([2459218.5]*2, 'JD')
        irb = spc.Coords(self.cvals.data, 'GSE', 'car', ticks=tt, use_irbem=True)
        non_irb = spc.Coords(self.cvals.data, 'GSE', 'car', ticks=tt, use_irbem=False)
        irb_got = irb.convert('SM', 'sph')
        nonirb_got = non_irb.convert('SM', 'sph')
        # Test is approx. as IRBEM systems are relative to TOD not MOD
        np.testing.assert_allclose(irb_got.data, nonirb_got.data, rtol=5e-3)
        self.assertEqual(irb_got.dtype, nonirb_got.dtype)
        self.assertEqual(irb_got.carsph, nonirb_got.carsph)
        self.assertEqual(irb_got.units, nonirb_got.units)

    def test_SM_to_MAG_vs_spacepy(self):
        """Compare outputs on GSE Cartesian to SM spherical conversion"""
        tt = Ticktock([2459218.5]*2, 'JD')
        irb = spc.Coords(self.cvals.data, 'SM', 'car', ticks=tt, use_irbem=True)
        non_irb = spc.Coords(self.cvals.data, 'SM', 'car', ticks=tt, use_irbem=False)
        irb_got = irb.convert('CDMAG', 'car')
        nonirb_got = non_irb.convert('CDMAG', 'car')
        # Test is approx. as IRBEM systems are relative to TOD not MOD
        np.testing.assert_allclose(irb_got.data, nonirb_got.data, rtol=1e-3)
        self.assertEqual(irb_got.dtype, nonirb_got.dtype)
        self.assertEqual(irb_got.carsph, nonirb_got.carsph)
        self.assertEqual(irb_got.units, nonirb_got.units)


class QuaternionFunctionTests(unittest.TestCase):
    """Test of quaternion-related functions"""

    def test_quaternionNormalize(self):
        """quaternionNormalize should have known results"""
        tst = spc.quaternionNormalize([0.707, 0, 0.707, 0.2])
        ans = [ 0.69337122,  0.        ,  0.69337122,  0.19614462]
        np.testing.assert_array_almost_equal(ans, tst)

    def test_quaternionNormalize_2(self):
        """quaternionNormalize should have known result and magnitude 1"""
        tst1 = spc.quaternionNormalize([1, 0, 1, 0], scalarPos='first')
        tst2 = spc.quaternionNormalize([1, 0, 1, 0], scalarPos='last')
        ans = [0.70710678, 0.0, 0.70710678, 0]
        np.testing.assert_array_almost_equal(ans, tst1)
        np.testing.assert_array_almost_equal(ans, tst2)
        np.testing.assert_almost_equal(1.0, np.linalg.norm(tst1))
        np.testing.assert_almost_equal(1.0, np.linalg.norm(tst2))

    def test_quaternionNormalize_small(self):
        """test quaternionNormalize for very small values"""
        tst1 = spc.quaternionNormalize([1e-15, 0, 1e-15, 0], scalarPos='first')
        tst2 = spc.quaternionNormalize([1e-15, 0, 1e-15, 0], scalarPos='last')
        ans1 = [1.0, 0.0, 0.0, 0.0]
        ans2 = [0.0, 0.0, 0.0, 1.0]
        np.testing.assert_array_almost_equal(ans1, tst1)
        np.testing.assert_array_almost_equal(ans2, tst2)

    def test_quaternionNormalize_leadingdim(self):
        """test quaternionNormalize with leading degenerate dimensions"""
        tst = spc.quaternionNormalize([[1, 0, 1, 0]])
        ans = np.array([[0.7071068, 0, 0.7071068, 0]])
        np.testing.assert_array_almost_equal(ans, tst)

    def test_quaternionMultiply(self):
        """quaternionMultiply should have known results"""
        q1 = [1.0, 0.0, 0.0, 0.0]
        q2 = [0.0, 0.0, 0.0, 1.0]
        ans = [0.0, 0.0, 0.0, 1.0]
        tst = spc.quaternionMultiply(q2, q1, scalarPos='first')
        np.testing.assert_array_almost_equal(ans, tst)

    def test_quaternionConjugate_last(self):
        tst = spc.quaternionConjugate([0.707, 0, 0.707, 0.2], scalarPos='last')
        ans = [ -0.707,  -0.        ,  -0.707,  0.2]
        np.testing.assert_array_almost_equal(ans, tst)

    def test_quaternionConjugate_first(self):
        tst = spc.quaternionConjugate([0.2, 0.707, 0, 0.707], scalarPos='first')
        ans = [ 0.2,  -0.707,  -0.        ,  -0.707]
        np.testing.assert_array_almost_equal(ans, tst)

    def test_quaternionRotateVector(self):
        """Simple vector rotations"""
        cos45 = 0.5 ** 0.5 # 1/sqrt(2), or cos/sin of 45 degrees
        # No rotation, 90 degrees around each of X, Y, and Z
        Qin = np.array([
            [0, 0, 0, 1],
            [cos45, 0, 0, cos45],
            [0, cos45, 0, cos45],
            [0, 0, cos45, cos45],
            ])
        invect = np.array([
            [1, 0, 0],
            [0, 0, 1],
            [1, 0, 0],
            [0, 1, 0]
        ])
        expected = np.array([
            [1, 0, 0],
            [0, -1, 0],
            [0, 0, -1],
            [-1, 0, 0]
        ])
        outvect = spc.quaternionRotateVector(Qin, invect)
        np.testing.assert_array_almost_equal(outvect, expected)

    def test_quaternionFromMatrix_simple(self):
        """Test several simple rotations"""
        # Identity, and rotations by 90 degrees around X, Y, and Z axis
        inputs = np.array([
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            [[1, 0, 0], [0, 0, -1], [0, 1, 0]],
            [[0, 0, 1], [0, 1, 0], [-1, 0, 0]],
            [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
            ])
        cos45 = 0.5 ** 0.5 # 1/sqrt(2), or cos/sin of 45 degrees
        expected = np.array([
            [0, 0, 0, 1],
            [cos45, 0, 0, cos45],
            [0, cos45, 0, cos45],
            [0, 0, cos45, cos45],
            ])
        # Test single rotation at a time
        for i in range(expected.shape[0]):
            np.testing.assert_array_almost_equal(
                spc.quaternionFromMatrix(inputs[i, ...]),
                expected[i, ...])
        # Whole array at once
        actual = spc.quaternionFromMatrix(inputs)
        np.testing.assert_array_almost_equal(actual, expected)
        # Put scalar on other side
        expected = np.array([
            [1, 0, 0, 0],
            [cos45, cos45, 0, 0],
            [cos45, 0, cos45, 0],
            [cos45, 0, 0, cos45],
            ])
        actual = spc.quaternionFromMatrix(inputs, scalarPos='first')
        np.testing.assert_array_almost_equal(actual, expected)

    def test_quaternionFromMatrix_4D(self):
        """Simple rotations with 4D input"""
        # This is identical to simple tests, but with a different shape
        # Identity, and rotations by 90 degrees around X, Y, and Z axis
        inputs = np.array([
            [[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
             [[1, 0, 0], [0, 0, -1], [0, 1, 0]]],
            [[[0, 0, 1], [0, 1, 0], [-1, 0, 0]],
             [[0, -1, 0], [1, 0, 0], [0, 0, 1]]],
            ])
        cos45 = 0.5 ** 0.5 # 1/sqrt(2), or cos/sin of 45 degrees
        expected = np.array([
            [[0, 0, 0, 1],
             [cos45, 0, 0, cos45]],
            [[0, cos45, 0, cos45],
             [0, 0, cos45, cos45]],
            ])
        actual = spc.quaternionFromMatrix(inputs)
        np.testing.assert_array_almost_equal(actual, expected)
        # Put scalar on other side
        expected = np.array([
            [[1, 0, 0, 0],
             [cos45, cos45, 0, 0]],
            [[cos45, 0, cos45, 0],
             [cos45, 0, 0, cos45]],
            ])
        actual = spc.quaternionFromMatrix(inputs, scalarPos='first')
        np.testing.assert_array_almost_equal(actual, expected)

    def test_quaternionFromMatrix_perturbed(self):
        """Add error to a rotation matrix"""
        # Rotation by 90 degrees around X axis
        matrix = np.array(
            [[1., 0, 0], [0, 0, -1], [0, 1, 0]],
        )
        cos45 = 0.5 ** 0.5 # 1/sqrt(2), or cos/sin of 45 degrees
        # Equivalent quaternion
        expected = np.array([cos45, 0, 0, cos45],)
        # Add error, make sure still comes up with something reasonable
        np.random.seed(0x0d15ea5e)
        err = np.random.rand(3, 3) / 50 - 0.01 #-0.01 to 0.01
        matrix += err
        actual = spc.quaternionFromMatrix(matrix)
        np.testing.assert_array_almost_equal(actual, expected, decimal=3)

    def test_quaternionFromMatrix_nasty(self):
        """Pick an arbitrary rotation and verify it works"""
        # https://csm.mech.utah.edu/content/wp-content/uploads/2011/08/orthList.pdf
        # Axis of rotation
        u = [12. / 41, -24. / 41, 31. / 41]
        # Rotation angle
        theta = np.radians(58)
        ux, uy, uz = u
        c = np.cos(theta)
        s = np.sin(theta)
        # Construct rotation matrix from axis and angle
        # This might be doable more nicely in matrix notation...
        matrix = np.array([
            [c + ux ** 2 * (1 - c),
             ux * uy * (1 - c) - uz * s,
             ux * uz * (1 - c) + uy * s],
            [uy * ux * (1 - c) + uz * s,
             c + uy ** 2 * (1 - c),
             uy * uz * (1 - c) - ux * s],
            [uz * ux * (1 - c) - uy * s,
             uz * uy * (1 - c) + ux * s,
             c + uz ** 2 * (1 - c)]
        ])
        Qout = spc.quaternionFromMatrix(matrix)
        # Sample inputs to rotate
        invect = np.array([[5, 3, 2], [1, 0, 0], [.2, 5, 20],
                              [0, 2, 2]])
        # Transform the row vectors into column vectors so the
        # numpy multiplication gives the right result (then
        # transform back to row vectors for comparison.)
        expected = np.dot(matrix, invect.transpose()).transpose()
        actual = spc.quaternionRotateVector(
            np.tile(Qout, (4, 1)), invect)
        np.testing.assert_array_almost_equal(
            actual, expected)

    def test_quaternionFromMatrix_rt(self):
        """Round-trip arbitrary rotation matrix to quaternion and back"""
        # Same matrix as test_quaternionFromMatrix_nasty
        u = [12. / 41, -24. / 41, 31. / 41]
        theta = np.radians(58)
        ux, uy, uz = u
        c = np.cos(theta)
        s = np.sin(theta)
        matrix = np.array([
            [c + ux ** 2 * (1 - c),
             ux * uy * (1 - c) - uz * s,
             ux * uz * (1 - c) + uy * s],
            [uy * ux * (1 - c) + uz * s,
             c + uy ** 2 * (1 - c),
             uy * uz * (1 - c) - ux * s],
            [uz * ux * (1 - c) - uy * s,
             uz * uy * (1 - c) + ux * s,
             c + uz ** 2 * (1 - c)]
        ])
        Qout = spc.quaternionFromMatrix(matrix)
        matrix_rt = spc.quaternionToMatrix(Qout)
        np.testing.assert_array_almost_equal(
            matrix_rt, matrix)

    def test_quaternionFromMatrix_errors(self):
        """Test bad input"""
        matrix = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]])
        with self.assertRaises(NotImplementedError) as cm:
            spc.quaternionFromMatrix(matrix, 'FOO')
        self.assertEqual(
            'quaternionFromMatrix: scalarPos must be set to "First" or "Last"',
            str(cm.exception))
        with self.assertRaises(ValueError) as cm:
            spc.quaternionFromMatrix([[1, 2, 3]])
        self.assertEqual(
            'Input does not appear to be 3D rotation matrix, wrong size.',
            str(cm.exception))
        with self.assertRaises(ValueError) as cm:
            spc.quaternionFromMatrix(
                [[1, 1, 1], [2, 2, 2], [3, 3, 3]])
        self.assertEqual(
            'Input rotation matrix not orthogonal.',
            str(cm.exception))
        with self.assertRaises(ValueError) as cm:
            spc.quaternionFromMatrix([
                [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                [[1, 1, 1], [2, 2, 2], [3, 3, 3]]
            ])
        self.assertEqual(
            'Input rotation matrix at (1,) not orthogonal.',
            str(cm.exception))
        with self.assertRaises(ValueError) as cm:
            spc.quaternionFromMatrix(
                [[1, 0, 0], [0, -1, 0], [0, 0, 1]])
        self.assertEqual(
            'Input rotation matrix at () not proper.',
            str(cm.exception))

    def test_quaternionToMatrix_simple(self):
        """Test several simple rotations"""
        # Rotations by 90 degrees around X, Y, and Z axis
        cos45 = 0.5 ** 0.5 # 1/sqrt(2), or cos/sin of 45 degrees
        inputs = np.array([
            [cos45, 0, 0, cos45],
            [0, cos45, 0, cos45],
            [0, 0, cos45, cos45],
            ])
        expected = np.array([
            [[1, 0, 0], [0, 0, -1], [0, 1, 0]],
            [[0, 0, 1], [0, 1, 0], [-1, 0, 0]],
            [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
            ])
        actual = spc.quaternionToMatrix(inputs)
        np.testing.assert_array_almost_equal(actual, expected)
        # Put scalar on other side
        inputs = np.array([
            [cos45, cos45, 0, 0],
            [cos45, 0, cos45, 0],
            [cos45, 0, 0, cos45],
            ])
        actual = spc.quaternionToMatrix(inputs, scalarPos='first')
        np.testing.assert_array_almost_equal(actual, expected)

    def test_quaternionToMatrix_nasty(self):
        """Pick an arbitrary rotation and verify it works"""
        # Numbers pulled out of air
        Qin = spc.quaternionNormalize(np.array([0.25, 0.5, 0.71, 0.25]))
        matrix = spc.quaternionToMatrix(Qin)
        # Verify it's a rotation matrix
        ortho_test = np.dot(matrix, matrix.transpose())
        np.testing.assert_array_almost_equal(
            ortho_test, np.identity(3))
        det = np.linalg.det(matrix)
        self.assertTrue(det > 0) # Proper?
        invect = np.array([[5, 3, 2], [1, 0, 0],
                           [0.2, 5, 20], [0, 2, 2]])
        # Test matrix vs. quaternion rotation, single vector
        expected = spc.quaternionRotateVector(Qin, invect[1, :])
        actual = np.dot(matrix, invect[1, :])
        np.testing.assert_array_almost_equal(
            actual, expected)
        # All vectors at once
        expected = spc.quaternionRotateVector(
            np.tile(Qin, (4, 1)), invect)
        # Transform the row vectors into column vectors so the
        # numpy multiplication gives the right result (then
        # transform back to row vectors for comparison.)
        actual = np.dot(matrix, invect.transpose()).transpose()
        np.testing.assert_array_almost_equal(
            actual, expected)

    def testQuaternionToMatrixRT(self):
        """Round-trip test quaternion to matrix and back"""
        # Numbers pulled out of air
        Qin = spc.quaternionNormalize(np.array([0.25, 0.5, 0.71, 0.25]))
        matrix = spc.quaternionToMatrix(Qin)
        Qrt = spc.quaternionFromMatrix(matrix)
        if np.sign(Qrt[-1]) != np.sign(Qin[-1]):
            Qrt *= -1 #Handle the sign ambiguity
        np.testing.assert_array_almost_equal(
            Qrt, Qin)

    def test_quaternionToMatrix_errors(self):
        """Test bad input"""
        # Rotation by 90 degrees around X axis
        Qin = np.array([0.5 ** 0.5, 0, 0, 0.5 ** 0.5])
        with self.assertRaises(NotImplementedError) as cm:
            spc.quaternionToMatrix(Qin, 'FOO')
        self.assertEqual(
            'quaternionToMatrix: scalarPos must be set to "First" or "Last"',
            str(cm.exception))
        with self.assertRaises(ValueError) as cm:
            spc.quaternionToMatrix([1, 2, 3])
        self.assertEqual(
            'Input does not appear to be quaternion, wrong size.',
            str(cm.exception))
        with self.assertRaises(ValueError) as cm:
            spc.quaternionToMatrix([1, 2, 3, 4], normalize=False)
        self.assertEqual(
            'Input quaternion not normalized.',
            str(cm.exception))
        actual = spc.quaternionToMatrix([1, 2, 3, 4])
        expected = spc.quaternionToMatrix(spc.quaternionNormalize([1, 2, 3, 4]))
        np.testing.assert_array_almost_equal(
            actual, expected)


if __name__ == "__main__":
    ## suite = unittest.TestLoader().loadTestsFromTestCase(coordsTest)
    ## unittest.TextTestRunner(verbosity=2).run(suite)

    unittest.main()
