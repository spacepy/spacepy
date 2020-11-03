
import unittest
import spacepy.coordinates as spc
import glob
import os
import datetime
from numpy import array
import numpy as np
from spacepy.time import Ticktock
try:
    import spacepy.irbempy as ib
except ImportError:
    pass #tests will fail, but won't bring down the entire suite
import spacepy.toolbox as tb

__all__ = ['coordsTest']


class coordsTest(unittest.TestCase):
    def setUp(self):
        #super(tFunctionTests, self).setUp()
        try:
            self.cvals = spc.Coords([[1,2,4],[1,2,2]], 'GEO', 'car')
        except ImportError:
            pass #tests will fail, but won't bring down the entire suite

    def tearDown(self):
        #super(tFunctionTests, self).tearDown()
        pass

    def test_coords(self):
        """Coords should create and do simple conversions"""
        np.testing.assert_equal([1,1], self.cvals.x)
        np.testing.assert_equal([2,2], self.cvals.y)
        np.testing.assert_equal([4,2], self.cvals.z)
        self.cvals.ticks = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO') # add ticktock
        newcoord = self.cvals.convert('GSM', 'sph')

    def test_append(self):
        c2 = spc.Coords([[6,7,8],[9,10,11]], 'GEO', 'car')
        actual = self.cvals.append(c2)
        expected = [[1,2,4],[1,2,2],[6,7,8],[9,10,11]]
        np.testing.assert_equal(expected, actual.data.tolist())

    def test_slice(self):
        expected = spc.Coords([1,2,4], 'GEO', 'car')
        np.testing.assert_equal(expected.data, self.cvals[0].data)

    def test_slice_with_ticks(self):
        self.cvals.ticks = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO')
        expected = spc.Coords([1,2,4], 'GEO', 'car')
        np.testing.assert_equal(expected.data, self.cvals[0].data)


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
