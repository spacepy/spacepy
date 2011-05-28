# -*- coding: utf-8 -*-

import unittest
import numpy as np
import spacepy.time as spt
import spacepy.empiricals as em


class empFunctionTests(unittest.TestCase):
    def setUp(self):
        super(empFunctionTests, self).setUp()
        self.ticks = spt.tickrange('2002-01-01T12:00:00','2002-01-04T00:00:00',.25)
        
    def tearDown(self):
        super(empFunctionTests, self).tearDown()
        
    def test_getPlasmaPause_regress(self):
        """regression test for plasmapause location"""
        real_ans = np.array([ 6.3475,  6.3475,  6.3475,  6.3475,  6.3475,  6.3475,  6.3475,
                           6.1775,  5.625 ,  5.4975,  5.4975])
        ans = em.getPlasmaPause(self.ticks, LT=12)
        for i, val in enumerate(real_ans):
            self.assertAlmostEqual(val, ans[i])
        
        real_ans = np.array([ 3.76 ,  3.76 ,  4.358,  4.358,  4.358,  4.358,  4.358,  4.358,
                            4.358,  4.542,  5.14])
        ans = em.getPlasmaPause(self.ticks, 'CA1992')
        for i, val in enumerate(real_ans):
            self.assertAlmostEqual(val, ans[i])
            
        real_ans = np.array([ 6.4214,  6.4214,  6.4214,  6.4214,  6.4214,  6.4214,  6.4214,
                            6.2686,  5.772 ,  5.6574,  5.6574])
        ans = em.getPlasmaPause(self.ticks)
        for i, val in enumerate(real_ans):
            self.assertAlmostEqual(val, ans[i])
            
    def test_getPlasmaPause(self):
        '''tests for exceptions in getPlasmaPause'''
        #check for fail on bad LT
        foo = lambda: em.getPlasmaPause(self.ticks, LT='bad')
        self.assertRaises(ValueError, foo)
        #check for fail on bad model
        bar = lambda: em.getPlasmaPause(self.ticks, model=3)
        self.assertRaises(ValueError, bar)
        #check for fail on LT out of legal range
        spam = lambda: em.getPlasmaPause(self.ticks, LT=25)
        self.assertRaises(IndexError, spam)
    

if __name__ == "__main__":
    ## suite = unittest.TestLoader().loadTestsFromTestCase(SimpleFunctionTests)
    ## unittest.TextTestRunner(verbosity=2).run(suite)

    ## suite = unittest.TestLoader().loadTestsFromTestCase(tFunctionTests)
    ## unittest.TextTestRunner(verbosity=2).run(suite)


    unittest.main()





