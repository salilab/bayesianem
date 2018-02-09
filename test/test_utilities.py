import IMP
import IMP.core
import IMP.algebra
import IMP.test

import IMP.bayesianem


import os
import operator
import math

class GaussianEMRestraintRigidBody(IMP.test.TestCase):
    def setUp(self):
        IMP.test.TestCase.setUp(self)
        self.m = IMP.Model()

    def test_sanity(self):
        """test that functions are called"""
        a=IMP.bayesianem.get_masked_map


if __name__ == '__main__':
    IMP.test.main()
