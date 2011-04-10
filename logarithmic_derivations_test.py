#!/usr/local/bin/sage

import unittest

import sage.all

from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

from logarithmic_derivations import LogarithmicDerivations

from singular_module import SingularModule

class TestLogaritmicDerivations(unittest.TestCase):

  def setUp(self):
    self.poly_ring = PolynomialRing(QQ,"x",3);
    self.x = self.poly_ring.gens()[0];
    self.y = self.poly_ring.gens()[1];
    self.z = self.poly_ring.gens()[2];
    

  def test_normal_crossing(self):
    zero = self.poly_ring.zero()
    log_der = LogarithmicDerivations(self.x*self.y*self.z)
    free =  SingularModule([[self.x,zero,zero],[zero,self.y,zero],[zero,zero,self.z]])
    self.assertTrue(log_der.equals(free))
  
  
if __name__=="__main__":
  unittest.main()
