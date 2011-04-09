#!/usr/local/bin/sage

import unittest

from logarithmic_form import convert_polynomial_to_symbolic
from logarithmic_form import convert_symbolic_to_polynomial


from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.symbolic.ring import var

from sage.tensor.differential_form_element import DifferentialForm

class TestConvertSymToPoly(unittest.TestCase):

  def setUp(self):
    self.poly_ring = PolynomialRing(QQ,"x",3);
    self.x = self.poly_ring.gens()[0];
    self.y = self.poly_ring.gens()[1];
    self.z = self.poly_ring.gens()[2];
    self.vars = var('x,y,z')
    
  def test_zero(self):
    zero = 0*self.vars[0]
    poly = convert_symbolic_to_polynomial(zero,self.poly_ring,self.vars)
    self.assertEqual(poly,self.poly_ring.zero())
    
  def test_convert(self):
    x = self.x
    y = self.y
    z = self.z
    sym_poly = 4*self.vars[0]**4 + self.vars[1]*self.vars[0]**12*self.vars[1]-self.vars[2]
    poly = convert_symbolic_to_polynomial(sym_poly,self.poly_ring,self.vars)
    self.assertEqual(poly,4*x**4+y*x**12*y-z)

  def test_convert_univarient(self):
    x = self.x
    y = self.y
    z = self.z
    sym_poly = 4*self.vars[1]**4 + self.vars[1]*self.vars[1]**12*self.vars[1]-self.vars[1]
    poly = convert_symbolic_to_polynomial(sym_poly,self.poly_ring,self.vars)
    self.assertEqual(poly,4*y**4+y*y**12*y-y)

  def test_convert_partial(self):
    x = self.x
    y = self.y
    z = self.z
    sym_poly = 4*self.vars[2]**4 + self.vars[1]*self.vars[1]**12*self.vars[1]-self.vars[1]
    poly = convert_symbolic_to_polynomial(sym_poly,self.poly_ring,self.vars)
    self.assertEqual(poly,4*z**4+y*y**12*y-y)
    
class TestConvertPolyToSym(unittest.TestCase):

  def setUp(self):
    self.poly_ring = PolynomialRing(QQ,"x",3);
    self.x = self.poly_ring.gens()[0];
    self.y = self.poly_ring.gens()[1];
    self.z = self.poly_ring.gens()[2];
    self.vars = var('x,y,z')
    
  def test_zero(self):
    zero = self.poly_ring.zero()
    poly = convert_polynomial_to_symbolic(zero,self.vars)
    self.assertEqual(poly,0)
    
  def test_convert(self):
    x = self.x
    y = self.y
    z = self.z
    poly = 4*x**4+y*x**12*y-z
    sym_poly = 4*self.vars[0]**4 + self.vars[1]*self.vars[0]**12*self.vars[1]-self.vars[2]
    con_poly = convert_polynomial_to_symbolic(poly,self.vars)
    self.assertEqual(sym_poly,con_poly)

if __name__=="__main__":
    unittest.main()
