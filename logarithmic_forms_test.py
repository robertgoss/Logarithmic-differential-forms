#!/home/robert/sage-4.7.alpha2/sage-python

import unittest

from logarithmic_forms import homogenous_wieghts
from logarithmic_forms import NotWieghtHomogeneousException
from logarithmic_forms import convert_polynomial_to_symbolic
from logarithmic_forms import convert_symbolic_to_polynomial
from logarithmic_forms import skew_iter
from logarithmic_forms import LogarithmicDifferentialForms


from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.symbolic.ring import var

from sage.tensor.differential_form_element import DifferentialForm

    
class TestHomogeneousWieghts(unittest.TestCase):

  def setUp(self):
    self.poly_ring = PolynomialRing(QQ,"x",3);
    self.x = self.poly_ring.gens()[0];
    self.y = self.poly_ring.gens()[1];
    self.z = self.poly_ring.gens()[2];
  
  def test_homogenous(self):
    x = self.x
    y = self.y
    z = self.z
    divisor = x*y*z + x**3 + z**3
    wieghts = homogenous_wieghts(divisor)
    #Test this works with a homogenous divisor
    self.assertEquals(wieghts,[3,1,1,1]);
    
  def test_weighted_homogenous(self):
    x = self.x
    y = self.y
    z = self.z
    divisor = x**2*y-z**2
    wieghts = homogenous_wieghts(divisor)
    #Test this works with a weighted homogenous divisor
    self.assertEquals(wieghts,[4,1,2,2])
  
  def test_not_homogenous(self):
    x = self.x
    y = self.y
    z = self.z
    divisor = x**2*y-x**3 + z +y**2;
    self.assertRaises(NotWieghtHomogeneousException,homogenous_wieghts,divisor)
    
class TestConvertSymToPoly(unittest.TestCase):

  def setUp(self):
    self.poly_ring = PolynomialRing(QQ,"x",3);
    self.x = self.poly_ring.gens()[0];
    self.y = self.poly_ring.gens()[1];
    self.z = self.poly_ring.gens()[2];
    self.vars = var('x,y,z')
    
  def test_zero(self):
    zero = 0*self.vars[0]
    poly = convert_symbolic_to_polynomial(zero,self.poly_ring)
    self.assertEqual(poly,self.poly_ring.zero())
    
  def test_convert(self):
    x = self.x
    y = self.y
    z = self.z
    sym_poly = 4*self.vars[0]**4 + self.vars[1]*self.vars[0]**12*self.vars[1]-self.vars[2]
    poly = convert_symbolic_to_polynomial(sym_poly,self.poly_ring)
    self.assertEqual(poly,4*x**4+y*x**12*y-z)
    
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
    
class TestSkewIter(unittest.TestCase):

  def test_skew_iter(self):
    vars = []
    for index in skew_iter(3,2):
      vars.append(index)
    self.assertEquals(vars,[[0,1],[0,2],[1,2]])
    
  def test_depth_zero(self):
    vars = []
    for index in skew_iter(10,1):
      vars.append(index)
    self.assertEquals(vars,[[i] for i in range(10)])
   
  def test_satisfied(self):
    for index in skew_iter(7,3):
      self.assertTrue(index[2] > index[1])
      self.assertTrue(index[1] > index[0])
      
class TestLogartihmicDifferentialForms(unittest.TestCase):

  def setUp(self):
    self.poly_ring = PolynomialRing(QQ,"x",3);
    self.x = self.poly_ring.gens()[0];
    self.y = self.poly_ring.gens()[1];
    self.z = self.poly_ring.gens()[2];
    self.vars = var('x,y,z')
    
  def test_p_forms_crossing_ngens(self):
    crossing = self.x*self.y*self.z
    logdf = LogarithmicDifferentialForms(crossing)
    self.assertEqual(len(logdf.p_form_generators(0)),1)
    self.assertEqual(len(logdf.p_form_generators(1)),3)
    self.assertEqual(len(logdf.p_form_generators(2)),3)
    self.assertEqual(len(logdf.p_form_generators(3)),1)

    
    
if __name__=="__main__":
  unittest.main()
