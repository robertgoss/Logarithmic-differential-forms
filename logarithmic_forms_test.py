#!/usr/local/bin/sage

import unittest

from logarithmic_forms import homogenous_wieghts
from logarithmic_forms import NotWieghtHomogeneousException
from logarithmic_forms import convert_polynomial_to_symbolic
from logarithmic_forms import convert_symbolic_to_polynomial
from logarithmic_forms import skew_iter
from logarithmic_forms import LogarithmicDifferentialForms

from singular_module import SingularModule

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
    
  def test_0_modules_crossing_ngens(self):
    crossing = self.x*self.y*self.z
    logdf = LogarithmicDifferentialForms(crossing)
    crossing_0_module = SingularModule([[crossing]])
    self.assertTrue(crossing_0_module.equals(logdf.p_module(0)))
    
  def test_1_modules_crossing_ngens(self):
    x = self.x
    y = self.y
    z = self.z
    zero = self.poly_ring.zero()
    crossing = x*y*z
    logdf = LogarithmicDifferentialForms(crossing)
    crossing_1_module = SingularModule([[y*z,zero,zero],[zero,x*z,zero],[zero,zero,x*y]])
    self.assertTrue(crossing_1_module.equals(logdf.p_module(1)))
    
  def test_2_modules_crossing_ngens(self):
    x = self.x
    y = self.y
    z = self.z
    zero = self.poly_ring.zero()
    crossing = self.x*self.y*self.z
    logdf = LogarithmicDifferentialForms(crossing)
    crossing_2_module = SingularModule([[z,zero,zero],[zero,y,zero],[zero,zero,x]])
    self.assertTrue(crossing_2_module.equals(logdf.p_module(2)))
    
  def test_3_modules_crossing_ngens(self):
    crossing = self.x*self.y*self.z
    logdf = LogarithmicDifferentialForms(crossing)
    crossing_3_module = SingularModule([[self.poly_ring.one()]])
    self.assertTrue(crossing_3_module.equals(logdf.p_module(3)))

  def test_zero_part_crossing_forms(self):
    crossing = self.x*self.y*self.z
    logdf = LogarithmicDifferentialForms(crossing)
    for i in range(4):
      self.assertEqual(len(logdf.p_forms_zero_basis(i)),len(logdf.p_form_generators(i)))
    
  def test_complement_homology(self):
    crossing = self.x*self.y*self.z
    logdf = LogarithmicDifferentialForms(crossing)
    hom = logdf.complement_homology()
    betti = {}
    for i,h in hom.iteritems():
      betti[i] = len(h)
    self.assertEqual(betti,{0:1,1:3,2:3,3:1})

  def test_equivarient_homology_whitney(self):
    crossing = self.x**2*self.y - self.z**2
    logdf = LogarithmicDifferentialForms(crossing)
    hom = logdf.equivarient_homology()
    betti = {}
    for i,h in hom.iteritems():
      betti[i] = len(h)
    self.assertEqual(betti,{0:1,1:0,2:0,3:0})

if __name__=="__main__":
  unittest.main()
