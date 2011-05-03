#!/usr/local/bin/sage

import unittest

from logarithmic_forms import homogenous_wieghts
from logarithmic_forms import NotWieghtHomogeneousException
from logarithmic_forms import skew_iter
from logarithmic_forms import orth_complement
from logarithmic_forms import LogarithmicDifferentialForms

from singular_module import SingularModule

from sage.rings.rational_field import QQ
from sage.modules.free_module import VectorSpace
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.symbolic.ring import var


    
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

class TestOrthagomalComplement(unittest.TestCase):
  
  def test_zero(self):
    v_space = VectorSpace(QQ,4)
    zero = v_space.subspace([v_space.zero()])
    comp = orth_complement(v_space,zero)
    self.assertEqual(v_space,comp)
    
  def test_full(self):
    v_space = VectorSpace(QQ,4)
    zero = v_space.subspace([v_space.zero()])
    comp = orth_complement(v_space,v_space)
    self.assertEqual(zero,comp)

  def test_comp(self):
    v_space = VectorSpace(QQ,4)
    sub = v_space.subspace([[1,1,1,1],[2,2,-3,-1]])
    comp = orth_complement(v_space,sub)
    true_comp = v_space.subspace([[19,-17,3,-5],[1,-1,0,0]])
    self.assertEqual(comp,true_comp)

  def test_satisfy_inter(self):
    v_space = VectorSpace(QQ,4)
    sub = v_space.subspace([[1,-1,1,1],[2,-3,4,5]])
    comp = orth_complement(v_space,sub)
    zero = v_space.subspace([v_space.zero()])
    inter = sub.intersection(comp)
    self.assertEqual(zero,inter)

  def test_satisfy_sum(self):
    v_space = VectorSpace(QQ,4)
    sub = v_space.subspace([[1,1,1,-1],[2,-3,41,5]])
    comp = orth_complement(v_space,sub)
    self.assertEqual(v_space,sub+comp)
    
    
class TestSkewIter(unittest.TestCase):

  def test_skew_iter(self):
    vars = []
    for index in skew_iter(3,2):
      vars.append(index)
    self.assertEquals(vars,[[0,1],[0,2],[1,2]])
  
  def test_skew_iter_7_3(self):
    vars = []
    for index in skew_iter(7,3):
      vars.append(index)
    self.assertEquals(vars,[[0,1,2],[0,1,3],[0,1,4],[0,1,5],[0,1,6],
                            [0,2,3],[0,2,4],[0,2,5],[0,2,6],[0,3,4],
                            [0,3,5],[0,3,6],[0,4,5],[0,4,6],[0,5,6],
                            [1,2,3],[1,2,4],[1,2,5],[1,2,6],[1,3,4],
                            [1,3,5],[1,3,6],[1,4,5],[1,4,6],[1,5,6],
                            [2,3,4],[2,3,5],[2,3,6],[2,4,5],[2,4,6],
                            [2,5,6],[3,4,5],[3,4,6],[3,5,6],[4,5,6]])
    
  def test_skew_iter_6_3(self):
    vars = []
    for index in skew_iter(6,3):
      vars.append(index)
    self.assertEquals(vars,[[0,1,2],[0,1,3],[0,1,4],[0,1,5],
                            [0,2,3],[0,2,4],[0,2,5],[0,3,4],[0,3,5],[0,4,5],
                            [1,2,3],[1,2,4],[1,2,5],[1,3,4],[1,3,5],[1,4,5],
                            [2,3,4],[2,3,5],[2,4,5],[3,4,5]])
    
  def test_depth_zero(self):
    vars = []
    for index in skew_iter(10,1):
      vars.append(index)
    self.assertEquals(vars,[[i] for i in range(10)])
   
  def test_satisfied(self):
    for index in skew_iter(7,3):
      self.assertTrue(index[2] > index[1])
      self.assertTrue(index[1] > index[0])

  def test_edge(self):
    #Based on a bug issue found
    vars = []
    for index in skew_iter(5,4):
      vars.append(index)
    self.assertEquals(vars,[[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],[1,2,3,4]])
      
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

  def test_p_module_n_crossing(self):
    #Make sure this doesnt throw an error - fix bug
    for i in range(4,5):
      p_ring = PolynomialRing(QQ,i,"z")
      crossing = p_ring.one()
      for g in p_ring.gens():
        crossing *= g
      logdf = LogarithmicDifferentialForms(crossing)
      logdf.p_module(i-1)

  def test_complement_complex_crossing(self):
    crossing = self.x*self.y*self.z
    logdf = LogarithmicDifferentialForms(crossing)
    complex = logdf.chain_complex("complement")
    complex_size = {}
    for i,c in complex.iteritems():
      complex_size[i] = len(c)
    self.assertEqual(complex_size,{0:1,1:3,2:3,3:1})

  def test_complement_complex_whitney(self):
    whitney = self.x**2*self.y - self.z**2
    logdf = LogarithmicDifferentialForms(whitney)
    complex = logdf.chain_complex("complement")
    complex_size = {}
    for i,c in complex.iteritems():
      complex_size[i] = len(c)
    self.assertEqual(complex_size,{0:1,1:1,2:0,3:0})

  def test_equi_complex_crossing(self):
    crossing = self.x*self.y*self.z
    logdf = LogarithmicDifferentialForms(crossing)
    complex = logdf.chain_complex("equivarient")
    complex_size = {}
    for i,c in complex.iteritems():
      complex_size[i] = len(c)
    self.assertEqual(complex_size,{0:1,1:3,2:4,3:4})

  def test_equi_complex_whitney(self):
    whitney = self.x**2*self.y - self.z**2
    logdf = LogarithmicDifferentialForms(whitney)
    complex = logdf.chain_complex("equivarient")
    complex_size = {}
    for i,c in complex.iteritems():
      complex_size[i] = len(c)
    self.assertEqual(complex_size,{0:1,1:1,2:1,3:1})

  def test_complement_homology_crossing(self):
    crossing = self.x*self.y*self.z
    logdf = LogarithmicDifferentialForms(crossing)
    homology = logdf.homology("complement")
    homology_size = {}
    for i,c in homology.iteritems():
      homology_size[i] = len(c)
    self.assertEqual(homology_size,{0:1,1:3,2:3,3:1})

  def test_equi_homology_whitney(self):
    whitney = self.x**2*self.y-self.z**2
    logdf = LogarithmicDifferentialForms(whitney)
    homology = logdf.homology("equivarient")
    homology_size = {}
    for i,c in homology.iteritems():
      homology_size[i] = len(c)
    self.assertEqual(homology_size,{0:1,1:0,2:0,3:0})

  def test_relative_complex_0_crossing(self):
    crossing = self.x*self.y*self.z
    logdf = LogarithmicDifferentialForms(crossing)
    homology = logdf.chain_complex("relative",None,0)
    homology_size = {}
    for i,c in homology.iteritems():
      homology_size[i] = len(c)
    self.assertEqual(homology_size,{0:1,1:2,2:1,3:0})

  def test_relative_complex_0_whitney(self):
    whitney = self.x**2*self.y-self.z**2
    logdf = LogarithmicDifferentialForms(whitney)
    homology = logdf.chain_complex("relative",None,0)
    homology_size = {}
    for i,c in homology.iteritems():
      homology_size[i] = len(c)
    self.assertEqual(homology_size,{0:1,1:0,2:0,3:0})

  def test_relative_homology_0_crossing(self):
    crossing = self.x*self.y*self.z
    logdf = LogarithmicDifferentialForms(crossing)
    homology = logdf.homology("relative",0)
    homology_size = {}
    for i,c in homology.iteritems():
      homology_size[i] = len(c)
    self.assertEqual(homology_size,{0:1,1:2,2:1,3:0})

  def test_relative_homology_0_whitney(self):
    whitney = self.x**2*self.y-self.z**2
    logdf = LogarithmicDifferentialForms(whitney)
    homology = logdf.homology("relative",0)
    homology_size = {}
    for i,c in homology.iteritems():
      homology_size[i] = len(c)
    self.assertEqual(homology_size,{0:1,1:0,2:0,3:0})


if __name__=="__main__":
  unittest.main()
