#!/usr/local/bin/sage

import unittest

from graded_module import wieghted_max_degree
from graded_module import wieghted_min_degree
from graded_module import monomials_of_order
from graded_module import GradedModule

from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.sets.set import Set

class TestWeightedMaxDegree(unittest.TestCase):

  def setUp(self):
    self.poly_ring = PolynomialRing(QQ,"x",3);
    self.x = self.poly_ring.gens()[0];
    self.y = self.poly_ring.gens()[1];
    self.z = self.poly_ring.gens()[2];

  def test_zero(self):
    self.assertEqual(wieghted_max_degree(self.poly_ring.zero(),[1,1,1]),0);
  
  def test_homogeneous(self):
    x = self.x;
    y = self.y;
    z = self.z;
    f = x**4 + 4*y*z**3 - z**2;
    self.assertEqual(wieghted_max_degree(f,[1,1,1]),4);

  def test_non_homogeneous(self):
    x = self.x;
    y = self.y;
    z = self.z;
    f = x**4 + 4*y*z**3 - z**2;
    self.assertEqual(wieghted_max_degree(f,[2,1,2]),8);

  def test_negative_wieghts(self):
    x = self.x;
    y = self.y;
    z = self.z;
    f = x**4 + 4*y*z**3 - z**2;
    self.assertEqual(wieghted_max_degree(f,[-1,-1,1]),2);

class TestWeightedMinDegree(unittest.TestCase):

  def setUp(self):
    self.poly_ring = PolynomialRing(QQ,"x",3);
    self.x = self.poly_ring.gens()[0];
    self.y = self.poly_ring.gens()[1];
    self.z = self.poly_ring.gens()[2];

  def test_zero(self):
    self.assertEqual(wieghted_min_degree(self.poly_ring.zero(),[1,1,1]),0);
  
  def test_homogeneous(self):
    x = self.x;
    y = self.y;
    z = self.z;
    f = x**4 + 4*y*z**3 - z**2;
    self.assertEqual(wieghted_min_degree(f,[1,1,1]),2);

  def test_non_homogeneous(self):
    x = self.x;
    y = self.y;
    z = self.z;
    f = x**4 + 4*y*z**3 - z**2;
    self.assertEqual(wieghted_min_degree(f,[2,1,2]),4);

  def test_negative_wieghts(self):
    x = self.x;
    y = self.y;
    z = self.z;
    f = x**4 + 4*y*z**3 - z**2;
    self.assertEqual(wieghted_min_degree(f,[-1,-1,1]),-4);

class TestMonomialsOfOrder(unittest.TestCase):
 
  def setUp(self):
    self.poly_ring = PolynomialRing(QQ,"x",3);
    self.x = self.poly_ring.gens()[0];
    self.y = self.poly_ring.gens()[1];
    self.z = self.poly_ring.gens()[2];
  
  def test_zero(self):
    mons = [mon for mon in monomials_of_order(0,self.poly_ring,[1,1,1])]
    self.assertEqual(mons,[self.poly_ring.one()])

  def test_homogeneous_3(self):
    x = self.x;
    y = self.y;
    z = self.z;
    true_mons = Set([x**3,y**3,z**3,x**2*y,x**2*z,y**2*x,y**2*z,z**2*x,z**2*y,x*y*z])
    mons = [mon for mon in monomials_of_order(3,self.poly_ring,[1,1,1])]
    self.assertEqual(true_mons,Set(mons))

  def test_non_homogeneous_4(self):
    x = self.x;
    y = self.y;
    z = self.z;
    true_mons = Set([x**4,x**2*y,x*z,y**2])
    mons = [mon for mon in monomials_of_order(4,self.poly_ring,[1,2,3])]
    self.assertEqual(true_mons,Set(mons))
   
class TestGradedModule(unittest.TestCase):
 
  def setUp(self):
    self.poly_ring = PolynomialRing(QQ,"x",3);
    self.x = self.poly_ring.gens()[0];
    self.y = self.poly_ring.gens()[1];
    self.z = self.poly_ring.gens()[2];

  def test_monomial_basis_zero(self):
    one = self.poly_ring.one()
    zero = self.poly_ring.zero()
    gm = GradedModule([[one,one,one,one]],[0,1,2,3],[1,2,3])
    self.assertEqual(gm.monomial_basis(0),[(one,zero,zero,zero)])

  def test_monomial_basis(self):
    x = self.x
    y = self.y
    z = self.z
    one = self.poly_ring.one()
    zero = self.poly_ring.zero()
    gm = GradedModule([[one,one]],[0,1],[1,2,3])
    true_basis  = [(x**2,zero),(y,zero),(zero,x)]
    self.assertEqual(Set(gm.monomial_basis(2)),Set(true_basis))

  def test_homogeneous_parts_A(self):
    x = self.x
    y = self.y
    z = self.z
    one = self.poly_ring.one()
    zero = self.poly_ring.zero()
    gm = GradedModule([[one,one]],[0,1],[1,2,3])
    parts = gm.get_homogeneous_parts([one,one])
    parts_true = { 0:[one,zero] , 1:[zero,one] }
    self.assertEqual(parts,parts_true)

  def test_homogeneous_parts_B(self):
    x = self.x
    y = self.y
    z = self.z
    one = self.poly_ring.one()
    zero = self.poly_ring.zero()
    gm = GradedModule([[one,one]],[0,1],[1,2,3])
    parts = gm.get_homogeneous_parts([x*y,x**3+z*y])
    self.assertEqual(parts,{3:[x*y,zero],4:[zero,x**3],6:[zero,z*y]})
  
  def test_homogeneous_part_basisA(self):
    x = self.x
    y = self.y
    z = self.z
    one = self.poly_ring.one()
    zero = self.poly_ring.zero()
    gm = GradedModule([[z,one,x**2 + y]],[0,1,2],[1,2,3])
    basis = gm.homogeneous_part_basis(6);
    self.assertEqual(len(basis),10)
 
  def test_homogeneous_part_basisB(self):
    #From bug found with ncd
    x = self.x
    y = self.y
    z = self.z
    one = self.poly_ring.one()
    zero = self.poly_ring.zero()
    gm = GradedModule([[zero, x*z, x*y], [zero, -x*z, zero], [y*z, zero, zero]],[1, 1, 1],[1, 1, 1])
    basis = gm.homogeneous_part_basis(3);
    self.assertEqual(len(basis),3)

if __name__=="__main__":
  unittest.main()
