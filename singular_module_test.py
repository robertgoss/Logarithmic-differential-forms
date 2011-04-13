#!/usr/local/bin/sage

import unittest

from singular_module import SingularModule

from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.symbolic.ring import var
from sage.rings.ideal import Ideal

class TestSingularModule(unittest.TestCase):

  def setUp(self):
    self.poly_ring = PolynomialRing(QQ,"x",3);
    self.x = self.poly_ring.gens()[0];
    self.y = self.poly_ring.gens()[1];
    self.z = self.poly_ring.gens()[2];
    
  def test_creation(self):
    x = self.x
    y = self.y
    z = self.z
    sm = SingularModule([[x,y+z,z**3-2*y],[x,y,z]])
    #Check this did cause an error
    self.assertTrue(True)
    
  def test_create_ring_str(self):
    x = self.x
    y = self.y
    z = self.z
    sm = SingularModule([[x,y+z,z**3-2*y],[x,y,z]])
    self.assertEqual(sm.create_ring_str(),"ring r=0,x(1..3),dp;\n");
    
  def test_create_module_str(self):
    x = self.x
    y = self.y
    z = self.z
    sm = SingularModule([[x,2*y**2,3*z**3]])
    self.assertEqual(sm.create_module_str("MT"),"module MT=[1*x(1)^1,2*x(2)^2,3*x(3)^3];\n groebner(MT);\n")
    
  def test_contains_zero(self):
    x = self.x
    y = self.y
    z = self.z
    sm = SingularModule([[x,y+z,z**3-2*y],[x,y,z]])
    zero = self.poly_ring.zero()
    #Assert we always contain zero
    self.assertTrue(sm.contains([zero,zero,zero]))
    
  def test_contains_gen(self):
    x = self.x
    y = self.y
    z = self.z
    sm = SingularModule([[x,y+z,z**3-2*y],[x,y,z]])
    #Assert we contain our generators
    self.assertTrue(sm.contains([x,y+z,z**3-2*y]))
    self.assertTrue(sm.contains([x,y,z]))
    
  def test_contains_combination(self):
    x = self.x
    y = self.y
    z = self.z
    sm = SingularModule([[x,y+z,z**3-2*y],[x,y,z]])
    #Assert we contain a combination of the generorators
    self.assertTrue(sm.contains([2*x,2*y+z,z**3-2*y+z]))
    
  def test_contains_fail(self):
    x = self.x
    y = self.y
    z = self.z
    one = self.poly_ring.one()
    sm = SingularModule([[x,y+z,z**3-2*y],[x,y,z]])
    #Detect that certain vectors are not contained
    self.assertFalse(sm.contains([one,x,z]))
    self.assertFalse(sm.contains([y,y,y]))
    self.assertFalse(sm.contains([z,y,x]))
    
  def test_contains_fail_multiple(self):
    x = self.x
    y = self.y
    z = self.z
    sm = SingularModule([[x**2,x*y,x*z]])
    #Detect that this is not contianed even though a multiple is
    self.assertFalse(sm.contains([x,y,z]))
    
  def test_contains_module(self):
    x = self.x
    y = self.y
    z = self.z
    zero = self.poly_ring.zero()
    smA = SingularModule([[x**2,x*y,y*z],[x,y,z]])
    smB = SingularModule([[x**2+x,x*y+y,y*z+z],[zero,zero,y*z-x*z]])
    self.assertTrue(smA.contains(smB))
    
  def test_contained_in_ambient(self):
    x = self.x
    y = self.y
    z = self.z
    smA = SingularModule([[x,y**3,z+x],[z,x,x**2]])
    am = smA.ambient_free_module()
    self.assertTrue(am.contains(smA))
    
  def test_intersect(self):
    x = self.x
    y = self.y
    z = self.z
    smA = SingularModule([[x,y**3,z+x],[z,x,x**2]])
    smB = SingularModule([[z,y,x],[-z,x**2,4*y + z]]);
    smI = smA.intersection(smB)
    gen_1_1 = x**5*z+x*y**3*z**2+4*y**4*z**2+y**3*z**3+x**3*y*z-x**3*z**2-x**2*z**3-x**3*z-4*x**2*y*z-x**2*z**2-x*y*z**2-y*z**3
    gen_2_1 = x**4*y**3*z-x**3*y**3*z+x**2*y**4*z+4*y**5*z+y**4*z**2+x**5-x**4*z-x**3*z**2-4*x**2*y**2-2*x**2*y*z-x*y*z**2
    gen_3_1 = x**3*y**3*z-x**3*y*z+4*x**2*y**4*z+x**2*y**3*z**2+x**6-4*x**3*y**2-x**4*z-x**3*z**2-x**3*z-4*x**2*y*z+4*x*y**2*z-2*x**2*z**2-3*x*y*z**2+4*y**2*z**2-x*z**3+y*z**3
    gens = [[gen_1_1,gen_2_1,gen_3_1]]
    #Test this example
    self.assertEqual(gens,smI.gens)

  def test_intersect_contains(self):
    x = self.x
    y = self.y
    z = self.z
    sm1 = SingularModule([[x,y+z,z**3-2*y],[x,y,z],[x,y,z**2]])
    sm2 = SingularModule([[x,y**2,z**3]])
    #Assert the intersection is contained in both modules
    gens = (sm1.intersection(sm2)).gens
    for gen in gens:
      self.assertTrue(sm1.contains(gen))
    for gen in gens:
      self.assertTrue(sm2.contains(gen))
    
  def test_intersection_symetry(self):
    x = self.x
    y = self.y
    z = self.z
    sm1 = SingularModule([[x,y+z,z**3-2*y],[x,y,z],[x,y,z**2]])
    sm2 = SingularModule([[x,y**2,z**3]])
    sm_1_2 = sm1.intersection(sm2)
    sm_2_1 = sm2.intersection(sm1)
    #Assert that equals is symetric
    self.assertTrue(sm_1_2.equals(sm_2_1))
    
  def test_reduce_lossless(self):
    x = self.x
    y = self.y
    z = self.z
    sm1 = SingularModule([[x,y+z,z**3-2*y],[x,y,z],[x,y,z**2]])
    sm2 = SingularModule([[x,y+z,z**3-2*y],[x,y,z],[x,y,z**2]])
    sm2.reduce_generators()
    self.assertTrue(sm1.equals(sm2))
    
  def test_standard_basis_irredundent(self):
    x = self.x
    y = self.y
    z = self.z
    sm1 = SingularModule([[x,y+z,z**3-2*y],[x,y,z],[x,y,z**2]])
    std_gens = sm1.standard_basis()
    std_mod = SingularModule(std_gens)
    self.assertTrue(sm1.equals(std_mod))
      
  def test_equals_A(self):
    x = self.x
    y = self.y
    z = self.z
    #Assert these are equal
    sm1 = SingularModule([[x,y+z,z**3-2*y],[x,y,z],[x,y,z**2]])
    sm2 = SingularModule([[x,y,z],[x,y+z,z**3-2*y],[x,y,z**2]])
    self.assertTrue(sm1.equals(sm2))
    
  def test_equals_B(self):
    x = self.x
    y = self.y
    z = self.z
    zero = self.poly_ring.zero()
    #Assert these are equal - from crossing divisor
    sm1 = SingularModule([[x,zero,zero],[zero,y,zero],[zero,zero,z]])
    sm2 = SingularModule([[zero,y,z],[zero,zero,z],[x,zero,zero]])
    self.assertTrue(sm1.equals(sm2))
    
  def test_equals_not(self):
    x = self.x
    y = self.y
    z = self.z
    sm1 = SingularModule([[x,y+z,z**3-2*y],[x,y,z],[x,y,z**2]])
    sm2 = SingularModule([[x,y,z]])
    #Assert these are not equal
    self.assertFalse(sm1.equals(sm2))
    
  def test_ambient_free_module(self):
    x = self.x
    y = self.y
    z = self.z
    one = self.poly_ring.one()
    zero = self.poly_ring.zero()
    sm1 = SingularModule([[x,y+z,z**3-2*y],[x,y,z],[x,y,z**2]])
    self.assertEqual(sm1.ambient_free_module().gens,[[one,zero,zero],[zero,one,zero],[zero,zero,one]])
    
  def test_is_free_trivial(self):
    free = SingularModule.create_free_module(3,self.poly_ring)
    self.assertTrue(free.is_free())
    
  def test_is_free(self):
    x = self.x
    y = self.y
    z = self.z
    zero = self.poly_ring.zero()
    free = SingularModule([[x,zero,zero],[zero,y,zero],[zero,zero,x**2]])
    self.assertTrue(free.is_free())
    
  def test_is_free_not(self):
    x = self.x
    y = self.y
    z = self.z
    sm1 = SingularModule([[x,x,x],[y,y,y]])
    self.assertFalse(sm1.is_free())
    
  def test_create_relationA(self):
    x = self.x
    y = self.y
    z = self.z
    one = self.poly_ring.one()
    zero = self.poly_ring.zero()
    relation = [x**2,one+z**2,y]
    ideal = Ideal(self.poly_ring,[x**2*y-z**2])
    mod = SingularModule.create_from_relation(relation,ideal)
    true_mod = SingularModule([[one,-x**2,x**4],[zero,y,-z**2-1],
                                [zero,z**2,-x**2*z**2-x**2],[zero,zero,x**2*y-z**2]])
    self.assertTrue(mod.equals(true_mod))

  def test_create_relationB(self):
    #From an error uncovered in log derivations
    x = self.x
    y = self.y
    z = self.z
    zero = self.poly_ring.zero()
    relation = [y*z,x*z,x*y]
    ideal = Ideal(self.poly_ring,[x*y*z])
    mod = SingularModule.create_from_relation(relation,ideal)
    true_mod = SingularModule([[x,zero,zero],[zero,y,zero],[zero,zero,z]])
    self.assertTrue(mod.equals(true_mod))
    
  def test_create_relation_satisfy_A(self):
    x = self.x
    y = self.y
    z = self.z
    zero = self.poly_ring.zero()
    relation = [x+x**4-y,(y+z)**3,-2*y+self.poly_ring.one()]
    ideal = Ideal(self.poly_ring,[x**2*y-z**2])
    mod = SingularModule.create_from_relation(relation,ideal)
    for gen in mod.gens:
      sum = zero
      for g,rel in zip(gen,relation):
        sum = sum +g*rel
    self.assertTrue(sum in ideal)
  
  def test_create_relation_satisfy_B(self):
    x = self.x
    y = self.y
    z = self.z
    one = self.poly_ring.one()
    zero = self.poly_ring.zero()
    relation = [x**2,one+z**2,y]
    ideal = Ideal(self.poly_ring,[x**2*y-z**2])
    mod = SingularModule.create_from_relation(relation,ideal)
    for gen in mod.gens:
      sum = zero
      for g,rel in zip(gen,relation):
        sum = sum +g*rel
    self.assertTrue(sum in ideal)
    
  def test_create_relations_satisfy(self):
    x = self.x
    y = self.y
    z = self.z
    zero = self.poly_ring.zero()
    relations = [[x**2,z*y],[y**2,-x]]
    ideals = [Ideal(self.poly_ring,[x**2*y-z**2]),Ideal(self.poly_ring,[x*y-z])]
    mod = SingularModule.create_from_relations(relations,ideals)
    for relation,ideal in zip(relations,ideals):
      for gen in mod.gens:
        sum = zero
        for g,rel in zip(gen,relation):
          sum = sum +g*rel
      self.assertTrue(sum in ideal)
      
  def  test_create_singular_free(self):
    free_matrix = "MM[1,1]=1\nMM[1,2]=0\nMM[1,3]=0\n"
    free_matrix = free_matrix+"MM[2,1]=0\nMM[2,2]=1\nMM[2,3]=0\n"
    free_matrix = free_matrix+"MM[3,1]=0\nMM[3,2]=0\nMM[3,3]=1\n"
    free_c = SingularModule.create_from_singular_matrix(self.poly_ring,free_matrix)
    free = SingularModule.create_free_module(3,self.poly_ring)
    self.assertTrue(free.equals(free_c))
    
  def test_create_singular(self):
    out = "mat_inter[1,1]=0\n"
    out = out + "mat_inter[1,2]=0\n"
    out = out + "mat_inter[1,3]=x(1)*x(2)*x(3)\n"
    out = out + "mat_inter[2,1]=x(2)\n"
    out = out + "mat_inter[2,2]=0\n"
    out = out + "mat_inter[2,3]=0\n"
    out = out + "mat_inter[3,1]=-x(3)\n"
    out = out + "mat_inter[3,2]=x(3)\n"
    out = out + "mat_inter[3,3]=0\n"
    mod = SingularModule.create_from_singular_matrix(self.poly_ring,out,"mat_inter")
    x = self.x
    y = self.y
    z = self.z
    zero = self.poly_ring.zero()
    mod_true= [[zero,y,-z],[zero,zero,z],[x*y*z,zero,zero]]
    self.assertEqual(mod.gens,mod_true)

  def test_lift_zero(self):
    x = self.x
    y = self.y
    z = self.z
    one = self.poly_ring.one()
    zero = self.poly_ring.zero()
    sm1 = SingularModule([[x,y+z,z**3-2*y],[x,y,z],[x,y,z**2],[one,x,y]])
    vec = sm1.lift([zero,zero,zero])
    self.assertEqual(vec,[zero for _ in range(4)])

  def test_lift_satisfy(self):
    x = self.x
    y = self.y
    z = self.z
    one = self.poly_ring.one()
    zero = self.poly_ring.zero()
    sm1 = SingularModule([[x,y+z,z**3-2*y],[x,y,z]])
    vec = sm1.lift([x**2+x*z+x,x*y+x*z+y*z+y,z**3*x-2*y*x+z**2+z])
    self.assertEqual(vec[0]*x+vec[1]*x,x**2+x*z+x)
    self.assertEqual(vec[0]*(y+z)+vec[1]*y,x*y+x*z+y*z+y)
    self.assertEqual(vec[0]*(z**3-2*y)+vec[1]*z,z**3*x-2*y*x+z**2+z)

  def test_lift_linear(self):
    x = self.x
    y = self.y
    z = self.z
    one = self.poly_ring.one()
    zero = self.poly_ring.zero()
    sm1 = SingularModule([[x,y+z,z**3-2*y],[x,y,z]])
    vec = sm1.lift([3*x,3*y+z,z**3-2*y+2*z],True)
    self.assertEqual(vec,[1,2])

if __name__=="__main__":
  unittest.main()
    
