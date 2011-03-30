#!/home/robert/sage-4.7.alpha2/sage-python

import unittest

from logarithmic_forms import SingularModule
from logarithmic_forms import homogenous_wieghts
from logarithmic_forms import NotWieghtHomogeneousException


from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

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
    self.assertEqual(sm.create_module_str("MT"),"module MT=[1*x(1)^1,2*x(2)^2,3*x(3)^3];\n")
    
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
      
  def test_equals(self):
    x = self.x
    y = self.y
    z = self.z
    #Assert these are equal
    sm1 = SingularModule([[x,y+z,z**3-2*y],[x,y,z],[x,y,z**2]])
    sm2 = SingularModule([[x,y,z],[x,y+z,z**3-2*y],[x,y,z**2]])
    self.assertTrue(sm1.equals(sm2))
    
  def test_equals_not(self):
    x = self.x
    y = self.y
    z = self.z
    sm1 = SingularModule([[x,y+z,z**3-2*y],[x,y,z],[x,y,z**2]])
    sm2 = SingularModule([[x,y,z]])
    #Assert these are not equal
    self.assertFalse(sm1.equals(sm2))
    
    
    
class Test_Homogeneous_Wieghts(unittest.TestCase):

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
    self.assertEquals(wieghts,[1,1,1]);
    
  def test_weighted_homogenous(self):
    x = self.x
    y = self.y
    z = self.z
    divisor = x**2*y-z**2
    wieghts = homogenous_wieghts(divisor)
    #Test this works with a weighted homogenous divisor
    self.assertEquals(wieghts,[1,2,2])
  
  def test_not_homogenous(self):
    x = self.x
    y = self.y
    z = self.z
    divisor = x**2*y-x**3 + z +y**2;
    self.assertRaises(NotWieghtHomogeneousException,homogenous_wieghts,divisor)
    
if __name__=="__main__":
  unittest.main()
