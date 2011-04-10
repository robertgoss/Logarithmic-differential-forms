#!/usr/local/bin/sage

import unittest

from logarithmic_form import convert_polynomial_to_symbolic
from logarithmic_form import convert_symbolic_to_polynomial
from logarithmic_form import LogarithmicDifferentialForm

from logarithmic_forms import LogarithmicDifferentialForms

from sage.tensor.differential_form_element import DifferentialForm

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

class TestLogarithmicDifferentialForm(unittest.TestCase):
    
  def setUp(self):
    self.poly_ring = PolynomialRing(QQ,"x",3);
    x = self.x = self.poly_ring.gens()[0];
    y = self.y = self.poly_ring.gens()[1];
    z = self.z = self.poly_ring.gens()[2];
    self.normal_logdf = LogarithmicDifferentialForms(x*y*z)
    self.whitney_logdf = LogarithmicDifferentialForms(x**2*y-z**2)

  def test_creationA(self):
    x = self.normal_logdf.form_vars[0]
    y = self.normal_logdf.form_vars[1]
    z = self.normal_logdf.form_vars[2]
    logdf = LogarithmicDifferentialForm(2,[self.x,self.y,self.z],self.normal_logdf)
    form = DifferentialForm(self.normal_logdf.form_space,2)
    form[0,1] = 1/(y*z)
    form[0,2] = 1/(x*z)
    form[1,2] = 1/(x*y)
    self.assertEqual(form,logdf.form)

  def test_creation_B(self):
    x = self.whitney_logdf.form_vars[0]
    y = self.whitney_logdf.form_vars[1]
    z = self.whitney_logdf.form_vars[2]
    whitney = x**2*y-z**2
    logdf = LogarithmicDifferentialForm(2,[self.x,self.y,self.z],self.whitney_logdf)
    form = DifferentialForm(self.normal_logdf.form_space,2)
    form[0,1] = x/whitney
    form[0,2] = y/whitney
    form[1,2] = z/whitney
    self.assertEqual(form,logdf.form)

  def test_creation_0_form(self):
    x = self.normal_logdf.form_vars[0]
    y = self.normal_logdf.form_vars[1]
    z = self.normal_logdf.form_vars[2]
    logdf = LogarithmicDifferentialForm(0,[self.x*self.y],self.normal_logdf)
    form = DifferentialForm(self.normal_logdf.form_space,0,1/z)
    self.assertEqual(form,logdf.form)

  def test_add(self):
    x = self.normal_logdf.form_vars[0]
    y = self.normal_logdf.form_vars[1]
    z = self.normal_logdf.form_vars[2]
    logdfA = LogarithmicDifferentialForm(2,[self.x,self.y,self.z],self.normal_logdf)
    logdfB = LogarithmicDifferentialForm(2,[self.z,self.y,self.x],self.normal_logdf)
    form = DifferentialForm(self.normal_logdf.form_space,2)
    form[0,1] = 1/(y*z) + 1/(x*y)
    form[0,2] = 2/(x*z)
    form[1,2] = 1/(x*y) + 1/(y*z)
    logdf_sum = LogarithmicDifferentialForm.create_from_form(form,self.normal_logdf)
    self.assertTrue(logdf_sum.equals(logdfA+logdfB))

  def test_sub(self):
    x = self.normal_logdf.form_vars[0]
    y = self.normal_logdf.form_vars[1]
    z = self.normal_logdf.form_vars[2]
    logdfA = LogarithmicDifferentialForm(2,[self.x,self.y,self.z],self.normal_logdf)
    logdfB = LogarithmicDifferentialForm(2,[self.z,self.y,self.x],self.normal_logdf)
    form = DifferentialForm(self.normal_logdf.form_space,2)
    form[0,1] = 1/(y*z) - 1/(x*y)
    form[0,2] = 0
    form[1,2] = 1/(x*y) - 1/(y*z)
    logdf_diff = LogarithmicDifferentialForm.create_from_form(form,self.normal_logdf)
    self.assertTrue(logdf_diff.equals(logdfA-logdfB))

  def test_mul(self):
    x = self.normal_logdf.form_vars[0]
    y = self.normal_logdf.form_vars[1]
    z = self.normal_logdf.form_vars[2]
    logdf = LogarithmicDifferentialForm(2,[self.x,self.y,self.z],self.normal_logdf)
    form = DifferentialForm(self.normal_logdf.form_space,2)
    form[0,1] = -5/(y*z)
    form[0,2] = -5/(x*z)
    form[1,2] = -5/(x*y)
    self.assertTrue((logdf*(-5)).equals((-5)*logdf))
    logdf_mul = LogarithmicDifferentialForm.create_from_form(form,self.normal_logdf)
    self.assertTrue((-5*logdf).equals(logdf_mul))

  def test_wedge(self):
    x = self.normal_logdf.form_vars[0]
    y = self.normal_logdf.form_vars[1]
    z = self.normal_logdf.form_vars[2]
    norm = self.x*self.y*self.z
    logdfA = LogarithmicDifferentialForm(1,[self.x*norm,self.y*norm,self.z*norm],self.normal_logdf)
    logdfB = LogarithmicDifferentialForm(1,[self.z,self.y,self.x],self.normal_logdf)
    form = DifferentialForm(self.normal_logdf.form_space,2)
    form[0,1] = 1/z - 1/x
    form[0,2] = (x/(y*z)) - (z/(x*y))
    form[1,2] = 1/z - 1/x
    logdf_wedge = LogarithmicDifferentialForm.create_from_form(form,self.normal_logdf)
    self.assertTrue(logdf_wedge.equals(logdfA.wedge(logdfB)))

  def test_derivative(self):
    x = self.x
    y = self.y
    z = self.z
    logdf = LogarithmicDifferentialForm(1,[y*z,x*z,x*y],self.normal_logdf)
    form = DifferentialForm(self.normal_logdf.form_space,2)
    logdf_der = LogarithmicDifferentialForm.create_from_form(form,self.normal_logdf)
    self.assertTrue(logdf_der.equals(logdf.derivative()))

  def test_unit(self):
    one = LogarithmicDifferentialForm.make_unit(self.normal_logdf)
    x = self.normal_logdf.form_vars[0]
    y = self.normal_logdf.form_vars[1]
    z = self.normal_logdf.form_vars[2]
    form = DifferentialForm(self.normal_logdf.form_space,0,1)
    one_form = LogarithmicDifferentialForm.create_from_form(form,self.normal_logdf)
    self.assertTrue(one.equals(one_form))

  def test_zero_p(self):
    size = [1,3,3,1]
    for s,i in zip(size,range(4)):
      zero = LogarithmicDifferentialForm.make_zero(i,self.normal_logdf)
      zero_vec = [self.normal_logdf.poly_ring.zero() for _ in range(s)]
      self.assertEqual(zero.vec,zero_vec)

  def test_create_from_form(self):
    x = self.normal_logdf.form_vars[0]
    y = self.normal_logdf.form_vars[1]
    z = self.normal_logdf.form_vars[2]
    xp = self.x
    yp = self.y
    zp = self.z
    logdf = LogarithmicDifferentialForm(2,[xp*yp-yp*zp,xp**2-zp**2,xp*yp-yp*zp],self.normal_logdf)
    form = DifferentialForm(self.normal_logdf.form_space,2)
    form[0,1] = 1/z - 1/x
    form[0,2] = (x/(y*z)) - (z/(x*y))
    form[1,2] = 1/z - 1/x
    logdf_form = LogarithmicDifferentialForm.create_from_form(form,self.normal_logdf)
    self.assertTrue(logdf.equals(logdf_form))

  def test_create_from_0_form(self):
    x = self.normal_logdf.form_vars[0]
    y = self.normal_logdf.form_vars[1]
    z = self.normal_logdf.form_vars[2]
    xp = self.x
    yp = self.y
    zp = self.z
    logdf = LogarithmicDifferentialForm(0,[xp+yp+zp],self.normal_logdf)
    form = DifferentialForm(self.normal_logdf.form_space,0,(x+y+z)/(x*y*z))
    logdf_form = LogarithmicDifferentialForm.create_from_form(form,self.normal_logdf)
    self.assertTrue(logdf.equals(logdf_form))


  def test_interior_0_form(self):
    x = self.x
    y = self.y
    z = self.z
    logdf = LogarithmicDifferentialForm(0,[x*y-y*z+x**2-z**2+x*x-y*z],self.normal_logdf)
    int_product = logdf.interior_product()
    zero = LogarithmicDifferentialForm.make_zero(0,self.normal_logdf)
    self.assertTrue(int_product.equals(zero))

  def test_interior_1_form(self):
    x = self.normal_logdf.form_vars[0]
    y = self.normal_logdf.form_vars[1]
    z = self.normal_logdf.form_vars[2]
    form = DifferentialForm(self.normal_logdf.form_space,0,(x+y+z)/(x*y*z))
    one = self.poly_ring.one()
    logdf = LogarithmicDifferentialForm(1,[one,one,one],self.normal_logdf)
    int_product = logdf.interior_product()
    logdf_form = LogarithmicDifferentialForm.create_from_form(form,self.normal_logdf)
    self.assertTrue(int_product.equals(logdf_form))

  def test_interior(self):
    xp = self.x
    yp = self.y
    zp = self.z
    logdf = LogarithmicDifferentialForm(2,[xp**2,xp+4*yp+zp,xp**3*zp**3],self.whitney_logdf)
    int_product = logdf.interior_product()
    x = self.normal_logdf.form_vars[0]
    y = self.normal_logdf.form_vars[1]
    z = self.normal_logdf.form_vars[2]
    whitney = x**2*y - z**2
    form = DifferentialForm(self.normal_logdf.form_space,1)
    form[0] = (-2*x**2*y-2*x*z-8*y*z-2*z**2)/whitney
    form[1] = (x**3-2*x**3*z**4)/whitney
    form[2] = (x**2+4*x*y+x*z+2*x**3*y*z**3)/whitney
    logdf_form = LogarithmicDifferentialForm.create_from_form(form,self.whitney_logdf)
    self.assertTrue(int_product.equals(logdf_form))

if __name__=="__main__":
    unittest.main()
