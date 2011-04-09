
import sage.all

from sage.structure.sage_object import SageObject

from sage.numerical.mip import MixedIntegerLinearProgram

from sage.sets.set import Set

from sage.rings.rational import Rational

from sage.tensor.coordinate_patch import CoordinatePatch
from sage.tensor.differential_forms import DifferentialForms
from sage.tensor.differential_form_element import DifferentialForm

from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.symbolic.ring import var
from sage.modules.free_module_element import vector
from sage.modules.free_module import VectorSpace
from sage.matrix.constructor import matrix

def skew_iter(ngens,degree=1,start=0):
  if degree == 1:
    for i in range(start,ngens):
      yield [i]
  else:
    for i in range(start,ngens):
      for v in skew_iter(ngens,degree-1,start+i+1):
        yield [i]+v

def convert_polynomial_to_symbolic(poly,sym_vars):
  sym_poly = 0
  coeff = poly.coefficients()
  expon = poly.exponents()
  for c,e in zip(coeff,expon):
    mon = c
    for e_i,e_m in enumerate(e):
      mon = mon * sym_vars[e_i]**e_m
    sym_poly = sym_poly + mon
  return sym_poly
  
def convert_symbolic_to_polynomial(symbolic_poly,poly_ring,sym_vars):
  try:
    c = Rational(symbolic_poly)
    return c*poly_ring.one()
  except:
    pass
  base_ring = poly_ring.base_ring()
  poly = symbolic_poly.polynomial(base_ring)
  p_vars = symbolic_poly.variables()
  var_map = {}
  for v in p_vars:
    for i,s_v in enumerate(sym_vars):
      if v==s_v:
        var_map[v] = poly_ring.gens()[i]
  final_poly = poly_ring.zero()
  if len(p_vars)==1:
    for e,c in zip(poly.exponents(),poly.coefficients()):
      final_poly = final_poly + c*var_map[p_vars[0]]**e
  else:
    for e,c in zip(poly.exponents(),poly.coefficients()):
      mon = c*poly_ring.one()
      for e_i,e_m in enumerate(e):
        mon = mon * var_map[p_vars[e_i]]**e_m
      final_poly = final_poly + mon
  return final_poly
    

class LogarithmicDifferentialForm(SageObject):
    
  def __init__(self,degree,vec,differential_forms):
    self.vec = vec
    self.diff_forms = differential_forms
    self.divisor = self.diff_forms.divisor
    sym_divisor = convert_polynomial_to_symbolic(self.divisor,self.diff_forms.form_vars)
    #Construct the p_forms version of self
    if degree==0:
      self.form = DifferentialForm(self.diff_forms.form_space,degree,self.diff_forms.poly_ring.one())
    else:
      self.form = DifferentialForm(self.diff_forms.form_space,degree)
      for i,v in enumerate(skew_iter(self.diff_forms.poly_ring.ngens(),degree)):
        sym_poly = convert_polynomial_to_symbolic(self.vec[i],self.diff_forms.form_vars)
        self.form[tuple(v)] = sym_poly/sym_divisor;

  def wedge(self,other):
    diff_form = DifferentialForm.wedge(self.form,other.form)
    total_degree = diff_form.degree()
    #Compute vec
    vec = []
    sym_divisor = convert_polynomial_to_symbolic(self.divisor,self.diff_forms.form_vars)
    for i,v in enumerate(skew_iter(self.diff_forms.poly_ring.ngens(),total_degree)):
      sym_poly = diff_form[tuple(v)] * sym_divisor
      vec.append(convert_symbolic_to_polynomial(sym_poly,self.diff_forms.poly_ring,self.diff_forms.form_vars))
    return LogarithmicDifferentialForm(total_degree,vec,self.diff_forms)

  @classmethod
  def make_unit(cls,diff_forms):
    unit = LogarithmicDifferentialForm(0,[diff_forms.divisor],diff_forms)
    return unit

    
