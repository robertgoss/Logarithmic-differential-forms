
import sage.all

from sage.structure.sage_object import SageObject


from sage.rings.rational import Rational

from sage.tensor.differential_form_element import DifferentialForm

def skew_iter(ngens,degree=1,start=0):
  if degree == 1:
    for i in range(start,ngens):
      yield [i]
  else:
    if degree == ngens-1:
      for i in range(ngens-1,start-1,-1):
        yield range(i) + range(i+1,ngens)
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
    
def _get_sym_from_0_form(form,diff_forms):
  dz = DifferentialForm(diff_forms.form_space,1)
  dz[(0,)] = 1
  return form.wedge(dz)[(0,)]

class LogarithmicDifferentialForm(SageObject):
    
  def __init__(self,degree,vec,differential_forms):
    self.vec = vec
    self.degree = degree
    self.diff_forms = differential_forms
    self.divisor = self.diff_forms.divisor
    sym_divisor = convert_polynomial_to_symbolic(self.divisor,self.diff_forms.form_vars)
    #Construct the p_forms version of self
    if degree==0:
      sym_poly = convert_polynomial_to_symbolic(self.vec[0],self.diff_forms.form_vars)
      self.form = DifferentialForm(self.diff_forms.form_space,degree,sym_poly/sym_divisor)
      return
    if degree==self.diff_forms.poly_ring.ngens():
      sym_poly = convert_polynomial_to_symbolic(self.vec[0],self.diff_forms.form_vars)
      self.form = DifferentialForm(self.diff_forms.form_space,degree)
      self.form[tuple(range(degree))] = sym_poly/sym_divisor
    else:
      self.form = DifferentialForm(self.diff_forms.form_space,degree)
      for i,v in enumerate(skew_iter(self.diff_forms.poly_ring.ngens(),degree)):
        sym_poly = convert_polynomial_to_symbolic(self.vec[i],self.diff_forms.form_vars)
        self.form[tuple(v)] = sym_poly/sym_divisor;

  def equals(self,other):
    if self.diff_forms==other.diff_forms:
      if self.degree==other.degree:
        for self_v,other_v in zip(self.vec,other.vec):
          if not self_v==other_v:
            return False
          return True
    return False


  def wedge(self,other):
    diff_form = DifferentialForm.wedge(self.form,other.form)
    return LogarithmicDifferentialForm.create_from_form(diff_form,self.diff_forms)

  def derivative(self):
    diff_form = self.form.derivative()
    return LogarithmicDifferentialForm.create_from_form(diff_form,self.diff_forms)

  def __add__(self,other):
    diff_form = self.form+other.form
    return LogarithmicDifferentialForm.create_from_form(diff_form,self.diff_forms)

  def __sub__(self,other):
    diff_form = self.form-other.form
    return LogarithmicDifferentialForm.create_from_form(diff_form,self.diff_forms)

  def __mul__(self,scalar):
    vec = []
    for v in self.vec:
      vec.append(v*Rational(scalar))
    return LogarithmicDifferentialForm(self.form.degree(),vec,self.diff_forms)

  def __rmul__(self,scalar):
    vec = []
    for v in self.vec:
      vec.append(v*Rational(scalar))
    return LogarithmicDifferentialForm(self.form.degree(),vec,self.diff_forms)

  def interior_product(self):
    if self.degree==0:
      return LogarithmicDifferentialForm.make_zero(0,self.diff_forms)
    prod_form = DifferentialForm(self.diff_forms.form_space,self.degree-1)
    for v in skew_iter(self.diff_forms.poly_ring.ngens(),self.degree):
      for e_i,e in enumerate(v):
        partial_var = self.diff_forms.form_vars[e]
        partial_w = self.diff_forms.wieghts[e]
        partial = ((-1)**(e_i)) * self.form[tuple(v)] * partial_var * partial_w
        if self.degree>1:
          rest = [ e_other for e_other in v if e_other != e]
          prod_form[tuple(rest)] = prod_form[tuple(rest)] + partial
        else:
          prod_form = prod_form + partial
    return LogarithmicDifferentialForm.create_from_form(prod_form,self.diff_forms)

  def __repr__(self):
    return self.form.__repr__()

  def __latex__(self):
    if self.degree==0:
      return "$"+self.form.__repr__()+"$"
    if self.degree == self.diff_forms.poly_ring.ngens():
      string =  "$"+self.form[tuple(range(self.degree))].__repr__()
      diff = []
      for v in self.diff_forms.form_vars:
        diff.append("d "+v.__repr__())
      return string+("\wedge ".join(diff))+"$"
    parts = []
    for i,v in enumerate(skew_iter(self.diff_forms.poly_ring.ngens(),self.degree)):
      if not self.vec[i].is_zero():
        parts.append("$"+self.form[tuple(v)].__repr__())
        diff = []
        for e_i in v:
          diff.append("d "+self.diff_forms.form_vars[e_i].__repr__())
        parts[-1] = parts[-1]+("\wedge ".join(diff)) + "$"
    return "+".join(parts)

  @classmethod
  def create_from_form(cls,form,diff_forms):
    sym_divisor = convert_polynomial_to_symbolic(diff_forms.divisor,diff_forms.form_vars)
    if form.degree()==0:
      sym_poly = _get_sym_from_0_form(form,diff_forms)*sym_divisor
      poly = convert_symbolic_to_polynomial(sym_poly,diff_forms.poly_ring,diff_forms.form_vars)
      return LogarithmicDifferentialForm(0,[poly],diff_forms)
    n = diff_forms.poly_ring.ngens()
    if form.degree()==n:
      sym_poly = form[tuple(range(n))]*sym_divisor
      poly = convert_symbolic_to_polynomial(sym_poly,diff_forms.poly_ring,diff_forms.form_vars)
      return LogarithmicDifferentialForm(n,[poly],diff_forms)
    #Compute vec
    vec = []
    for _,v in enumerate(skew_iter(diff_forms.poly_ring.ngens(),form.degree())):
      sym_poly = form[tuple(v)] * sym_divisor
      vec.append(convert_symbolic_to_polynomial(sym_poly,diff_forms.poly_ring,diff_forms.form_vars))
    return LogarithmicDifferentialForm(form.degree(),vec,diff_forms)
  @classmethod
  def make_unit(cls,diff_forms):
    unit = LogarithmicDifferentialForm(0,[diff_forms.divisor],diff_forms)
    return unit

  @classmethod
  def make_zero(cls,p,diff_forms):
    zero_form = DifferentialForm(diff_forms.form_space,p)
    return LogarithmicDifferentialForm.create_from_form(zero_form,diff_forms)

    
