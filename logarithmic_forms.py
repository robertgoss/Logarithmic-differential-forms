#!/home/robert/sage-4.7.alpha2/sage-python

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

#TODO - log forms class

from singular_module import SingularModule

class NotImplementedException(Exception):
  pass

class NotWieghtHomogeneousException(Exception):
  pass
  
class SymbolicNotPolynomialException(Exception):
  pass
  
  
def homogenous_wieghts(divisor):
  hw = []
  milp = MixedIntegerLinearProgram(maximization=False)
  #Use cts to see if solvable - sometimes using intergers in an unsolvable causes a unhandled sigfault
  wieghts = milp.new_variable(real=True)
  milp.add_constraint(wieghts[0]>=1)
  for i in range(divisor.parent().ngens()):
    milp.add_constraint(wieghts[i+1]>=1)
  for ex in divisor.exponents():
    rel = -wieghts[0]
    for ex_i,ex_m in enumerate(ex):
      rel = rel + ex_m*wieghts[ex_i+1]
    milp.add_constraint(rel==0)
  try:
    milp.solve()
    #Redo with integers to get actual solution
    milp_i = MixedIntegerLinearProgram(maximization=False)
    wieghts_i = milp_i.new_variable(integer=True)
    milp_i.add_constraint(wieghts_i[0]>=1)
    for i in range(divisor.parent().ngens()):
      milp_i.add_constraint(wieghts_i[i+1]>=1)
    for ex in divisor.exponents():
      rel = -wieghts_i[0]
      for ex_i,ex_m in enumerate(ex):
        rel = rel + ex_m*wieghts_i[ex_i+1]
      milp_i.add_constraint(rel==0)
    milp_i.solve()
    for i,v in milp_i.get_values(wieghts_i).iteritems():
      hw.append(v)
  except:
    raise NotWieghtHomogeneousException
  return hw
  
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
  
def convert_symbolic_to_polynomial(symbolic_poly,poly_ring):
  try:
    c = Rational(symbolic_poly)
    return c*poly_ring.one()
  except:
    pass
  base_ring = poly_ring.base_ring()
  poly = symbolic_poly.polynomial(base_ring)
  coeff = poly.coefficients()
  expon = poly.exponents()
  final_poly = poly_ring.zero()
  for c,e in zip(coeff,expon):
    mon = c*poly_ring.one()
    for e_i,e_m in enumerate(e):
      mon = mon * poly_ring.gens()[e_i]**e_m
    final_poly = final_poly + mon
  return final_poly
  
def _log_1_form_rels(divisor):
  poly_ring = divisor.parent()
  n = poly_ring.ngens()
  rels = []
  for i in range(n):
    for j in range(i+1,n):
      rel = []
      for k in range(n):
        if k==i:
          rel.append(divisor.derivative(poly_ring.gens()[i]))
        if k==j:
          rel.append(-divisor.derivative(poly_ring.gens()[j]))
        if k!=i and k!=j:
          rel.append(poly_ring.zero())
      rels.append(rel)
  return rels
  
def _make_poly_1_form(polys,differential_forms,sym_vars):
  form = DifferentialForm(differential_forms,1)
  for i,poly in enumerate(polys):
    form[i] = convert_polynomial_to_symbolic(poly,sym_vars)
  return form
  
def skew_iter(n,depth=1,start=0):
  if depth == 1:
    for i in range(start,n):
      yield [i]
  else:
    for i in range(start,n):
      for v in skew_iter(n,depth-1,start+i+1):
        yield [i]+v
      
class LogarithmicDifferentialForms(SageObject):
  def __init__(self,divisor,var_name="z"):
    self.divisor = divisor
    self.poly_ring = divisor.parent()
    hw = homogenous_wieghts(divisor)
    self.wieghts = hw[1:]
    self.degree = hw[0]
    #Setup patch for differential form
    var_names = [ var_name + str(i) for i in range(self.poly_ring.ngens())]
    self.form_vars = var(",".join(var_names))
    self.form_patch = CoordinatePatch(self.form_vars)
    self.form_space = DifferentialForms(self.form_patch)
    #compute the generators of the logarithmic p-forms
    self._p_modules = {} # modules of logarithmic differential p-forms
    self._p_gens = {} # the generators of the module of logarithmic differential p-forms
    
  def _compute_1_form_generators(self):
    rels = _log_1_form_rels(self.divisor)
    ideals = [[self.divisor]*self.poly_ring for _ in range(len(rels))]
    self._p_modules[1] = SingularModule.create_from_relations(rels,ideals);
    
  def _compute_p_form_generators(self,p):
    if p==0:
      self._p_modules[0] = SingularModule([[self.divisor]])
      return
    n = self.poly_ring.ngens()
    if p==n:
      self._p_modules[n] = SingularModule([[self.poly_ring.one()]])
      return
    if p==1:
      self._compute_1_form_generators()
    else:
      if not 1 in self._p_modules.keys():
        self._compute_1_form_generators()
      #compute wedges of 1 forms
      gens_1_forms = self._p_modules[1].gens
      diff_1_forms = [_make_poly_1_form(polys,self.form_space,self.form_vars) for polys in gens_1_forms]
      diff_p_forms = []
      for s in Set(range(len(gens_1_forms))).subsets(p):
        p_form = DifferentialForm(self.form_space,0,1)
        for i in s:
          p_form = p_form.wedge(diff_1_forms[i])
        diff_p_forms.append(p_form)
      gens_p_forms = []
      for p_form in diff_p_forms:
        gen_p_form = []
        for v in skew_iter(self.poly_ring.ngens(),p):
          poly = convert_symbolic_to_polynomial(p_form[tuple(v)],self.poly_ring)
          gen_p_form.append(poly)
        gens_p_forms.append(gen_p_form)
      self._p_modules[p] = SingularModule(gens_p_forms)
    
  def p_form_generators(self,p):
    #the generators of the module of logarithmic differential p-forms
    if p in self._p_gens.keys():
      return self._p_gens[p]
    if not p in self._p_modules.keys():
      self._compute_p_form_generators(p)
    if p==0:
      zero_form = DifferentialForm(self.form_space,0,1);
      self._p_gens[0] = [zero_form]
      return self._p_gens[0]
    n = self.poly_ring.ngens()
    if p==n:
      top_form = DifferentialForm(self.form_space,n);
      sym_div = convert_polynomial_to_symbolic(self.divisor,self.form_vars)
      top_form[tuple(range(n))] = 1/sym_div
      self._p_gens[n] = [top_form]
      return self._p_gens[n]
    self._p_gens[p] = []
    for gen in self._p_modules[p].gens:
      p_form = DifferentialForm(self.form_space,p);
      for i,v in enumerate(skew_iter(self.poly_ring.ngens(),p)):
        p_form[tuple(v)] = gen[i]/self.divisor
      self._p_gens[p].append(p_form)
    return self._p_gens[p]
    
if __name__=="__main__":
  C = PolynomialRing(QQ,"x,y,z")
  x = C.gens()[0]
  y = C.gens()[1]
  z = C.gens()[2]
  logdf = LogarithmicDifferentialForms(x*y*z)
  print logdf.p_form_generators(0)
  print logdf.p_form_generators(1)
  print logdf.p_form_generators(2)
  print logdf.p_form_generators(3)
