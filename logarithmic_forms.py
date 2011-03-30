#!/home/robert/sage-4.7.alpha2/sage-python

import sage.all

from sage.structure.sage_object import SageObject

from singular_module import SingularModule

from sage.numerical.mip import MixedIntegerLinearProgram

from sage.symbolic.ring import var

from sage.tensor.coordinate_patch import CoordinatePatch
from sage.tensor.differential_forms import DifferentialForms

#TODO - log forms class

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
    for j in range(i,n):
      rel = []
      for k in range(n):
        if k==i:
          rel.append(divisor.derivative(poly_ring.gens()[j]))
        if k==j
          rel.append(divisor.derivative(-poly_ring.gens()[i]))
        if k!=i and k!=j:
          rel.append(poly_ring.zero())
      rels.append(rel)
  return rels
    
class LogarithmicDifferentialForms(SageObject):
  def __init__(self,divisor,var_name="z"):
    self.divisor = divisor
    self.poly_ring = divisor.parent()
    hw = homogenous_wieghts(divisor)
    self.wieghts = hw[1:]
    self.degree = hw[0]
    #Setup patch for differential form
    var_names = [ var_name + str(i) for i in range(self.poly_ring.ngens())]
    self.form_vars = var(join(var_names,","))
    self.form_patch = CoordinatePatch(self.form_vars)
    self.form_space = DifferentialForms(self.form_patch)
    #compute the generators of the logarithmic p-forms
    self._p_modules = {} # modules of logarithmic differential p-forms
    self._p_gens = {} # the generators of the module of logarithmic differential p-forms
    
  def _compute_1_form_generators(self):
    rels = _log_1_form_rels(self.divisor)
    self._p_modules[1] = SingularModule.create_from_relations(rels);
    
  def _compute_p_form_generators(self,p):
    if p==1:
      self._compute_1_form_generators(self)
    else:
      raise NotImplementedException();
    
  def p_form_generators(self,p):
    #the generators of the module of logarithmic differential p-forms
    if p in self._p_gens.keys():
      return self._p_gens[p]
    if not p in self._p_modules.keys():
      _compute_p_form_generators(p)
    self._p_gens[p] = []
    for gen in self._p_modules[p].gens:
      pass;