
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
from graded_module import GradedModule

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
          rel.append(divisor.derivative(poly_ring.gens()[j]))
        if k==j:
          rel.append(-divisor.derivative(poly_ring.gens()[i]))
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
    self._p_zero_part_basis = {}

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
        p_vec = self._convert_p_form_to_p_vec(p,p_form)
        gens_p_forms.append(p_vec)
      #Normaize out - important
      gens = []
      for gen in gens_p_forms:
        norm_gen = [g//(self.divisor**(p-1)) for g in gen]
        gens.append(norm_gen)
      self._p_modules[p] = SingularModule(gens)

  def _convert_p_vec_to_p_form(self,p,vec):
    p_form = DifferentialForm(self.form_space,p);
    for i,v in enumerate(skew_iter(self.poly_ring.ngens(),p)):
      p_form[tuple(v)] = vec[i]/self.divisor
    return p_form
    
  def _convert_p_form_to_p_vec(self,p,p_form):
    p_vec = []
    for i,v in enumerate(skew_iter(self.poly_ring.ngens(),p)):
      poly = convert_symbolic_to_polynomial(p_form[tuple(v)],self.poly_ring)
      p_vec.append(poly)
    return p_vec

  def p_form_generators(self,p):
    #the generators of the module of logarithmic differential p-forms
    if p in self._p_gens.keys():
      return self._p_gens[p]
    if not p in self._p_modules.keys():
      self._compute_p_form_generators(p)
    self._p_gens[p] = []
    for gen in self._p_modules[p].gens:
      self._p_gens[p].append(self._convert_p_vec_to_p_form(p,gen))
    return self._p_gens[p]
    
  def p_module(self,p):
    if not p in self._p_modules.keys():
      self._compute_p_form_generators(p)
    return self._p_modules[p]

  def _compute_zero_part_basis(self,p):
    p_mod = self.p_module(p)
    if p==0:
      self._p_zero_part_basis[0] = self.p_module(0).gens
      return
    column_wieghts = []
    for v in skew_iter(self.poly_ring.ngens(),p):
      wsum = 0
      for i in v:
        wsum = wsum + self.wieghts[i]
      column_wieghts.append(wsum)
    if p==self.poly_ring.ngens():
      column_wieghts = [sum(self.wieghts)]
    g_mod = GradedModule(p_mod.gens,column_wieghts,self.wieghts)
    basis = g_mod.homogeneous_part_basis(self.degree)
    self._p_zero_part_basis[p] = basis
    
  def p_forms_zero_basis(self,p):
    if not p in self._p_zero_part_basis.keys():
      self._compute_zero_part_basis(p)
    basis = self._p_zero_part_basis[p]
    basis_forms = []
    for vec in basis:
      basis_forms.append(self._convert_p_vec_to_p_form(p,vec))
    return basis_forms

  def p_module_zero_basis(self,p):
    if not p in self._p_zero_part_basis.keys():
      self._compute_zero_part_basis(p)
    return self._p_zero_part_basis[p]
