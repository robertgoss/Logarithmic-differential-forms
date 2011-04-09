
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
#TODO - log forms class

from singular_module import SingularModule
from graded_module import GradedModule

from logarithmic_form import LogarithmicDifferentialForm
from logarithmic_form import convert_symbolic_to_polynomial
from logarithmic_form import convert_polynomial_to_symbolic

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
    self._p_gens[1] = []
    for g in self._p_modules[1].gens:
      self._p_gens[1].append(LogarithmicDifferentialForm(1,g,self))
    
  def _compute_p_form_generators(self,p):
    if p==0:
      self._p_modules[0] = SingularModule([[self.divisor]])
      self._p_gens[0] = [LogarithmicDifferentialForm(0,[self.divisor],self)]
      return
    n = self.poly_ring.ngens()
    if p==n:
      self._p_modules[n] = SingularModule([[self.poly_ring.one()]])
      self._p_gens[n] = [LogarithmicDifferentialForm(0,[self.poly_ring.one()],self)]
      return
    if p > n:
      return
    if p==1:
      self._compute_1_form_generators()
    else:
      if not 1 in self._p_modules.keys():
        self._compute_1_form_generators()
      #compute wedges of 1 forms
      self._p_gens[p] = []
      for s in Set(range(len(self._p_gens[1]))).subsets(p):
        p_form = LogarithmicDifferentialForm.make_unit(self)
        for i in s:
          p_form = p_form.wedge(self._p_gens[1][i])
        self._p_gens[p].append(p_form)
      gens = []
      for p_form in self._p_gens[p]:
        gens.append(p_form.vec)
      self._p_modules[p] = SingularModule(gens)

  def _convert_p_vec_to_p_form(self,p,vec):
    p_form = DifferentialForm(self.form_space,p);
    div_sym = convert_polynomial_to_symbolic(self.divisor,self.form_vars)
    for i,v in enumerate(skew_iter(self.poly_ring.ngens(),p)):
      vec_sym = convert_polynomial_to_symbolic(vec[i],self.form_vars)
      p_form[tuple(v)] = vec_sym/div_sym
    return p_form
    
  def _convert_p_form_to_p_vec(self,p,p_form):
    p_vec = []
    for i,v in enumerate(skew_iter(self.poly_ring.ngens(),p)):
      poly = convert_symbolic_to_polynomial(p_form[tuple(v)],self.poly_ring,self.form_vars)
      p_vec.append(poly)
    return p_vec

  def p_form_generators(self,p):
    if p > self.poly_ring.ngens():
      return []
    #the generators of the module of logarithmic differential p-forms
    if p in self._p_gens.keys():
      return self._p_gens[p]
    self._compute_p_form_generators(p)
    return self._p_gens[p]
    
  def p_module(self,p):
    if p > self.poly_ring.ngens():
      return None
    if not p in self._p_modules.keys():
      self._compute_p_form_generators(p)
    return self._p_modules[p]

  def _compute_zero_part_basis(self,p):
    if p==0:
      self._p_zero_part_basis[0] = self.p_module(0).gens
      return
    g_mod = self._p_graded_module(p)
    basis = g_mod.homogeneous_part_basis(self.degree)
    self._p_zero_part_basis[p] = basis
    
  def p_forms_zero_basis(self,p):
    if p > self.poly_ring.ngens():
      return []
    if not p in self._p_zero_part_basis.keys():
      self._compute_zero_part_basis(p)
    basis = self._p_zero_part_basis[p]
    basis_forms = []
    for vec in basis:
      basis_forms.append(LogarithmicDifferentialForm(p,vec,self))
    return basis_forms

  def p_module_zero_basis(self,p):
    if p > self.poly_ring.ngens():
      return []
    if not p in self._p_zero_part_basis.keys():
      self._compute_zero_part_basis(p)
    return self._p_zero_part_basis[p]

  def _p_graded_module(self,p):
    if p > self.poly_ring.ngens():
      return None
    p_mod = self.p_module(p)
    column_wieghts = []
    for v in skew_iter(self.poly_ring.ngens(),p):
      wsum = 0
      for i in v:
        wsum = wsum + self.wieghts[i]
      column_wieghts.append(wsum)
    if p==self.poly_ring.ngens():
      column_wieghts = [sum(self.wieghts)]
    return GradedModule(p_mod.gens,column_wieghts,self.wieghts)

  def _form_in_terms_of_basis(self,p,form,basis,gm):
    poly_form = form.vec
    poly_basis = []
    for b in basis:
      poly_basis.append(form.vec)
    lift_basis = []
    for pb in poly_basis:
      lift_basis.append(gm.lift(pb))
    lift_form = gm.lift(poly_form)
    mat = matrix(QQ,lift_basis)
    lift = mat.solve_left(vector(lift_form))
    return lift

  def _p_complement_homology(self,p):
    p_forms_0 = self.p_forms_zero_basis(p-1)
    p_forms_1 = self.p_forms_zero_basis(p)
    p_forms_2 = self.p_forms_zero_basis(p+1)
    if len(p_forms_1)==0:
      #In this case img and ker must  be trivial
      return []
    if len(p_forms_0)==0:
      p_space = VectorSpace(QQ,len(p_forms_1)) 
      img = p_space.subspace([p_space.zero()])
    else:
      g_mod_1 = self._p_graded_module(p)
      d_p_rows = []
      for form in p_forms_0:
        d_form = form.derivative()
        d_p_rows.append(self._form_in_terms_of_basis(p,d_form,p_forms_1,g_mod_1))
      mat_p = matrix(QQ,d_p_rows)
      img = mat_p.image() # Sage computes the image by left multiplication!
    if len(p_forms_2)==0:
      p_space = VectorSpace(QQ,len(p_forms_1))
      ker = p_space
    else:
      g_mod_2 = self._p_graded_module(p+1)
      d_p_1_rows = []
      for form in p_forms_1:
        d_form = form.derivative()
        d_p_1_rows.append(self._form_in_terms_of_basis(p+1,d_form,p_forms_2,g_mod_2))
      mat_p_1 = matrix(QQ,d_p_1_rows).transpose()
      ker = mat_p_1.right_kernel()
    hom = ker.quotient(img)
    hom_forms = []
    for b in hom.basis():
      lift = hom.lift(b)
      form = LogarithmicDifferentialForm.make_zero(p,self)
      for l,f in zip(lift,p_forms_1):
        form = form + (l*f)
      hom_forms.append(form)
    return hom_forms;

  def complement_homology(self):
    homology = {}
    homology[0] = [1]
    for i in range(1,self.poly_ring.ngens()+1):
      homology[i] = self._p_complement_homology(i)
    return homology
