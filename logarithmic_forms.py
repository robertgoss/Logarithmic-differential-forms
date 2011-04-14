
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

import singular_module
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
  hw = singular_module.wieghts(divisor)
  for w in hw:
    if w!=0:
      degree =  sum([w*e_m for w,e_m in zip(hw,divisor.exponents()[0])])
      for ex in divisor.exponents():
        if degree != sum([w*e_m for w,e_m in zip(hw,ex)]):
          raise NotWieghtHomogeneousException
      return [degree]+hw
  raise NotWieghtHomogeneousException
  
  
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

def lift_to_basis(form,basis,make_linear=False):
  b_mod = SingularModule([b.vec for b in basis])
  return b_mod.lift(form.vec,make_linear)
  
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
  def __init__(self,divisor):
    self.divisor = divisor
    self.poly_ring = divisor.parent()
    hw = homogenous_wieghts(divisor)
    self.wieghts = hw[1:]
    self.degree = hw[0]
    #Setup patch for differential form
    var_names = [ g.__repr__() for g in self.poly_ring.gens()]
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
  
  def _latex_p_part_0_basis(self,p):
    string = "\subsubsection{$\Omega_0^"+str(p)+"(log D)$}"
    string = string + "\\begin{enumerate}\n"
    for form in self.p_forms_zero_basis(p):
      string = string + "\item "+form.__latex__()
    string = string + "\end{enumerate}"
    return string
  
  def latex_basis(self):
    return "\n".join([self._latex_p_part_0_basis(p) for p in range(self.poly_ring.ngens())])
    
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
        d_p_rows.append(lift_to_basis(d_form,p_forms_1,True))
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
        d_p_1_rows.append(lift_to_basis(d_form,p_forms_2,True))
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

  def _p_equivarient_homology(self,p):
    eqi_p_0_space = [self.p_forms_zero_basis(i) for i in range(p-1,-1,-2)]
    eqi_p_1_space = [self.p_forms_zero_basis(i) for i in range(p,-1,-2)]
    eqi_p_2_space = [self.p_forms_zero_basis(i) for i in range(p+1,-1,-2)]
    eqi_dim_0 = sum([len(b) for b in eqi_p_0_space])
    eqi_dim_1 = sum([len(b) for b in eqi_p_1_space])
    eqi_dim_2 = sum([len(b) for b in eqi_p_2_space])
    if eqi_dim_1==0:
      return []
    if eqi_dim_0==0:
      p_space = VectorSpace(QQ,eqi_dim_1)
      img = p_space.subspace([p_space.zero()])
    else:
      g_mods_1 = [self._p_graded_module(i) for i in range(p,-1,-2)]
      d_p_rows = []
      for level_i,level in enumerate(eqi_p_0_space):
        level_dim = p-(level_i*2)
        for form in level:
          d_form = form.derivative()
          if len(eqi_p_1_space[level_i])!=0:
            d_vec = lift_to_basis(d_form,eqi_p_1_space[level_i],True)
          if level_i+1<len(eqi_p_1_space) and len(eqi_p_1_space[level_i+1])!=0:
            i_form = form.interior_product()
            i_vec = lift_to_basis(i_form,eqi_p_1_space[level_i+1],True)
          row = []
          for lev_i,lev in enumerate(eqi_p_1_space):
            if len(eqi_p_1_space[lev_i])!=0:
              if lev_i == level_i:
                row.extend(d_vec)
              else:
                if lev_i == level_i+1:
                  row.extend(i_vec)
                else:
                  row.extend([self.poly_ring.zero() for _ in eqi_p_1_space[lev_i]])
          d_p_rows.append(row)
      mat_p = matrix(QQ,d_p_rows)
      img = mat_p.image() # Sage computes the image by left multiplication!    
    if eqi_dim_2==0:
      p_space = VectorSpace(QQ,eqi_dim_1)
      ker = p_space
    else:
      g_mods_2 = [self._p_graded_module(i) for i in range(p+1,-1,-2)]
      d_p_1_rows = []
      for level_i,level in enumerate(eqi_p_1_space):
        for form in level:
          d_form = form.derivative()
          i_form = form.interior_product()
          if not len(eqi_p_2_space[level_i])==0:
            d_vec = lift_to_basis(d_form,eqi_p_2_space[level_i],True)
          if level_i+1<len(eqi_p_2_space) and len(eqi_p_2_space[level_i+1])!=0:
            i_form = form.interior_product()
            i_vec = lift_to_basis(i_form,eqi_p_2_space[level_i+1],True)
          row = []
          for lev_i,lev in enumerate(eqi_p_2_space):
            if len(eqi_p_2_space[lev_i])!=0:
              if lev_i == level_i:
                row.extend(d_vec)
              else:
                if lev_i == level_i+1:
                  row.extend(i_vec)
                else:
                  row.extend([self.poly_ring.zero() for _ in eqi_p_2_space[lev_i]])
          d_p_1_rows.append(row)
      mat_p_1 = matrix(QQ,d_p_1_rows).transpose()
      ker = mat_p_1.right_kernel()   
    hom = ker.quotient(img)
    return hom.rank()

  def equivarient_homology(self):
    homology = {}
    homology[0] = 1
    for i in range(1,self.poly_ring.ngens()+1):
      homology[i] = self._p_equivarient_homology(i)
    return homology

  def complement_homology(self):
    homology = {}
    homology[0] = [LogarithmicDifferentialForm.make_unit(self)]
    for i in range(1,self.poly_ring.ngens()+1):
      homology[i] = self._p_complement_homology(i)
    return homology

  def complement_homology_latex(self):
    string = "\\begin{description}\n"
    c_h = self.complement_homology()
    for i in range(self.poly_ring.ngens()+1):
      string = string+"\item[$H^"+str(i)+"(log D)$] rank: $"+str(len(c_h[i]))+"$\n"
      string = string + "\\begin{enumerate}\n"
      for b in c_h[i]:
        string = string +"\item "+b.__latex__()+"\n"
      string = string + "\end{enumerate}"
    return string + "\end{description}\n"
      
