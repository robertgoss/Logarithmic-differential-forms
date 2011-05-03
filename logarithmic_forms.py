
import sage.all

from sage.structure.sage_object import SageObject

from sage.sets.set import Set

from sage.tensor.coordinate_patch import CoordinatePatch
from sage.tensor.differential_forms import DifferentialForms
from sage.tensor.differential_form_element import DifferentialForm

from sage.rings.rational_field import QQ
from sage.symbolic.ring import var
from sage.modules.free_module_element import vector
from sage.modules.free_module import VectorSpace
from sage.matrix.constructor import matrix
#TODO - log forms class

import singular_module
from singular_module import SingularModule
from graded_module import GradedModule

from logarithmic_form import LogarithmicDifferentialForm
from logarithmic_form import skew_iter
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
  
def _weighted_sum(weights,forms,diff_forms):
  part = LogarithmicDifferentialForm.make_zero(forms[0].degree,diff_forms)
  for w,f in zip(weights,forms):
    part = part + w*f
  return part

def orth_complement(space,subspace):
  if subspace.dimension()==0:
    return space
  else:
    ker = subspace.basis_matrix().right_kernel()
    return space.intersection(ker)
  
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
  graded_basis = {}
  poly_ring = basis[0].diff_forms.poly_ring
  lift = [poly_ring.zero() for _ in basis]
  for i,b in enumerate(basis):
    if b.degree in graded_basis.keys():
      graded_basis[b.degree].append((i,b))
    else:
      graded_basis[b.degree] = [(i,b)]
  for f in form:
    if f.degree in graded_basis.keys():
      g_part = graded_basis[f.degree]
      b_mod = SingularModule([b.vec for _,b in g_part])
      res = b_mod.lift(f.vec,make_linear)
      for i,c in enumerate(res):
        lift[g_part[i][0]] = c
  return lift
  
def _make_poly_1_form(polys,differential_forms,sym_vars):
  form = DifferentialForm(differential_forms,1)
  for i,poly in enumerate(polys):
    form[i] = convert_polynomial_to_symbolic(poly,sym_vars)
  return form
      
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
      self._p_modules[p].reduce_generators()

  def _convert_p_vec_to_p_form(self,p,vec):
    p_form = DifferentialForm(self.form_space,p);
    div_sym = convert_polynomial_to_symbolic(self.divisor,self.form_vars)
    for i,v in enumerate(skew_iter(self.poly_ring.ngens(),p)):
      vec_sym = convert_polynomial_to_symbolic(vec[i],self.form_vars)
      p_form[tuple(v)] = vec_sym/div_sym
    return p_form
    
  def _convert_p_form_to_p_vec(self,p,p_form):
    p_vec = []
    for _,v in enumerate(skew_iter(self.poly_ring.ngens(),p)):
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
    if p==0:
      column_wieghts = [0]
    if p==self.poly_ring.ngens():
      column_wieghts = [sum(self.wieghts)]
    return GradedModule(p_mod.gens,column_wieghts,self.wieghts)
  
  def _latex_p_part_0_basis(self,p):
    string = "\subsubsection{$\Omega_0^"+str(p)+"(log D)$}"
    string = string + "\\begin{enumerate}\n"
    zero_part = self._p_graded_module(p).homogeneous_part_basis(self.degree)
    for form in zero_part:
      string = string + "\item "+form.__latex__()
    string = string + "\end{enumerate}"
    return string
  
  def latex_basis(self):
    return "\n".join([self._latex_p_part_0_basis(p) for p in range(self.poly_ring.ngens())])
    
  def _complex_complement(self,n,*args):
    if n<0 or n>self.poly_ring.ngens():
      return []
    basis = self._p_graded_module(n).homogeneous_part_basis(self.degree)
    complex = []
    for b in basis:
      complex.append(LogarithmicDifferentialForm(n,b,self))
    return complex

  def _complex_equivarient(self,n,*args):
    equi_complex = []
    for level in range(n,-1,-2):
      equi_complex.extend(self._complex_complement(level,*args))
    return equi_complex

  def _complex_relative(self,n,*args):
    if len(args)==2:
      complex = []
      for i in range(args[0],args[1]):
        complex.append(self._complex_relative(n,i))
      return complex
    if n<0 or n>self.poly_ring.ngens():
      return []
    if n==0:
      hom_basis = self._p_graded_module(n).homogeneous_part_basis(self.degree+args[0])
      return [LogarithmicDifferentialForm(n,b,self) for b in hom_basis]
    base = self._p_graded_module(n).homogeneous_part_basis(self.degree+args[0])
    if len(base)==0:
      return []
    vs_base = VectorSpace(QQ,len(base))
    df_base = [LogarithmicDifferentialForm(n,b,self) for b in base]
    pre_base = self._p_graded_module(n-1).homogeneous_part_basis(self.degree+args[0])
    if len(pre_base)==0:
      return df_base
    dh = [self.divisor.derivative(g) for g in self.poly_ring.gens()]
    dh = LogarithmicDifferentialForm(1,dh,self)
    rel_gens = []
    for b in pre_base:
      b_form = LogarithmicDifferentialForm(n-1,b,self)
      w = dh.wedge(b_form)
      rel_gens.append(lift_to_basis([w],df_base))
    rel = vs_base.subspace(rel_gens)
    comp = orth_complement(vs_base,rel)
    #Lift
    rel_complex = []
    for vec in comp.basis():
      rel_complex.append(_weighted_sum(vec,df_base,self))
    return rel_complex
    

  def _differential_complement(self,form,*args):
    return [form.derivative()]

  def _differential_equivarient(self,form,*args):
    der = form.derivative()
    i = form.interior_product()
    if form.degree<1:
      return [der]
    else:
      return [der,i]

  def _differential_relative(self,form,*args):
    n = form.degree
    deg = self._p_graded_module(n).total_degree(form.vec)
    deg -= self.degree
    der = form.derivative()
    full = self._p_graded_module(n+1).homogeneous_part_basis(self.degree+deg)
    full_forms = [LogarithmicDifferentialForm(n+1,b,self) for b in full]
    full_space = VectorSpace(QQ,len(full))
    target_comp = self._complex_relative(n+1,deg)
    comp_vecs = []
    for b in target_comp:
      comp_vecs.append(lift_to_basis([b],full_forms))
    comp_vecs = full_space.subspace(comp_vecs)
    lift_full = lift_to_basis([der],full_forms)
    lift_prog = comp_vecs.basis_matrix()*vector(lift_full)
    return [_weighted_sum(lift_prog,target_comp,self)]
    
  
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

  def chain_complex(self,varient="complement",return_range=None,*args):
    cc = {}
    if return_range==None:
      return_range = (0,self.poly_ring.ngens())
    for i in range(return_range[0],return_range[1]+1):
      cc[i] = getattr(self,"_complex_"+varient)(i,*args)
    return cc

  def differential(self,form,varient="complement",*args):
    return getattr(self,"_differential_"+varient)(form,*args)

  def homology(self,varient="complement",*args):
    hom = {}
    cc = self.chain_complex(varient,(-1,self.poly_ring.ngens()+1),*args)
    #Some need wider return range
    for i in range(self.poly_ring.ngens()+1):
      vs_am = VectorSpace(QQ,len(cc[i]))
      if len(cc[i])==0:
        hom[i] = []
        continue
      if len(cc[i-1])!=0:
        d_im = []
        for b in cc[i-1]:
          d_b = self.differential(b,varient,*args)
          d_im.append(lift_to_basis(d_b,cc[i]))
        img = vs_am.subspace(d_im)
      else:
        img = vs_am.subspace([vs_am.zero()])
      if len(cc[i+1])!=0:
        d_ker = []
        for b in cc[i]:
          d_b = self.differential(b,varient,*args)
          d_ker.append(lift_to_basis(d_b,cc[i+1]))
        ker = (matrix(QQ,d_ker)).left_kernel()
      else:
        ker = vs_am
      quo = ker.quotient(img)
      hom[i] = []
      for b in quo.basis():
        vec = quo.lift(b)
        part_sum = LogarithmicDifferentialForm.make_zero(i,self)
        for c,f in zip(vec,cc[i]):
          part_sum = part_sum + c*f
        hom[i].append(part_sum)
    return hom
      
