
import sage.all

from sage.matrix.constructor import matrix

from sage.rings.rational_field import QQ

from singular_module import SingularModule

def _wieghted_degree(mon,var_wieghts):
  sum = 0
  for w,ex in zip(var_wieghts,mon.exponents()[0]):
    sum += w*ex
  return sum

def wieghted_max_degree(poly,var_wieghts):
  max_wieght = 0;
  for ex in poly.exponents():
    wieght = sum([e*w for e,w in zip(ex,var_wieghts)])
    max_wieght = max(max_wieght,wieght)
  return max_wieght

def wieghted_min_degree(poly,var_wieghts):
  min_wieght = poly.degree()*max(var_wieghts);
  for ex in poly.exponents():
    wieght = sum([e*w for e,w in zip(ex,var_wieghts)])
    min_wieght = min(min_wieght,wieght)
  return min_wieght
  
def _inter_monomials_upto_order(k,poly_ring,var_wieghts,depth):
  if depth>=poly_ring.ngens():
    yield poly_ring.one()
  else:
    for i in range(k/var_wieghts[depth]+1):
      for mon in  _inter_monomials_upto_order(k,poly_ring,var_wieghts,depth+1):
        new_mon = poly_ring.gens()[depth]**i * mon
        if _wieghted_degree(new_mon,var_wieghts) <= k:
          yield new_mon
          
def monomials_upto_order(k,poly_ring,var_wieghts):
  #Catch wieghts giving infinite recursion?
  mons = []
  for mon in _inter_monomials_upto_order(k,poly_ring,var_wieghts,0):
    mons.append(mon)
  return mons
  
def monomials_between_orders(n,k,poly_ring,var_wieghts):
  #Catch wieghts giving infinite recursion?
  mons = []
  for mon in _inter_monomials_upto_order(k,poly_ring,var_wieghts,0):
    if _wieghted_degree(mon,var_wieghts) >= n:
      mons.append(mon)
  return mons


class GradedModule(SingularModule):

  def __init__(self,gens,column_wieghts,var_wieghts):
    SingularModule.__init__(gens)
    self.column_wieghts = column_wieghts
    self.var_wieghts = var_wieghts
    
  def _max_degree(self,vector):
    max_deg = 0
    for p,w in zip(vector,self.column_wieghts):
      max_deg = max(max_deg,w+wieghted_max_degree(p,self.var_wieghts))
    return max_deg
    
     
  def _min_degree(self,vector):
    min_deg = self._max_degree(vector)
    for p,w in zip(vector,self.column_wieghts):
      min_deg = max(max_deg,w+wieghted_max_degree(p,self.var_wieghts))
    return max_deg
    
  def homogeneous_part_basis(self,k):
