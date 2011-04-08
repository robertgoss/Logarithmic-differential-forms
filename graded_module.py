
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
  if poly==poly.parent().zero():
    return 0
  max_wieght = 0;
  for ex in poly.exponents():
    wieght = sum([e*w for e,w in zip(ex,var_wieghts)])
    max_wieght = max(max_wieght,wieght)
  return max_wieght

def wieghted_min_degree(poly,var_wieghts):
  if poly==poly.parent().zero():
    return 0;
  min_wieght = poly.degree()*max(var_wieghts);
  for ex in poly.exponents():
    wieght = sum([e*w for e,w in zip(ex,var_wieghts)])
    min_wieght = min(min_wieght,wieght)
  return min_wieght

def monomials_of_order(k,poly_ring,var_wieghts,start=0):
  if k==0:
    yield poly_ring.one()
    return
  if start < len(var_wieghts)-1:
    for i in range(k/var_wieghts[start]+1):
      for mon in monomials_of_order(k-(i*var_wieghts[start]),poly_ring,var_wieghts,start+1):
        new_mon = poly_ring.gens()[start]**i * mon
        if wieghted_max_degree(new_mon,var_wieghts)==k:
          yield new_mon
  else:
    new_mon  = poly_ring.gens()[start]**(k/var_wieghts[start])
    if wieghted_max_degree(new_mon,var_wieghts)==k:
      yield new_mon

class GradedModule(SingularModule):

  def __init__(self,gens,column_wieghts,var_wieghts):
    SingularModule.__init__(self,gens)
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

  def monomial_basis(self,k):
    basis = []
    zero = [self.poly_ring.zero() for _ in range(self.rank)]
    for i in range(self.rank):
      for mon in monomials_of_order(k-self.column_wieghts[i],self.poly_ring,self.var_wieghts):
        mon_vec = list(zero)
        mon_vec[i] = mon
        basis.append(tuple(mon_vec))
    return basis

  def get_homogeneous_parts(self,vector):
    parts = {}
    for i,shift in enumerate(self.column_wieghts):
      poly = vector[i]
      for c,mon in poly:
        deg = _wieghted_degree(mon,self.var_wieghts)+shift
        try:
          parts[deg][i] = parts[deg][i] + c*mon
        except:
          parts[deg] = [self.poly_ring.zero() for _ in range(self.rank)]
          parts[deg][i] = c*mon
    return parts
    
  def homogeneous_part_basis(self,k):
    pass
