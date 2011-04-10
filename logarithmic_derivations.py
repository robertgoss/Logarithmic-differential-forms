
import sage.all

from singular_module import SingularModule

from sage.rings.ideal import Ideal

class LogarithmicDerivations(SingularModule):
  
  def __init__(self,divisor):
    self.divisor = divisor
    poly_ring = divisor.parent()
    #Compute gen set
    rel = [poly_ring.zero() for _ in range(poly_ring.ngens())]
    for c,mon in divisor:
      exp = mon.exponents()[0]
      for i in range(poly_ring.ngens()):
        rel[i] = rel[i] + exp[i]*mon//poly_ring.gens()[i]
    ideal = Ideal(poly_ring,[divisor])
    base_mod = SingularModule.create_from_relation(rel,ideal)
    SingularModule.__init__(self,base_mod.gens)
    self.relation = rel
