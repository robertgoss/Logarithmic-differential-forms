
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

  def latex_basis(self):
    string = "\\begin{enumerate}\n"
    for g in self.gens:
      string = string + "\item "
      for g_i,g_p in enumerate(g):
        parts = []
        if not g_p.is_zero():
          parts.append("$"+g_p.__repr__())
          parts[-1] = parts[-1] + "\\frac{\\partial}{\\partial "+self.poly_ring.gens()[g_i].__repr__() + "}$\n"
        string = string + "+".join(parts)
    string = string + "\\end{enumerate}"
    return string
