from sage.interfaces.singular import Singular

from sage.rings.integer import Integer
from sage.rings.rational import Rational

class NotImplementedException(Exception):
  pass

def _sing_ring(ring,var_name="x"):
  return "ring r = 0,(x(1.."+str(n_vars)+")),dp;\n"
  
def _sing_mod(gens,mod_name="MD",vars_name="x"):
  poly_ring = gens[0][0].parent()
  mod = "module "+str(mod_name)+"="
  first_gen = True
  for gen in gens:
    if first_gen:
      first_gen = False
    else:
      mod = mod + ","
    mod += "["
    first_poly = True
    for poly in gen:
      if first_poly:
        first_poly = False;
      else:
        mod = mod + ","
      if poly==poly_ring.zero():
        mod = mod + "0"
      first_mon = True
      for mon in poly:
        if first_mon:
          first_mon = False;
        else:
          mod = mod + "+"
        mod = mod + str(mon[0])
        ex = mon[1].exponents()[0]
        for ex_i,ex_m in enumerate(ex):
          if ex_m != 0:
            mod = mod + "*"+str(vars_name)+"("+str(ex_i+1)+")^"+str(ex_m);
    mod = mod + "]"
  mod = mod + ";\n"
  return mod

def _gens_from_sing_mod(poly_ring,mod_string,mat_name='MM'):
  mat = {}
  for line in mod_string.split('\n'):
    if line.find(mat_name)==0:
      indices = line.split('[')[1]
      indices = indices.split(']')[0]
      indices = indices.split(',')
      indices = [ int(i) for i in indices ]
      eqn = line.split('=')[1]
      eq = poly_ring.zero()
      for mon_pos in eqn.split('+'):
        for i,mon in enumerate(mon_pos.split('-')):
          if i==0:
            coeff = poly_ring.one()
          else:
            coeff = -poly_ring.one()
          for part in mon.split('*'):
            try:
              coeff = coeff * Rational(part)
            except:
              num = part.split('(')[1]
              num = num.split(')')[0]
              num = Integer(num)
              try:
                pow = Integer(part.split('^')[1])
              except:
                pow = 1
              coeff = coeff * poly_ring.gens()[num-1]**pow;
          eq = eq + coeff;
      mat[(indices[0],indices[1])] = eq
  ngens = max([ k[1] for k in mat.iterkeys() ])
  nel = max([ k[0] for k in mat.iterkeys() ])
  gens = []
  for i in range(ngens):
    gen = []
    for j in range(nel):
      gen.append(mat[(j+1,i+1)])
    gens.append(gen)
  return gens;
 

class SingularModule(SageObject):

  def __init__(self,gens):
    self.gens = gens
    self.poly_ring = gens[0][0].parent()
    
  def create_ring_str(self):
    return "ring r=0,x(1.."+str(self.poly_ring.ngens())+"),dp;\n"
    
  def create_module_str(self,module_name="MD"):
    return _sing_mod(self.gens,module_name);
	
  def contains(self,vector):
    #based on our code we need to deal with zero seperatly
    is_zero = True
    for coord in vector:
      if not coord.is_zero():
        is_zero = False
    if is_zero:
      return True;
    ring = self.create_ring_str()
    moduleA = self.create_module_str("modA")
    vec_module = SingularModule([vector])
    moduleB = vec_module.create_module_str("modB")
    contain_code = "module cont = intersect(modA,modB);\n cont;\n"
    all_code = ring+moduleA+moduleB+contain_code;
    singular = Singular()
    output = singular.eval(all_code)
    lines = output.split('\n')
    last_line = lines[-1]
    if last_line == "cont[1]=0":
      return False
    return True;
  
  def intersection(self,module):
    ring = self.create_ring_str()
    moduleA = self.create_module_str("modA")
    moduleB = module.create_module_str("modB")
    intersect_code = "module inter=intersect(modA,modB);\n matrix mat_inter = inter;\n mat_inter;\n";
    all_code = ring+moduleA+moduleB+intersect_code;
    singular = Singular()
    output = singular.eval(all_code)
    gens = _gens_from_sing_mod(self.poly_ring,output,"mat_inter");
    return SingularModule(gens)
    
  def equals(self,module):
    equal = True
    for gen in self.gens:
      if not module.contains(gen):
        equal = False
        break
    if equal:
      for gen in module.gens:
        if not self.contains(gen):
          equal = False
          break
      if equal:
        return True
    return False
    
  @classmethod
  def create_from_relation(cls,relation):
    raise NotImplementedException()
    
  @classmethod
  def create_from_relations(cls,relations):
    raise NotImplementedException()