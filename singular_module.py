import sage.all

from sage.structure.sage_object import SageObject

from sage.interfaces.singular import Singular

from sage.rings.integer import Integer
from sage.rings.rational import Rational

from sage.rings.ideal import Ideal

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
  mod = mod + ";\n "
  mod = mod + "groebner("+str(mod_name)+");\n"
  return mod

 
def _possible_gens(coeffs,ideal,var=0):
  coeff_gens = list(ideal.gens())
  poly_ring = coeffs[var].parent()
  if coeffs[var]==poly_ring.zero():
    return []
  for i,coeff in enumerate(coeffs):
    if not i==var:
      coeff_gens.append(coeff)
  poss_ideal = (Ideal(poly_ring,coeff_gens)).intersection(Ideal(poly_ring,coeffs[var]))
  poss_gens = []
  for g in poss_ideal.gens():
    poss_gens.append(g//coeffs[var])
  return poss_gens

def _generators_from_relation(rel_coeffs,ideal):
  poly_ring = rel_coeffs[0].parent()
  #posible first vaules
  poss_gens_0 = _possible_gens(rel_coeffs,ideal,0)
  #Particulars for thise values
  part_ideal = Ideal(poly_ring,rel_coeffs[1:]+list(ideal.gens()))
  poss_gens = []
  if not rel_coeffs[0].is_zero():
    for g in poss_gens_0:
      part_lift = (-g*rel_coeffs[0]).lift(part_ideal)
      part_lift = part_lift[:(len(rel_coeffs)-1)] # Drop ideal stuff
      poss_gens.append([g]+part_lift)
  else:
    poss_gens = [ [poly_ring.one()] + [poly_ring.zero() for _ in rel_coeffs[1:]] ]
  #Solve for when first value is zero
  if len(rel_coeffs)>1:
    red_coeffs = rel_coeffs[1:]
    red_gens = _generators_from_relation(red_coeffs,ideal)
    for red_gen in red_gens:
      poss_gens.append([poly_ring.zero()]+red_gen)
  return poss_gens

class SingularModule(SageObject):
  #TODO refactor to spawn only 1 singular process (per module?)

  def __init__(self,gens):
    self.gens = gens
    self.poly_ring = gens[0][0].parent()
    self.rank = len(gens[0])
    
  def create_ring_str(self):
    return "ring r=0,x(1.."+str(self.poly_ring.ngens())+"),dp;\n"
    
  def create_module_str(self,module_name="MD"):
    return _sing_mod(self.gens,module_name);
    
  def _contains_vector(self,vector):
    ring = ring = self.create_ring_str()
    module = self.create_module_str("modA")
    vector  = SingularModule([vector]).create_module_str("vec")
    vector.replace("module","vector")#Bad hack
    #See if reduction of vec is 0
    contain_code = "reduce(vec,std(modA));\n"
    all_code = ring+module+vector+contain_code;
    singular = Singular()
    output = singular.eval(all_code)
    last_line = output.split('\n')[-1];
    if last_line=="_[1]=0":
      return True
    return False
    
  def _contains_module(self,module):
    for g in module.gens:
      if not self.contains(g):
        return False
    return True
	
  def contains(self,object):
    if hasattr(object,"gens"):
      return self._contains_module(object)
    return self._contains_vector(object)
  
  def intersection(self,module):
    ring = self.create_ring_str()
    moduleA = self.create_module_str("modA")
    moduleB = module.create_module_str("modB")
    intersect_code = "module inter=intersect(modA,modB);\n matrix mat_inter = inter;\n mat_inter;\n";
    all_code = ring+moduleA+moduleB+intersect_code;
    singular = Singular()
    output = singular.eval(all_code)
    mod_iter = SingularModule.create_from_singular_matrix(self.poly_ring,output,"mat_inter");
    return mod_iter
    
  def equals(self,module):
    if self.contains(module):
      if module.contains(self):
        return True
    return False
    
  def reduce_generators(self):
    for g in self.gens:
      other_gens = [ h for h in self.gens if not h==g ]
      other_mod = SingularModule(other_gens)
      if other_mod.contains(g):
        self.gens.remove(g)
        
  def standard_basis(self):
    ring = self.create_ring_str()
    moduleA = self.create_module_str("modA")
    std_code = "module std_mod=std(modA);\n matrix mat_std = std_mod;\n mat_std;\n"
    all_code = ring+moduleA+std_code
    singular = Singular()
    output = singular.eval(all_code)
    std_mod = SingularModule.create_from_singular_matrix(self.poly_ring,output,"mat_std");
    return std_mod.gens
    
  def ambient_free_module(self):
    return SingularModule.create_free_module(self.rank,self.poly_ring)
    
  def is_free(self):
    if len(self.gens)==1:
      return True
    std_gens = self.standard_basis()
    for i,g in enumerate(std_gens):
      other_gens = std_gens[:i]+std_gens[i+1:]
      g_mod = SingularModule([g])
      other_mod = SingularModule(other_gens)
      inter_mod = other_mod.intersection(g_mod)
      zero = SingularModule.create_zero_module(1,self.poly_ring)
      if not inter_mod.equals(zero):
        return False
    return True
    
  def __repr__(self):
    return "A submodule of the rank "+str(self.rank)+" free module over "+repr(self.poly_ring)+" with generators "+repr(self.gens)

  def lift(self,vector,make_linear=False):
    ring = self.create_ring_str()
    moduleA = self.create_module_str("modA")
    vec_mod = SingularModule([vector])
    vec_module = vec_mod.create_module_str("modV")
    std_code = "module lift_mod=lift(modA,modV);\n matrix mat_lift = lift_mod;\n mat_lift;\n"
    all_code = ring+moduleA+vec_module+std_code
    singular = Singular()
    output = singular.eval(all_code)
    lift_mod = SingularModule.create_from_singular_matrix(self.poly_ring,output,"mat_lift");
    lift = lift_mod.gens[0]
    if make_linear:
      for i,l in enumerate(lift):
        lift[i] = l.constant_coefficient()
    return lift
    
  @classmethod
  def create_free_module(cls,n,poly_ring):
    gens = []
    for i in range(n):
      gen = [ poly_ring.zero() for _ in range(n) ]
      gen[i]=poly_ring.one()
      gens.append(gen)
    return SingularModule(gens);
    
  @classmethod
  def create_zero_module(cls,n,poly_ring):
    zero = poly_ring.zero()
    zero_gen = [zero for _ in range(n)]
    return SingularModule([zero_gen])
    
  @classmethod
  def create_from_relation(cls,relation,ideal):
    gens = _generators_from_relation(relation,ideal)
    return SingularModule(gens)
    
  @classmethod
  def create_from_singular_matrix(cls,poly_ring,singular_output,matrix_name="MM"):
    mat = {}
    for line in singular_output.split('\n'):
      if line.find(matrix_name)==0:
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
                try:
                  num = part.split('(')[1]
                  num = num.split(')')[0]
                  num = Integer(num)
                  try:
                    pow = Integer(part.split('^')[1])
                  except:
                    pow = 1
                  coeff = coeff * poly_ring.gens()[num-1]**pow;
                except:
                  pass
            if mon!='':
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
    return SingularModule(gens)
    
  @classmethod
  def create_from_relations(cls,relations,ideals):
    rank = len(relations[0])
    poly_ring = relations[0][0].parent()
    module = SingularModule.create_free_module(rank,poly_ring)
    for rel,ideal in zip(relations,ideals):
      rel_mod = SingularModule.create_from_relation(rel,ideal)
      module = module.intersection(rel_mod)
    return module
