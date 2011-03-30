#!/home/robert/sage-4.7.alpha2/sage

#Todo add,test intersection
#Todo add in
#Todo add reduce
  
#Checked
def _possible_gens(coeffs,ideal,var=0):
  coeff_gens = list(ideal.gens())
  poly_ring = coeffs[var].parent()
  if coeffs[var]==poly_ring.zero():
    return []
  for i,coeff in enumerate(coeffs):
    if not i==var:
      coeff_gens.append(coeff)
  poss_ideal = (coeff_gens*poly_ring).intersection((coeffs[var])*poly_ring)
  poss_gens = []
  for g in poss_ideal.gens():
    poss_gens.append(g//coeffs[var])
  return poss_gens
  
#Checked 
def _relation_generators(coeffs,ideal):
  poly_ring = coeffs[0].parent()
  #posible first vaules
  poss_gens_0 = _possible_gens(coeffs,ideal,0)
  #Particulars for thise values
  part_ideal = (coeffs[1:]+list(ideal.gens()))*poly_ring
  poss_gens = []
  for g in poss_gens_0:
    part_lift = (g*coeffs[0]).lift(part_ideal)
    part_lift = part_lift[:(len(coeffs)-1)] # Drop ideal stuff
    poss_gens.append([g]+part_lift)
  #Solve for when first value is zero
  if len(coeffs)>1:
    red_coeffs = coeffs[1:]
    red_gens = _relation_generators(red_coeffs,ideal)
    for red_gen in red_gens:
      poss_gens.append([poly_ring.zero()]+red_gen)
  return poss_gens
  
def _sing_ring(ring,var_name="x"):
  return "ring r = 0,(x(1.."+str(n_vars)+")),dp;\n"

  
#Checked
def _make_sing_mod_from_gens(gens,mod_name="MD"):
  poly_ring = gens[0][0].parent()
  n_vars = poly_ring.ngens();
  #create ring
  mod = _sing_ring(poly_ring)
  mod = mod + "module "+str(mod_name)+" = "
  first_gen = true
  for gen in gens:
    if first_gen:
      first_gen = false
    else:
      mod = mod + ","
    mod += "["
    first_poly = true
    for poly in gen:
      if first_poly:
        first_poly = false;
      else:
        mod = mod + ","
      if poly==poly_ring.zero():
        mod = mod + "0"
      first_mon = true
      for mon in poly:
        if first_mon:
          first_mon = false;
        else:
          mod = mod + "+"
        mod = mod + str(mon[0])
        ex = mon[1].exponents()[0]
        for ex_i,ex_m in enumerate(ex):
          if ex_m != 0:
            mod = mod + "*x("+str(ex_i+1)+")^"+str(ex_m);
    mod = mod + "]"
  mod = mod + ";\n"
  return mod
  
#Checked
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
      for mon in eqn.split('+'):
        coeff = poly_ring.one()
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
            coeff = coeff * poly_ring.gens()[num-1]^pow;
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
  
#checked ish
def _simplify_gens_by_div(gens,div):
  new_gens = []
  for gen in gens:
    new_gen = []
    for c in gen:
      new_gen.append(c%div)
    new_gens.append(new_gen)
  return new_gens
  
class SingularModule():
  def __init__(self,gens):
    self.gens = gens
    
  def intersection(self,module):
    mod = _make_sing_mod_from_gens(gens,"self_mod")
    other_mod = _make_sing_mod_from_gens(module.gens(),"other_mod")
    sing = mod + ";\n" + other_mod + ";\n"
    sing  = sing + "int_mod = intersection(self_mod,other_mod); matrix int_mat = int_mod; int_mat"
    poly_ring = gens[0][0].parent()
    gens = _gens_from_sing_mod(poly_ring,singular.eval(sing),"int_mat")
    return SingularModule(gens);
    
  def create_ring_str(self):
    return _sing_ring
    
  def reduce(self,div):
    pass;
    
  def contains(self,g):
    pass;
    
  @classmethod
  def from_relations(self,rels,ideal):
    mod = SingularModule(_relation_generators(rels[0],ideal))
    for rel in rels[1:]:
      other_mod = SingularModule(_relation_generators(rel,ideal))
      mod = mod.intersection(other_mod)
    return mod;
    
  @classmethod
  def from_relation(self,rel,ideal):
    return SingularModule(_relation_generators(rel,ideal));
  
#Checked
def _compute_p_relations(divisor,p):
  diffs = []
  for g in divisor.parent().gens():
    diffs.append(divisor.derivative(g))
  rels = []
  for s in Set(range(divisor.parent().ngens())).subsets(p+1):
    rel = {}
    for i in range(divisor.parent().ngens()):
      if i in s:
        sign = sum([1 for j in s if j<i])%2
        if sign==1:
          rel[tuple(s.difference([i]))] = diffs[i]
        else:
          rel[tuple(s.difference([i]))] = - diffs[i]
    rels.append(rel)
  #reform relations
  new_rels = []
  for rel in rels:
    new_rel = []
    for s in Set(range(divisor.parent().ngens())).subsets(p):
      t = tuple(s)
      if t in rel.keys():
        new_rel.append(rel[t])
      else:
        new_rel.append(divisor.parent().zero())
    new_rels.append(new_rel)
  return new_rels;
  
if __name__=="__main__":  
  C.<x,y,z> = PolynomialRing(QQ,3);
  h = x^2*y - z^2
  print "Testing Generators"
  print  "0: ",_compute_p_relations(h,0)
  print  "1: ",_compute_p_relations(h,1)
  print  "2: ",_compute_p_relations(h,2)
  
  print "Testing Possible generators"
  rel = [x^3,y+x,0,z]
  ideal = (x^2*y-z^2)*C
  print "Rel: ",rel
  print "Ideal: ",ideal
  print "Pos 0:",_possible_gens(rel,ideal,0)
  print "Pos 1:",_possible_gens(rel,ideal,1)
  print "Pos 2:",_possible_gens(rel,ideal,2)
  print "Pos 3:",_possible_gens(rel,ideal,3)
  
  print "Testing Relation generator"
  rel = [x^3,y+x,C.zero(),z]
  ideal = (x^2*y-z^2)*C
  print "Rel: ",rel
  print "Ideal: ",ideal
  print "Gens: ",_relation_generators(rel,ideal)
  
  
  print "Testing mod string generator"
  gens = [[x^4,C.zero(),x+y+z^3],[4*x*y*z,C.one(),y*x]]
  print "Gens: ",gens
  string = _make_sing_mod_from_gens(gens)
  print singular.eval(string+";\n MD ;\n")
  
  print "Testing geting gens from singular modular"
  print "Gens: ",gens
  string = _make_sing_mod_from_gens(gens)
  string = singular.eval(string+";\n matrix MM = MD;\n MM;")
  gens = _gens_from_sing_mod(C,string,"MM")
  print gens

  