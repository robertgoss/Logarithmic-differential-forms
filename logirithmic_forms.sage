#!/home/robert/sage-4.7.alpha2/sage

def _relations_generators(relations,ideal):
  gens = _relation_generators(relations[0],ideal)
  for rel in rels[1:]:
    new_rel = _apply_relation(gens,rel)
    new_gens = _relation_generators(new_rel,ideal)
    gens = _unapply_relation(gens,new_gens)
  return gens
  
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
  
#checked
def _apply_relation(generators,relation):
  n = len(relation)
  m = len(generators)
  poly_ring = relation[0].parent()
  coeffs = []
  for i in range(m):
    coeff = poly_ring.zero()
    for j in range(n):
      coeff += relation[j]*generators[i][j]
    coeffs.append(coeff)
  return coeffs
  
def _unapply_relation(old_gens,new_gens):
  gens = []
  poly_ring = old_gens[0][0].parent()
  for t in new_gens:
    gen = []
    for i in range(len(old_gens[0])):
      coeff = poly_ring.zero()
      for j in range(len(old_gens)):
        coeff += old_gens[j][i] * t[j]
      gen.append(coeff)
    gens.append(gen)
  return gens
  

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
  
  print "Testing relations generator"
  rels = [[x^2,x*y,z],[x,z^2,x+y]]
  print "Rels: ",rels
  print "Ideal: ",ideal
  print "Gens: ",_relations_generators(rels,ideal)