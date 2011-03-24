#!/home/robert/sage-4.7.alpha2/sage


def _compute_gens(relations,divisor):
  #Deal with 1 var base case
  if(len(relations[0])==1):
    return _posible_gens(relations,divisor)
  #Compute posible solutions to the first var
  # and other solutuions based on this
  r_1 = [rel for rel in relations if rel[0]!=0]
  P_1 = _posible_gens(r_1,divisor)
  part = []
  for gen in P_1:
    part.append(_particular_solution(gen,relations,divisor))
  #Compute the general solutions where the first var is zero
  r_gen = copy(r_1)
  for r in r_gen:
    r[0]=0
  general_sol = _compute_gens(r_gen,divisor)
  first_gens = _define_gens(part,general_sol)
  #Get the gen set to satify remaining relations
  r_2 = [rel for rel in relations if rel[0]==0]
  r_sub = _subs_gens(first_gens,r_2)
  subs_gen = _compute_gens(r_sub,divisor)
  #Un sub
  return _unsub_gens(subs_gens,first_gen)
  
def _posible_gens(relations,divisor,var):
  posible = (1)*divisor.parent()
  for rel in relations:
    if not var in rel.keys():
      continue;
    other_gens = [coeff for v,coeff in rel.iteritems() if v!=var]
    other_gens.append(divisor)
    I = other_gens*divisor.parent()
    I = I.intersection((rel[var],)*divisor.parent())
    new_pos = []
    for g in I.gens():
      new_pos.append(g/rel[var])
    posible = posible.intersection(new_pos*divisor.parent())
  return posible
  
def _particular_solution(val,relations,divisor,var):
  pass;
  
def _define_gens(particular,general):
  pass

def _subs_gens(gens,relations):
  pass;

def _unsub_gens(subs_gen,initial_gen):
  pass;

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
  return rels;
  
if __name__=="__main__":  
  C.<x,y,z> = PolynomialRing(QQ,3);
  h = x^2*y - z^2
  print  "0: ",_compute_p_relations(h,0)
  print  "1: ",_compute_p_relations(h,1)
  print  "2: ",_compute_p_relations(h,2)
  
  rel = _compute_p_relations(h,0)
  print "p: poss: ",_posible_gens(rel,h,tuple([]))
  rel = _compute_p_relations(h,1)
  print "p_x: poss: ",_posible_gens(rel,h,(0,))
  print "p_y: poss: ",_posible_gens(rel,h,(1,))
  print "p_z: poss: ",_posible_gens(rel,h,(2,))
  rel = _compute_p_relations(h,2)
  print "p_xy: poss: ",_posible_gens(rel,h,(0,1))
  print "p_xz: poss: ",_posible_gens(rel,h,(0,2))
  print "p_yz: poss: ",_posible_gens(rel,h,(1,2))

