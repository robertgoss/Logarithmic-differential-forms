import sage.all

from sage.rings.rational_field import QQ
from sage.modules.free_module import VectorSpace
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.symbolic.ring import var

from random import randrange
from random import random

def crossing_divisor(n=3,var="z"):
  poly_ring = PolynomialRing(QQ,n,var)
  div = poly_ring.one()
  for g in poly_ring.gens():
    div *= g
  return div

def whitney_divisor(var=None):
  if var:
    poly_ring = PolynomialRing(QQ,3,var)
  else:
    poly_ring = PolynomialRing(QQ,3,"xyz")
  gens = poly_ring.gens()
  return gens[0]**2*gens[1]-gens[2]**2

def breiskorn_pham_divisor(n,var="z"):
  poly_ring = PolynomialRing(QQ,n,var)
  for i,g in enumerate(poly_ring.gens()):
    if i==0:
      div = g**3
    else:
      div += g**2
  return div

def rand_pham_divisor(n,var="z"):
  poly_ring = PolynomialRing(QQ,n,var)
  div = poly_ring.zero()
  for g in poly_ring.gens():
    div += g**(randrange(2,9))
  return div

def braid_divisor(n,var="z"):
  poly_ring = PolynomialRing(QQ,n,var)
  div = poly_ring.one()
  gens = poly_ring.gens()
  for i in range(n):
    for j in range(i+1,n):
      div *= (gens[i]-gens[j])
  return div

def rand_w_hom_divisor(n,degs=None,mon_num=None,var="z"):
  if degs==None:
    degs = [randrange(2,6) for _ in range(n)]
  deg = sum(degs)
  if mon_num==None:
    mon_num = randrange(2,8)
  poly_ring = PolynomialRing(QQ,n,var)
  div = poly_ring.zero()
  min_w = min(degs)
  for i in mon_num:
    expo = [0]*n
    cur_deg = 0
    while cur_deg!=deg:
      if cur_deg>deg:
        expo = [0]*n
        cur_deg = 0
      if deg-cur_deg<min_w:
        expo = [0]*n
        cur_deg = 0
      next_g = randrange(0,n)
      expo[next_g] += 1
      cur_deg += degs[next_g]
    coeff = randrange(-n,n)/n
    mon = poly_ring.one()
    for i,e in enumerate(expo):
      mon *= poly_ring.gens()[i]**e
    div += coeff*mon
  return div

def rand_hom_divisor(n,deg=None,mon_num=None,var="z"):
  if deg==None:
    deg = randrange(2,6)
  return rand_w_hom_divisor(n,[1]*deg,mon_num,var)
