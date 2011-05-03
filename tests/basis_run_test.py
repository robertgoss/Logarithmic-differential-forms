#!/usr/local/bin/sage

import nose

from logarithmic_forms import LogarithmicDifferentialForms

from examples import *

def _crossing_divisor(n):
  div = examples.crossing_divisor(n)
  logdf = LogarithmicDifferentialForms(div)
  logdf.homology("complement")
  
def test_crossing_divisor():
  for i in range(2,8):
    yield _crossing_divisor, i

def _braid_divisor(n):
  div = examples.braid_divisor(n)
  logdf = LogarithmicDifferentialForms(div)
  logdf.homology("complement")
  
def test_braid_divisor():
  for i in range(2,8):
    yield _braid_divisor, i

def _breiskorn_pham_divisor(n):
  div = examples.breiskorn_pham_divisor(n)
  logdf = LogarithmicDifferentialForms(div)
  logdf.homology("complement")
  
def test_breiskorn_pham_divisor():
  for i in range(2,8):
    yield _breiskorn_pham_divisor, i

def _rand_pham_divisor(n):
  div = examples.rand_pham_divisor(n)
  logdf = LogarithmicDifferentialForms(div)
  logdf.homology("complement")
  
def test_rand_pham_divisor():
  for i in range(2,8):
    for _ in range(20):
      yield _rand_pham_divisor, i

def _rand_w_hom_divisor():
  div = examples.rand_w_hom_divisor(n)
  logdf = LogarithmicDifferentialForms(div)
  logdf.homology("complement")

def test_rand_w_hom_divisor():
  for i in range(2,8):
    for _ in range(20):
      yield _rand_w_hom_divisor, i

if __name__=="__main__":
  nose.main()
