#!/usr/local/bin/sage

from logarithmic_forms import LogarithmicDifferentialForms

import examples


def _crossing_divisor(n):
  print "Crossing divisor: ",n
  div = examples.crossing_divisor(n)
  logdf = LogarithmicDifferentialForms(div)
  logdf.homology("complement")
  
def test_crossing_divisor():
  for i in range(2,7):
    yield _crossing_divisor, i

def _braid_divisor(n):
  print "Braid Divisor: ",n
  div = examples.braid_divisor(n)
  logdf = LogarithmicDifferentialForms(div)
  logdf.homology("complement")
  
def test_braid_divisor():
  for i in range(2,5):
    yield _braid_divisor, i

def _breiskorn_pham_divisor(n):
  print "Breiskorn Pham Divisor: ",n
  div = examples.breiskorn_pham_divisor(n)
  logdf = LogarithmicDifferentialForms(div)
  logdf.homology("complement")
  
def test_breiskorn_pham_divisor():
  for i in range(2,6):
    yield _breiskorn_pham_divisor, i
