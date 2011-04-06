#!/usr/local/bin/sage

import unittest

from singular_module_test import TestSingularModule

from logarithmic_forms_test import TestHomogeneousWieghts
from logarithmic_forms_test import TestConvertSymToPoly
from logarithmic_forms_test import TestConvertPolyToSym
from logarithmic_forms_test import TestSkewIter
from logarithmic_forms_test import TestLogartihmicDifferentialForms

from logarithmic_derivations_test import TestLogaritmicDerivations

cases = [ TestSingularModule , TestHomogeneousWieghts , TestConvertSymToPoly,
          TestConvertPolyToSym , TestSkewIter, TestLogartihmicDifferentialForms,
          TestLogaritmicDerivations]

if __name__=="__main__":
  suites = []
  for case in cases:
    suite = unittest.TestLoader().loadTestsFromTestCase(case)
    suites.append(suite)
  suite = unittest.TestSuite(suites)
  unittest.TextTestRunner(verbosity=2).run(suite)