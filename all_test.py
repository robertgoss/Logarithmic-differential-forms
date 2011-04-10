#!/usr/local/bin/sage

import unittest

from singular_module_test import TestSingularModule

from logarithmic_forms_test import TestHomogeneousWieghts

from logarithmic_form_test import TestConvertSymToPoly
from logarithmic_form_test import TestConvertPolyToSym
from logarithmic_form_test import TestLogarithmicDifferentialForm

from logarithmic_forms_test import TestSkewIter
from logarithmic_forms_test import TestLogartihmicDifferentialForms

from logarithmic_derivations_test import TestLogaritmicDerivations

from graded_module_test import TestWeightedMaxDegree
from graded_module_test import TestWeightedMinDegree
from graded_module_test import TestGradedModule
from graded_module_test import TestMonomialsOfOrder

cases = [ TestSingularModule , TestHomogeneousWieghts , TestConvertSymToPoly,
          TestConvertPolyToSym , TestSkewIter, TestLogartihmicDifferentialForms,
          TestLogaritmicDerivations, TestWeightedMinDegree, TestWeightedMaxDegree,
          TestMonomialsOfOrder, TestGradedModule, TestLogarithmicDifferentialForm]

if __name__=="__main__":
  suites = []
  for case in cases:
    suite = unittest.TestLoader().loadTestsFromTestCase(case)
    suites.append(suite)
  suite = unittest.TestSuite(suites)
  unittest.TextTestRunner(verbosity=2).run(suite)
