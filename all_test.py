import unittest

from singular_module_test import TestSingularModule

from logarithmic_forms_test import TestHomogeneousWieghts
from logarithmic_forms_test import TestConvertSymToPoly
from logarithmic_forms_test import TestConvertPolyToSym

if __name__=="__main__":
  suite = unittest.TestLoader()
  
  suite.loadTestsFromTestCase(TestSingularModule)
  suite.loadTestsFromTestCase(TestHomogeneousWieghts)
  suite.loadTestsFromTestCase(TestConvertSymToPoly)
  suite.loadTestsFromTestCase(TestConvertPolyToSym)
  
  unittest.TextTestRunner(verbosity=2).run(suite)