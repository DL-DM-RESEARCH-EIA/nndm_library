import unittest
from unittest import result

from nndm_library import ReadFile

class ReadFileTestCase(unittest.TestCase):

    def setUp(self) -> None:
        name = "data/eta_decay_events_mk_0.014980679431428716_eps2_3.883209914996183e-10.lhe"
        self.read_file = ReadFile(name)

    def test_extract_params_from_name(self):
        type_particle, mk, eps2 = self.read_file.extract_params_from_name()
        self.assertAlmostEqual(type_particle, 0)  # testing for eta particle
        self.assertAlmostEqual(mk, 0.014980679431428716)  # testing for eta particle
        self.assertAlmostEqual(eps2, 3.883209914996183e-10)  # testing for eta particle
        
