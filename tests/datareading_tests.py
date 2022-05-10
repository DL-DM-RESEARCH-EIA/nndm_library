import unittest
from nndm_library import ReadLhe
from nndm_library import Constants

class ReadFileTestCase(unittest.TestCase):
    def setUp(self) -> None:
        name = "tests/data_tests/signal/eta_decay_events_mk_0.014980679431428716_eps2_3.883209914996183e-10.lhe"
        self.read_file = ReadLhe(name)

    def test_extract_params_from_path(self):
        expected_params = {"particle_type" : Constants.ETA_ID, "mk": 0.014980679431428716, "eps2" : 3.883209914996183e-10}
        res_parameters = self.read_file.extract_params_from_path()
        for param_name in res_parameters.keys():
            self.assertAlmostEqual(res_parameters[param_name], expected_params[param_name])
