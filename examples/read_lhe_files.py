import unittest
from nndm_library import ReadLhe
from nndm_library import Constants
import numpy as np

class LheReadTestCase(unittest.TestCase):
    def main():        
        ## Testing reading of lhe files

        # Single lhe file
        filename = "tests/data_tests/signal/eta_decay_events_mk_0.014980679431428716_eps2_3.8832099149961855e-09.lhe"
        print("Reading %s filetring by particle_ids, var_of_interest, and outgoing particles" % filename)
        print("*" * 50)
        file = ReadLhe(path=filename, particle_ids=[Constants.ELECTRON_ID], var_of_interest=['px'], outgoing=True)
        print(file.data)

        print("Reading %s filtering by particle_ids, and outgoing particles (all available variables are read" % filename)
        print("*" * 50)
        file = ReadLhe(path=filename, particle_ids=[Constants.ELECTRON_ID], outgoing=True)
        print(file.data)

        print("Reading %s filtering by particle_ids (all available variables for particles going in and out are read" % filename)
        print("*" * 50)
        file = ReadLhe(path=filename, particle_ids=[Constants.ELECTRON_ID])
        print(file.data)
        print("in out out states")
        print(np.unique(file.data["status"]))

        print("Reading from %s all particles for all available variables for particles going in and out" % filename)
        print("*" * 50)
        file = ReadLhe(path=filename)
        print(file.data)
        print("particle ids")
        print(np.unique(file.data["id"]))

        # files in a directory
        print("\n\nRead all files like the one before, but in a given directory")
        filename = "tests/data_tests/signal/"
        print("Reading directory %s" % filename)
        print("*" * 50)
        file = ReadLhe(path=filename)
        print(file.data)
        print(file.files_dir)
        res_dict = file.extract_params_from_path()
        print(res_dict)

        # files in a directory recursively
        print("\n\nRead all basic files like the one before, but all found recursively in a given directory")
        filename = "tests/data_tests/signal/"
        print("Reading all files found recursively in directory %s" % filename)
        print("*" * 50)
        file = ReadLhe(filename, recursive=True, relabel_events=True)
        print(file.data)
        print(file.files_dir)

        res_dict = file.extract_params_from_path()
        print(res_dict)


LheReadTestCase.main()
