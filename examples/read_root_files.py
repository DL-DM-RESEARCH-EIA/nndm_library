import unittest
from nndm_library import ReadRoot

class RootReadTestCase(unittest.TestCase):
    def main():        
        ## Testing reading of root files
        print("\n" * 2)
        print("*" * 50)
        print("Read Root file")
        file = ReadRoot("tests/data_tests/background/v_e_scattering/onantinuelepton10125.root")
        print(file.data)

        print("\n" * 2)
        print("*" * 50)
        print("Read background from all the roots in a given directory")
        file = ReadRoot("tests/data_tests/background/v_e_scattering/", relabel_events=True)
        print(file.data)
        print(file.files_dir)

        last_index = file.data.index.tolist()[-1]
        if type(last_index) != int:
            last_index = last_index[0]

        print('number of index times two: %d, lenght of df: %d' % (last_index * 2, len(file.data)) )

        print("\n" * 2)
        print("*" * 50)
        print("Read background all the root found walking through all the directories in the directory specified")
        file = ReadRoot("tests/data_tests/background/", recursive=True)
        print(file)
        print(file.data)

        last_index = file.data.index.tolist()[-1]
        if type(last_index) != int:
            last_index = last_index[0]

        print('number of index times two: %d, lenght of df: %d' % (last_index * 2, len(file.data)) )


RootReadTestCase.main()
