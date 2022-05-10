from nndm_library import ReadFileBase


class BasicReadTestCase():
    def main():        
        # Read standard basic file
        print("Read single basic file that contains two rows: particle id and particle name")
        file = "tests/data_tests/basic_type/some_particle_names/quarks.txt"
        print("Reading %s" % file)
        print("*" * 50)
        file = ReadFileBase("tests/data_tests/basic_type/some_particle_names/quarks.txt")
        print(file.data)

        print("\n\nRead all basic files like the one before, but in a given directory")
        file = "tests/data_tests/basic_type/some_particle_names/"
        print("Reading directory %s" % file)
        print("*" * 50)
        file = ReadFileBase("tests/data_tests/basic_type/some_particle_names/")
        print(file.data)

        print("\n\nRead all basic files like the one before, but all found recursively in a given directory")
        file = "tests/data_tests/basic_type/"
        print("Reading all files found recursively in directory %s" % file)
        print("*" * 50)
        file = ReadFileBase("tests/data_tests/basic_type/", recursive=True, relabel_events=True)
        print(file.data)

BasicReadTestCase.main()
