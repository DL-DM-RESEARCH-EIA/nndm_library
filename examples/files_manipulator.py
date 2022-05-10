from nndm_library import FilesManipulator

class BasicFileManipulatorTestCase():
    def main():
        # Reading files from signal
        print("\n" * 2)
        print("*" * 50)
        print("Read momenta of set of files")
        # Scan data
        file_manipulator = FilesManipulator("tests/data_tests/signal/*.lhe", var_of_interest=['e', 'px', 'py', 'pz'])
        file_manipulator.fill_up_scan()
        print(file_manipulator.scan)

        print("\n" * 2)
        print("*" * 50)
        print("Save scan of files")
        from os.path import exists
        save_file = "tests/data_tests/signal/signal.pickle"
        file_manipulator.save_scan(save_file)
        print("file " + save_file + " exists?", "yes" if exists(save_file) else "no")

BasicFileManipulatorTestCase.main()