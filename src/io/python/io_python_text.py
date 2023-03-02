import os.path

def file_exists (file_name):
    """
    Checks whether the file exists.
    """
    return os.path.isfile(file_name)


def  get_length (input_file):
    """
    Return the number of lines in the input file.
    """
    with open(f'{input_file}.txt') as file:
        for i, line in enumerate(file):
            pass
    return i + 1


def read_list (input_file):
    """
    Return the contents of the line as a list.
    """
    with open(f'{input_file}.txt') as file:
        list = [int(line) for line in file]
    return list
