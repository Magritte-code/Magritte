import numpy as np
import h5py  as hp


def read_length (io_file, file_name):
    """
    Return the number of lines in the input file.
    """
    with hp.File (io_file, 'r') as file:

        # print('In the file? here are the keys:')
        # print(file.keys())

        try:
            # print('fine')
            # Try to open the object
            object = file [file_name]
            # Check if it is a Dataset

            # print ('length = ', object.len())
            if isinstance (object, hp.Dataset):
                return object.len()
        except:

            # print ('not supposed to be here')
            # Get name of object we need to count
            object_name = file_name.split('/')[-1]
            # Get containing group
            group_name = file_name[:-len(object_name)]
            # Count occurences
            length = 0
            for key in file[group_name].keys():
                if (object_name in key):
                    length += 1
            return length
    # Error if not yet returned
    raise ValueError ('file_name is no Group or Dataset.')


def read_width (io_file, file_name):
    """
    Return the number of columns in the input file.
    """
    with hp.File (io_file, 'r') as file:
        try:
            # Try to open the object
            object = file [file_name]
            # Check if it is a Dataset
            if isinstance (object, hp.Dataset):
                return object.shape[1]
        except:
            # Get name of object we need to count
            object_name = file_name.split('/')[-1]
            # Get containing group
            group_name = file_name[:-len(object_name)]
            # Count occurences
            length = 0
            for key in file[group_name].keys():
                if (object_name in key):
                    length += 1
            return length
    # Error if not yet returned
    raise ValueError ('file_name is no Group or Dataset.')


def read_attribute (io_file, file_name):
    """
    Return the contents of the attribute
    """
    with hp.File (io_file, 'r') as file:
        object    = file_name.split('.')[0]
        attribute = file_name.split('.')[1]
        if object != '':
            return file[object].attrs[attribute]
        else:
            return file.attrs[attribute]


def write_attribute (io_file, file_name, data):
    """
    Write the data to the attribute
    """
    with hp.File (io_file) as file:
        object    = file_name.split('.')[0]
        attribute = file_name.split('.')[1]
        # Make sure all groups exists, if not create them
        # NOTE: ASSUMES THAT WORD IS WRITTEN TO A GROUP !
        group = ''
        for g in object.split('/'):
            group += f'/{g}'
            file.require_group (group)
        if object != '':
            file[object].attrs[attribute] = data
        else:
            file.attrs[attribute] = data


def read_number (io_file, file_name):
    """
    Return the contents of the attribute
    """
    return read_attribute(io_file, file_name)


def write_number (io_file, file_name, data):
    """
    Write the data to the attribute
    """
    write_attribute (io_file, file_name, data)


def get_element_type (l):
    if isinstance(l, list):
        return get_element_type(l[0])
    else:
        return type(l)


def read_array (io_file, file_name):
    """
    Return the contents of the data array.
    """
    with hp.File (io_file, 'r') as file:
        if (file_name in file):
            return np.array (file.get (file_name))


def write_array (io_file, file_name, data):
    """
    Write the contents to the data array.
    """
    with hp.File (io_file) as file:
        # print ('Writing array to HDF5 file...')
        # print (io_file, file_name)
        # Delete if dataset already exists
        try:
            # print("deleting ", file_name)
            del file[file_name]
        except:
            # print("Nothing to delete")
            pass
        # Make sure all groups exists, if not create them
        # NOTE: ASSUMES THAT DATA IS WRITTEN TO A DATASET
        group = ''
        for g in file_name.split('/')[:-1]:
            group += f'/{g}'
            file.require_group (group)
            # print("required ", group)
        # Write dataset
        try:
            # print('Creating dataset...')
            if (get_element_type(data) == str):
                # print(data)
                file.create_dataset (name=file_name, data=np.array(data, dtype='S'))
            else:
                file.create_dataset (name=file_name, data=data)
        except:
            print ("failed to write ", file_name)
