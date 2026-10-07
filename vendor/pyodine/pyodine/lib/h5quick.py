# A few simple tools for quick reading of HDF5-data

import h5py
import sys
import logging


def h5print(h, level=0):
    """Recursively print the structure of a HDF5 file or group
    
    :param h: The h5py object.
    :type h: h5py file or group
    :param level: The level within the structure at which to start (needed to 
        print the whole structure recursively).
    :type level: int
    """
    # Check if there are already logging handles initialized somewhere - if no,
    # initialize a basic logger so that the output is actually printed
    if not logging.getLogger().hasHandlers():
        logging.basicConfig(stream=sys.stdout, level=logging.INFO, 
                            format='%(message)s')
    
    if level == 0:
        logging.info(h.name)
    for k in h.keys():
        logging.info('    ' * level + '└ ' + k)
        if isinstance(h[k], h5py.Group):
            h5print(h[k], level + 1)


def h5data(h):
    """Retrieve a HDF5 dataset or group
    
    If the handle is a group, return value is a dict of datasets (recursive).
    
    :param h: The h5py to return.
    :type h: h5py dataset or group
    
    :return: The data packed into a dictionary resembling the structure of the 
        h5py object.
    :rtype: dict or dataset
    """
    if isinstance(h, h5py.Dataset):
        return h[()]
    else:
        return {k: h5data(h[k]) for k in h.keys()}


def h5get(filename, item):
    """Retrieve a named item from a HDF5 file
    
    If the item is a group, return value is a dict of datasets (recursive).
    
    :param filename: The path to the HDF5 file.
    :type filename: str
    :param item: A key to the dataset or group of interest.
    :type item: str
    
    :return: The data, either as dictionary or a dataset.
    :rtype: dict or dataset
    """
    with h5py.File(filename, 'r') as h:
        return h5data(h[item])


# Helper function
def dict_to_group(my_dict, base_group, new_group_name):
    """Pack data into a h5py object
    
    Create a hdf5 group with name `new_group_name` in the existing
    group handle `base_group` (could be the root) and fill in
    named datasets from `my_dict`
    
    :param my_dict: The data to pack.
    :type my_dict: dict
    :param base_group: The group handle to pack the data into.
    :type base_group: h5py group
    :param new_group_name: The name of the group.
    :type new_group_name: str
    """
    group = base_group.create_group(new_group_name)
    for k in my_dict:
        # Only add the item if it is not None (otherwise errors are thrown)
        if my_dict[k] is not None:
            if isinstance(my_dict[k], dict):
                dict_to_group(my_dict[k], group, k)
            else:
                group[k] = my_dict[k]