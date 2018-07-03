'''
Library information for RAD-Seq analysis.
'''

import re
import sys
import itertools
import os.path
from glob import glob
from collections import defaultdict


class Library(object):
    '''Class for library information'''
    def __init__(self, name, lib_dir, r1, r2, barcodes_filename):
        self.name = name
        self.lib_dir = lib_dir
        self.lib_id = os.path.basename(self.lib_dir)
        self.r1 = os.path.join(lib_dir, r1)
        self.r2 = os.path.join(lib_dir, r2)
        self.barcodes_filename = os.path.join(lib_dir, barcodes_filename)
        self.files = [self.r1, self.r2, self.barcodes_filename]
        self.samples = self.parse_sample_list()

    def parse_sample_list(self):
        '''Get sample list from barcodes file'''
        with open(self.barcodes_filename, "r") as f:
            sample_barcodes = f.read().strip().split("\n")
        samples = [x.split("\t")[2] for x in sample_barcodes]
        return samples

    def __str__(self):
        # TODO
        str = ""
        return str

def parse_libraries(libraries):
    '''
    Create Library object for each library and return list of libraries
    '''
    libraries_list = []
    for lib_name, d in libraries.items():
        new_lib = Library(name=lib_name,
                          lib_dir=d["lib_dir"],
                          r1=d["r1"],
                          r2=d["r2"],
                          barcodes_filename=d["barcodes"])
        libraries_list.append(new_lib)
    return libraries_list
