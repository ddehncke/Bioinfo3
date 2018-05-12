#!/usr/bin/env python

import argparse
import math
import sys
import os
from collections import OrderedDict
from operator import itemgetter
from Bio.PDB import *
import matplotlib.pyplot as plt
import numpy as np
import plotly.plotly as py


def main():

    parser = argparse.ArgumentParser(description='pdb files')
    parser.add_argument('-i', dest="folder", metavar="FILE", help='input folder', required=False)

    args = parser.parse_args()

    all_paths = os.listdir(args.folder)

    # for filename in all_paths:
    #     path_to_current_file = os.path.join(args.folder, filename)
    #     parser = PDBParser()
    #     structure = parser.get_structure('temp', path_to_current_file)
    #
    #

    # phi
    C1 = np.array([8.326, 10.351, 0.000])
    N = np.array([9.000, 9.000, 0.000])
    Ca = np.array([10.325, 9.000, 0.000])
    C2 = np.array([11.096, 7.766, 0.000])

    b_1 = N - C1
    b_2 = Ca - N
    b_3 = C2 - Ca

    q1 =np.cross(b_1, b_2)
    n_1 = np.cross(b_1, b_2) / np.linalg.norm(np.cross(b_1, b_2))
    n_2 = np.cross(b_2, b_3) / np.linalg.norm(np.cross(b_2, b_3))

    m_1 = np.cross(n_1, b_2/np.linalg.norm(b_2))

    x = np.dot(n_1, n_2)
    y = np.dot(m_1, n_2)


    print(np.degrees(-math.atan2(y, x)))

    print(calc_dihedral(C1,N,Ca,C2))

main()