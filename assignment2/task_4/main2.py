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

    for filename in all_paths:
        path_to_current_file = os.path.join(args.folder, filename)
        parser = PDBParser()
        structure = parser.get_structure('temp', path_to_current_file)

        phi = []
        psi = []

        temp_residue = None

        for model in structure:
            for chain in model:
                for residue in chain:

                    # need to skip residues that are incorrect
                    if 'O' not in residue or 'N' not in residue:
                        temp_residue = None
                        continue
                    if temp_residue is not None:
                        # TODO Es müssen drei residuen betrachtet werden. Einfach eine temp_residue mehr, da eh alles
                        # auf stack
                        phi.append(getAngle(temp_residue['C'], temp_residue['N'], residue['Ca'], residue['C']))
                        psi.append(getAngle(temp_residue['N'], temp_residue['Ca'], residue['Ca'], residue['C']))



                    # temp_dist.append(np.linalg.norm(stack_O.pop() - residue['N'].get_vector()))

                    print(residue['N'].get_vector())







    # A1 = np.array([8.326, 10.351, 0.000])
    # A2 = np.array([9.000, 9.000, 0.000])
    # A3 = np.array([10.325, 9.000, 0.000])
    # A4 = np.array([11.096, 7.766, 0.000])

    phi = getAngle(C1, N, Ca, C2)
    psi = getAngle(N, Ca, C, N)



def getAngle(A1, A2, A3, A4):


    b_1 = A2 - A1
    b_2 = A3 - A2
    b_3 = A4 - A3

    q1 =np.cross(b_1, b_2)
    n_1 = np.cross(b_1, b_2) / np.linalg.norm(np.cross(b_1, b_2))
    n_2 = np.cross(b_2, b_3) / np.linalg.norm(np.cross(b_2, b_3))

    m_1 = np.cross(n_1, b_2/np.linalg.norm(b_2))

    x = np.dot(n_1, n_2)
    y = np.dot(m_1, n_2)


    print(np.degrees(-math.atan2(y, x)))

main()