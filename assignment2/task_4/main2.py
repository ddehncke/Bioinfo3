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
    parser.add_argument('-i', dest="folder", metavar="FILE", help='input folder', required=True)
    parser.add_argument('-o', dest="outputfile", metavar="outputFile", help='save to output file', required=True)

    args = parser.parse_args()


    all_paths = os.listdir(args.folder)
    output_file = args.outputfile

    for filename in all_paths:
        path_to_current_file = os.path.join(args.folder, filename)
        parser = PDBParser()
        stack_residues = []
        structure = parser.get_structure('temp', path_to_current_file)


        phi = []
        psi = []
        phi_psi_tuple = []
        wait_three_residues = 1
        i = 0

        for model in structure:

            for chain in model:

                for residue in chain:

                    # need to skip residues that are incorrect
                    if 'O' not in residue or 'N' not in residue or 'CA' not in residue or 'C' not in residue:
                        # print(residue['CA'])
                        i += 1
                        # print(i)
                        stack_residues = []
                        wait_three_residues = 1
                        continue

                    stack_residues.append(residue)

                    if wait_three_residues >= 3:
                        # TODO Es müssen drei residuen betrachtet werden. Einfach eine temp_residue mehr, da eh alles
                        # auf stack
                        phi.append(getAngle(stack_residues[0]['C'].get_vector(), stack_residues[1]['N'].get_vector(), stack_residues[1]['CA'].get_vector(),
                                            stack_residues[2]['C'].get_vector()))

                        psi.append(getAngle(stack_residues[0]['N'].get_vector(), stack_residues[1]['CA'].get_vector(),
                                            stack_residues[2]['C'].get_vector(),
                                            stack_residues[2]['N'].get_vector()))

                        stack_residues.pop(0)
                    wait_three_residues += 1
                        # psi.append(getAngle(temp_residue['N'], temp_residue['Ca'], residue['Ca'], residue['C']))

        plt.scatter(phi, psi)
        plt.savefig(output_file)

        # print(phi)






    # A1 = np.array([8.326, 10.351, 0.000])
    # A2 = np.array([9.000, 9.000, 0.000])
    # A3 = np.array([10.325, 9.000, 0.000])
    # A4 = np.array([11.096, 7.766, 0.000])

    # phi = getAngle(C1, N, Ca, C2)
    # psi = getAngle(N, Ca, C, N)



def getAngle(A1, A2, A3, A4):


    br_1 = A2 - A1
    br_2 = A3 - A2
    br_3 = A4 - A3

    b_1 = np.array([br_1[0], br_1[1], br_1[2]])
    b_2 = np.array([br_2[0], br_2[1], br_2[2]])
    b_3 = np.array([br_3[0], br_3[1], br_3[2]])

    q1 =np.cross(np.array(b_1), np.array(b_2))
    n_1 = np.cross(b_1, b_2) / np.linalg.norm(np.cross(b_1, b_2))
    n_2 = np.cross(b_2, b_3) / np.linalg.norm(np.cross(b_2, b_3))

    m_1 = np.cross(n_1, b_2/np.linalg.norm(b_2))

    x = np.dot(n_1, n_2)
    y = np.dot(m_1, n_2)


    return np.degrees(-math.atan2(y, x))

main()