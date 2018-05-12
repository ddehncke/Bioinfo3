import argparse
import math
from Bio.PDB import *
import matplotlib.pyplot as plt
import numpy as np



def main():

    parser = argparse.ArgumentParser(description='pdb files')
    parser.add_argument('-i', dest="file", metavar="FILE", help='input file', required=True)
    parser.add_argument('-o', dest="outputfile", metavar="outputFile", help='save to output file', required=True)

    args = parser.parse_args()

    filename = args.file
    output_file = args.outputfile

    parser = PDBParser()
    stack_residues = []
    structure = parser.get_structure('temp', filename)

    phi = []
    psi = []

    wait_three_residues = 1
    i = 0

    for model in structure:

        for chain in model:

            for residue in chain:

                # need to skip residues that are missing some atoms, because else it is not possible to compute the
                # angles. Maybe own parsing would lead to better results.
                if 'O' not in residue or 'N' not in residue or 'CA' not in residue or 'C' not in residue:
                    # print(residue['CA'])
                    i += 1
                    # print(i)
                    stack_residues = []
                    wait_three_residues = 1
                    continue

                stack_residues.append(residue)

                # first have to read in 3 residues, to compute the angles.
                if wait_three_residues >= 3:
                    phi.append(getAngle(stack_residues[0]['C'].get_vector(), stack_residues[1]['N'].get_vector(), stack_residues[1]['CA'].get_vector(),
                                        stack_residues[2]['C'].get_vector()))

                    psi.append(getAngle(stack_residues[0]['N'].get_vector(), stack_residues[1]['CA'].get_vector(),
                                        stack_residues[2]['C'].get_vector(),
                                        stack_residues[2]['N'].get_vector()))

                    stack_residues.pop(0)
                wait_three_residues += 1

        plt.scatter(phi, psi)
        plt.xlabel('phi-angle')
        plt.ylabel('psi-angle')
        plt.title(filename.split('/')[1][:-4])
        plt.savefig(output_file)

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