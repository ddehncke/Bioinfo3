import argparse
import sys
import os

# initialize the DP matrix

def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg  # return an open file handle

# helper function so i don't mess up the indices
rangeInclusive = lambda start, end: range(start, end + 1)
rangeExclusive = lambda start, end: range(start + 1, end)

def main():

    parser = argparse.ArgumentParser(description='imput fasta file')
    parser.add_argument('-i', dest="filename", metavar="FILE", help='input file', required=True,type=lambda x: is_valid_file(parser, x))

    args = parser.parse_args()

    score_AU = 1
    score_CG = 1
    score_GU = 1
    input_file = args.filename

    sequence = ''
    dic_sequence = {}
    with open('input_files/supplement01/test.fasta', 'r') as file_handle:
        for line in file_handle:
            line = line.strip()

            if line[0] == '>':
                continue
            sequence = line.strip()

    # transform seuence to dict because of indices
    for i in rangeInclusive(1,len(sequence)):
        dic_sequence[i] = sequence[i-1]

    dp_dic = {}
    traceback_dic = {}



    # initialize with 0
    for i in rangeInclusive(2,len(sequence)):
        dp_dic[i][i-1] = 0

    for i in rangeInclusive(1,len(sequence)):
        dp_dic[i][i] = 0

    dic_score = {}
    dic_score['A','U'] = score_AU
    dic_score['U', 'A'] = score_AU
    dic_score['C', 'G'] = score_CG
    dic_score['G', 'C'] = score_CG
    dic_score['G', 'U'] = score_GU
    dic_score['U', 'G'] = score_GU

    # initialization recursion
    for I in rangeInclusive(2,len(sequence)):
        for j in rangeInclusive(I,len(sequence)):
            i = j - I + 1

            case4_max = -1
            temp_trace = {}
            for k in rangeExclusive(i,j):
                if dp_dic[i,k] + dp_dic[k+1,j] > case4_max:
                    temp_trace = ((i,k),(k+1,j))

            case1 = dp_dic[i + 1, j]
            case2 = dp_dic[i, j - 1]
            case3 = dp_dic[i + 1, j - 1] + dic_score[dic_sequence[i], dic_sequence[j]]




main()