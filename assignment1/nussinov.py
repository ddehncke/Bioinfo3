import argparse
import sys
import os

def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg

# helper function so i don't mess up the indices
rangeInclusive = lambda start, end: range(start, end + 1)
rangeExclusive = lambda start, end: range(start + 1, end)

# reads the fasta file. The dp dictionary is initialised, together with a dictionary to store the scores.
def computeNussinov(score_AU,score_CG,score_GU,min_loop_length,input_file):
    sequence = ''

    with open(input_file, 'r') as file_handle:
        for line in file_handle:
            line = line.strip()
            if line[0] == '>':
                continue
            sequence = line.strip()

    # transform seuence to dict because of indices
    dic_sequence = {}
    for i in rangeInclusive(1,len(sequence)):
        dic_sequence[i] = sequence[i-1]

    dp_dic = {}
    traceback_dic = {}

    # create 2 dimensional dictionatries for dp and traceback
    for i in rangeInclusive(1,len(sequence)):
        dp_dic[i] = {}

    # fill with 0
    for i in rangeInclusive(2,len(sequence)):
        dp_dic[i][i-1] = 0
    for i in rangeInclusive(1,len(sequence)):
        dp_dic[i][i] = 0

    dic_BaseScore = {('A','U'):score_AU, ('U','A'):score_AU,
                 ('C','G'):score_CG, ('G','C'):score_CG,
                 ('G','U'):score_GU, ('U','G'):score_GU,
                 ('U', 'U'): 0, ('G', 'G'): 0,
                 ('A', 'A'): 0, ('C', 'C'): 0,
                 ('U', 'C'): 0, ('C', 'U'): 0,
                 ('A', 'C'): 0, ('C', 'A'): 0,
                 ('A', 'G'): 0, ('G', 'A'): 0,
                     'minimumLoopLength':min_loop_length}

    return fillMatrices(dp_dic, dic_BaseScore, dic_sequence)

# filling the dp matrix and then computing the traceback
def fillMatrices(dp_dic, dic_BaseScore, dic_sequence):

    def getPairScore(i, j):
        return dic_BaseScore[(dic_sequence[i], dic_sequence[j])]

    # initialization of dp matrix
    for I in rangeInclusive(2, len(dic_sequence)):
        for j in rangeInclusive(I, len(dic_sequence)):
            i = j - I + 1

            case1 = dp_dic[i + 1][j]
            case2 = dp_dic[i][j - 1]
            case3 = dp_dic[i + 1][j - 1] + getPairScore(i, j) if abs(i - j) > dic_BaseScore['minimumLoopLength'] else 0

            case4 = 0
            for k in rangeExclusive(i,j):
                if dp_dic[i][k] + dp_dic[k+1][j] > case4:
                    case4 = dp_dic[i][k] + dp_dic[k+1][j]

            # if case1 > case2 and case1 > case3 and case1 > case4:
            #     dp_dic[i][j] = case1
            # elif case2 > case3 and case2 > case4:
            #     dp_dic[i][j] = case2
            # elif case3 > case4:
            #     dp_dic[i][j] = case3
            # else:
            #     dp_dic[i][j] = case4

            dp_dic[i][j] = max(case1, case2, case3, case4)


    matrix_BasePairs = [-1 for x in range(len(dic_sequence)+1)]

    # recursively go through the matrix and add the values to a traceback array
    def traceback(i,j):
        if i <= j:
            if dp_dic[i][j] == dp_dic[i+1][j]:
                matrix_BasePairs[i] = 0
                traceback(i + 1, j)
            elif dp_dic[i][j] == dp_dic[i][j-1]:
                matrix_BasePairs[j] = 0
                traceback(i, j - 1)
            elif dp_dic[i][j] == dp_dic[i + 1][j - 1] + getPairScore(i, j):
                matrix_BasePairs[j] = i
                matrix_BasePairs[i] = j
                traceback(i + 1, j - 1)
            else:
                for k in rangeInclusive(i+1,j-1):
                    if dp_dic[i][j] == dp_dic[i][k] + dp_dic[k+1][j]:
                        traceback(i, k)
                        traceback(k + 1,j)
                        break

    traceback(1, len(dic_sequence))
    bpseq = computeBPSequence(dic_sequence, matrix_BasePairs)

    return bpseq, dic_sequence, matrix_BasePairs

# print the matrix for error testing
def printMatrix(dp_dic):
    for k, v in dp_dic.items():
        print(v)

# constructs the bpseq format
def computeBPSequence(dic_sequence, matrix_BasePairs):
    bpseq = ''

    for i in range(1, len(matrix_BasePairs)):
        if matrix_BasePairs[i] == 0:
            bpseq+='.'
        elif i < matrix_BasePairs[i]:
            bpseq += '('
        else:
            bpseq += ')'

    return bpseq

# handle the input arguments
def inputHandler():
    parser = argparse.ArgumentParser(description='imput fasta file')
    parser.add_argument('-i', dest="filename", metavar="FILE", help='input file', required=True,
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("--min-loop-length", type=int, dest='min_loop_length', default=5)
    parser.add_argument("-score_AU", type=int, dest='score_AU', default=10)
    parser.add_argument("-score_CG", type=int, dest='score_CG', default=1)
    parser.add_argument("-score_GU", type=int, dest='score_GU', default=0)

    args = parser.parse_args()

    score_AU = args.score_AU
    score_CG = args.score_CG
    score_GU = args.score_GU
    min_loop_length = args.min_loop_length
    input_file = args.filename
    return score_AU,score_CG,score_GU,min_loop_length,input_file

# tests some parameters to find the ones for which the prediction fits the original structure
def findBestParameters(input_file):
    # method for finding best values for tRNA structure prediction
    best_score, best_score_AU, best_score_CG, best_score_GU, best_min_loop_length = 0,0,0,0,0
    correctSecondaryStructure = '(((((((..((((.........)))).(((((.......))))).....(((((.......)))))))))))).'
    bestSecondaryStructure = ''
    min,max = 1,8
    score = 0

    for au in range(min,max):
        for cg in range(min,max):
            for gu in range(min,max):
                for minLoopLength in range(min,max):
                    predictedStructure,dic_sequence, matrix_BasePairs = computeNussinov(au, cg, gu, minLoopLength, input_file)

                    score = 0

                    for i in range(len(predictedStructure)):
                        if predictedStructure[i] == correctSecondaryStructure[i]:
                            score = score + 1
                    if score > best_score:
                        best_score, best_score_AU, best_score_CG, best_score_GU, best_min_loop_length = score, au, cg, gu, minLoopLength
                        bestSecondaryStructure = predictedStructure

    print('The best score in comparisson to the original structure was achieved with a identity of',best_score)
    print('with the parameters au=', best_score_AU,'gc= ', best_score_CG,'gu=', best_score_GU, 'and a minimum loop length of', best_min_loop_length)
    print('first one is the correct structure')
    print(correctSecondaryStructure)
    print(predictedStructure)


def main():

    score_AU, score_CG, score_GU, min_loop_length, input_file = inputHandler()

    bpseq, dic_sequence, matrix_BasePairs  = computeNussinov(score_AU, score_CG, score_GU, min_loop_length, input_file)
    # printing the dict that contains the sequence
    print(''.join([dic_sequence[i] for i in rangeInclusive(1, len(dic_sequence))]))
    print(bpseq)
    for i in range(1, len(matrix_BasePairs)):
        print(i, dic_sequence[i], matrix_BasePairs[i])

    #findBestParameters(input_file)


main()