'''
AUROC analysis
'''

import sys
sys.path.insert(1, '../RWR/')
sys.path.insert(1, '../DiffusionKernel/')
sys.path.insert(1, '../PageRank/')
sys.path.insert(1, '../Scripts/')
import RandomWalk as rwr
import DiffusionKernel as dk
import PageRank as pr
import loader
import time
import matplotlib as plt

def roc_curve(result_vec, ground_truth_vec):
    TPR = []
    FPR = []
    for threshhold in range(len(result_vec)):
        tp = 0
        fp = 0
        fn = 0
        tn = 0
        for i in range(len(result_vec)):
            item = result_vec[i]
            if i < threshhold and item in ground_truth_vec:
                tp += 1
            elif i < threshhold:
                fp += 1
            elif item in ground_truth_vec:
                fn += 1
            else:
                tn += 1
        TPR.append(tp/(tp + fn))
        FPR.append(fp/(fp + tn))
    plt.pyplot.plot(FPR, TPR)

def main():
    # Get output vectors from each algorithm
    priors_file_path= '../Data/LymphomaProteinsPriors.tsv'
    PPI_Network = load_PPI_Network('../Data/9606.protein.links.v11.0.txt') # load network
    startVector = loader.load_start_vector('../Data/LymphomaProteins.tsv', PPI_Network)
    priors_vector = pr.load_priors(priors_file_path, PPI_Network)

    output_RWR = rwr.random_walk(PPI_Network, startVector)
    output_PR = pr.PageRank(PPI_Network, startVector, priors_vector)
    output_DK = dk.diffusion_kernel(PPI_Network, startVector)



if __name__ == '__main__':
    main()
