'''
AUROC analysis
'''

import sys
sys.path.insert(1, '../Algorithms/')
sys.path.insert(1, 'Algorithms/')
sys.path.insert(1, '../Imports/')
import RandomWalk as rwr
import DiffusionKernel as dk
import PageRank as pr
import loader
import time
import numpy as np
import matplotlib.pyplot as plt

def roc_curve(result_vec, ground_truth_vec, name):
    TPR = []
    FPR = []
    for threshhold in range(len(result_vec)):
        tp = fp = fn = tn = 0
        for i in range(len(result_vec)):
            item = result_vec[i][0]
            if i <= threshhold and item in ground_truth_vec:
                tp += 1
            elif i <= threshhold and item not in ground_truth_vec:
                fp += 1
            elif item in ground_truth_vec:
                fn += 1
            else:
                tn += 1

        TPR.append(tp/(tp + fn))
        FPR.append(fp/(fp + tn))
    file_path = name + '2.png'
    area = np.trapz(TPR, FPR)

    file_path = "../Results/" + name + '.png'
    c = 'b'
    if 'pr' in name:
        c = 'r'
    elif 'dk' in name:
        c = 'g'
    plt.plot(FPR, TPR, c)

    plt.show()

def main():
    #Get file path choices
    #pathToPPINetworkFile = "../" + sys.argv[1]
    pathToPPINetworkFile = '../Data/9606.protein.links.v11.0.txt'

    # Get output vectors from each algorithm

    PPI_Network = loader.load_graph(pathToPPINetworkFile)  # load network
    ground_truth_files = ['../Data/MalaCard-protein-Endometriosis.diseasegenes.tsv', '../Data/MalaCard-protein-ischaemic-stroke.diseasegenes.tsv','../Data/MalaCard-protein-lymphoma.diseasegenes.tsv']
    file_paths = ['../Data/endometriosis-proteins.diseasegenes.tsv','../Data/lymphoma-proteins.diseasegenes.tsv', '../Data/ischaemic-stroke-proteins.diseasegenes.tsv']
    prior_paths = ['../Data/endometriosis-proteins-priors.diseasegenes.tsv','../Data/lymphoma-proteins-priors.diseasegenes.tsv', '../Data/ischaemic-stroke-proteins-priors.diseasegenes.tsv']
    names = ['endometriosis', 'lymphoma', 'ischaemic-stroke']

    for i in range(3):
        # building ground truth
        ground_truth_vec = []
        with open(ground_truth_files[i], 'r') as input:
            input = input.readlines()
            for line in input:
                protein = line.rstrip('\n')
                ground_truth_vec.append(protein)
        gene_file = open(file_paths[i], 'r')
        file_contents = gene_file.readLines()
        for line in input:
            protein = line.rstrip('\n')
            if protein not in ground_truth_vec:
                ground_truth_vec.append(protein)
        print(ground_truth_vec)
        # building start and priors vector

        start_vector = loader.load_start_vector(file_paths[i], PPI_Network)
        priors_vector = pr.load_priors(prior_paths[i], PPI_Network)

        #getting output from algorithms
        start_time = time.time()
        output_RWR = rwr.random_walk(PPI_Network, start_vector)
        end_time = time.time()
        print("time for rwr:", end_time - start_time)
        start_time = time.time()
        output_PR = pr.page_rank(PPI_Network, start_vector, priors_vector)
        end_time = time.time()
        print("time for pr:", end_time - start_time)

        start_time = time.time()
        output_DK = dk.diffusion_kernel(PPI_Network, start_vector)
        end_time = time.time()
        print("time for dk:", end_time - start_time)

        #building roc curves

        start_time = time.time()
        name = "rwr-" + names[i]
        roc_curve(output_RWR, ground_truth_vec, name)
        end_time = time.time()
        print("time for roc curve, rwr:", end_time -start_time)
        start_time = time.time()
        name = "pr-" + names[i]
        roc_curve(output_PR, ground_truth_vec, name)
        end_time = time.time()
        print("time for roc curve, pr:", end_time - start_time)
        start_time = time.time()

        start_time = time.time()
        name = "dk-" + names[i]
        roc_curve(output_DK, ground_truth_vec, name)
        end_time = time.time()
        print("time for roc curve, dk:", end_time - start_time)
        file_path = '../Results/' + names[i] + 'roc_curve.png'
        plt.savefig(file_path) #moved from roc_curve
        plt.clf() #moved from roc_curve
        print(colored("Done. ", "green") + "Plots have been saved as png files in the Results folder.")


if __name__ == '__main__':
    main()
