'''
AUROC analysis
'''

import sys
sys.path.insert(1, '../Algorithms/')
sys.path.insert(1, 'Algorithms/')
sys.path.insert(1, '../Scripts/')
import RandomWalk as rwr
import DiffusionKernel as dk
import PageRank as pr
import loader
import time
import matplotlib.pyplot as plt

def roc_curve(result_vec, ground_truth_vec, name):
    TPR = []
    FPR = []
    for threshhold in range(len(result_vec)):
        tp = 0
        fp = 0
        fn = 0
        tn = 0
        for i in range(len(result_vec)):
            item = result_vec[i]
            if i <= threshhold and item in ground_truth_vec:
                tp += 1
            elif i <= threshhold:
                fp += 1
            elif item in ground_truth_vec:
                fn += 1
            else:
                tn += 1
        print("true positive:", tp)
        print("false positive:", fp)
        print("true negative:", tn)
        print("false negative:", fn)
        TPR.append(tp/(tp + fn))
        FPR.append(fp/(fp + tn))
    file_path = name + '.png'
    plot = plt.plot(FPR, TPR)
    plot.savefig(file_path)

def main():
    # Get output vectors from each algorithm
    priors_file_path= '../Data/LymphomaProteinsPriors.tsv'
    PPI_Network = loader.load_graph('../Data/9606.protein.links.v11.0.txt') # load network
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
        print(ground_truth_vec)
        # building start and priors vector
        #
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

       # start_time = time.time()
       # output_DK = dk.diffusion_kernel(PPI_Network, start_vector)
       # end_time = time.time()
       # print("time for dk:", end_time - start_time)
        # building roc curves

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


if __name__ == '__main__':
    main()
