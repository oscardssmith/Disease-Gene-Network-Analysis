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
