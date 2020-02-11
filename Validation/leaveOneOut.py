"""
This program assumes the following:
Each algorithm can be called with only two parameters:
    - a networkx graph containing the full PPI network
    - a list of the names of known disease genes

Each algorithm produces results in a 2 dimensional array where each item is a tuple containing a gene name and it's ranking/probability,
 and the list is in order of probability, high to low, so that the index is equivalent to the ranking:
    - [[geneName, probability], [geneName, probability], [geneName, probability], [geneName, probability]]

"""
import sys
import os
sys.path.insert(1, '../Algorithms/')
sys.path.insert(1, 'Algorithms/')
sys.path.insert(1, '../Imports/')
import RandomWalk as rwr
import DiffusionKernel as dk
import PageRank as pr
from CacheUtils import compute_if_not_cached
from loader import load_graph, load_start_vector
import time
import numpy as np
import networkx as nx




def leave_one_out(function, diseaseGeneFilePath, PPI_Network, param):
    print("Starting leaveOneOut function")

    # building list of disease genes
    diseaseGeneFile = open(diseaseGeneFilePath, 'r')
    allDiseaseGenes = diseaseGeneFile.read().splitlines()
    diseaseGeneFile.close()

    numDiseaseGenes = len(allDiseaseGenes)
    rankThreshhold = 150

    numGenesNotFound = 0

    graph_nodes = list(PPI_Network.nodes())
    startVector = load_start_vector(diseaseGeneFilePath, PPI_Network)
    startVector = (numDiseaseGenes/(numDiseaseGenes - 1)) * startVector
    # skipping 
    for index, skipGene in enumerate(allDiseaseGenes):
        
        # find the skip gene in the start vector, make it zero
        index = graph_nodes.index(skipGene)
        startVector[index] = 0
        priors_vector = np.zeros(PPI_Network.number_of_nodes())
        if function == pr.page_rank:
            priors_file_path = find_priors_file(diseaseGeneFilePath)
            priors_vector = pr.load_priors(priors_file_path, PPI_Network)
            priors_vector[index] = 0

        #run algorithm using modified disease gene file
        startTime = time.time()
        output= []
        if function == pr.page_rank:
            output = function(PPI_Network, startVector, priors_vector, param)
        else:
            output = function(PPI_Network, startVector, param)
        endTime = time.time()
        print("finished algorithm. Time elapsed:", endTime - startTime)

        #find the predicted probability of the omitted gene and add it to the current sum
        startTime = time.time()
        foundGene = False
        for i in range(rankThreshhold):
            if output[i][0] == skipGene:
                foundGene = True
                print("Found the gene: ", skipGene, "at rank: ", i)
                break
        if not foundGene:
            numGenesNotFound +=1

        endTime = time.time()

    print("------------------------\nFinished running algorithm with all disease genes left out\nCalculating mean squared difference")
    print("Num genes not found for this run of leave one out: ", numGenesNotFound)
    #Find average of all squared differences
    percentCorrectlyRankedGenes = 1 - numGenesNotFound/numDiseaseGenes
    return percentCorrectlyRankedGenes




def get_files_in_directory(path):
    return [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]

def find_priors_file(diseaseGeneFilePath):
    targetName = diseaseGeneFilePath.split(".")[0]
    for f in get_files_in_directory("Data/"):
        if 'priors' in f.split('.') and targetName in f.split('.'):
            return f
    cprint("No priors file exists for your specified disease gene set.\nCannot run leave-one-out with PageRank on this disease gene set.", "red")
    sys.exit(0)


def main():

    algorithm = sys.argv[1]
    pathToPPINetworkFile = sys.argv[2]
    pathToDiseaseGeneFile = sys.argv[3]
    param = float(sys.argv[4])
    outputFile = sys.argv[5]

    print("loading data from files..")
    ppiGraph = compute_if_not_cached(load_graph, pathToPPINetworkFile, fileName="ppiGraph")
    

    if algorithm == "Algorithms/DiffusionKernel.py":
        function = dk.diffusion_kernel
    elif algorithm == "Algorithms/PageRank.py":
        function = pr.page_rank
    elif algorithm == "Algorithms/RandomWalk.py":
        function = rwr.random_walk
    else:
        function = None


    result = leave_one_out(function, pathToDiseaseGeneFile, ppiGraph, param)

    print("Saving results to:", outputFile)
    with open(outputFile, "w") as of:
        of.write("Leave-one-out Validation Results:\n\nAlgorithm:\t\t{0}\nPPI Graph:\t\t{1}\nDisease Genes:\t\t{2}\nPercentage Correctly Found Genes:\t\t{3}%\n\n".format(algorithm, pathToPPINetworkFile, pathToDiseaseGeneFile, result*100))


if __name__ == '__main__':
    main()
