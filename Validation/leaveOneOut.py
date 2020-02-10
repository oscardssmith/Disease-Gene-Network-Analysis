"""
This program assumes the following:
Each algorithm can be called with only two parameters:
    - a networkx graph containing the full PPI network
    - a list of the names of known disease genes

Each algorithm produces results in a 2 dimensional array where each item is a tuple containing a gene name and it's ranking/probability,
 and the list is in order of probability, high to low, so that the index is equivalent to the ranking:
    - [[geneName, probability], [geneName, probability], [geneName, probability], [geneName, probability]]
    - TALK ABOUT THIS -- algorithms need to output pairings of gene names and rankings, it is not the job of validation to interpret output and match up numbers to genes


"""
import sys
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

    #print("finished initialization, starting disease gene loop")
    graph_nodes = list(PPI_Network.nodes())
    startVector = load_start_vector(diseaseGeneFilePath, PPI_Network)
    startVector = (numDiseaseGenes/(numDiseaseGenes - 1)) * startVector
    # skipping 
    for index, skipGene in enumerate(allDiseaseGenes):
        #print("looping! skipping gene: ", skipGene)
        # loading the start vector
        
        # find the skip gene in the start vector, make it zero
        index = graph_nodes.index(skipGene)
        startVector[index] = 0
        priors_vector = np.zeros(PPI_Network.number_of_nodes())
        if function == pr.page_rank:
            priors_file_path = find_priors_file(diseaseGeneFilePath)
            priors_vector = pr.load_priors(priors_file_path, PPI_Network)
            priors_vector[index] = 0

        #run algorithm using modified disease gene file
        #print("calling algorithm")
        startTime = time.time()
        output= []
        if function == pr.page_rank: #is this right?
            #print("Using PageRank now")
            output = function(PPI_Network, startVector, priors_vector, param)
        else:
            #print("Using RWR now")
            output = function(PPI_Network, startVector, param)
        endTime = time.time()
        print("finished algorithm. Time elapsed:", endTime - startTime)

        #find the predicted probability of the omitted gene and add it to the current sum
        #print("finding skipgene predicted probability")
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
        #print("added skipGene predictedProbability to squaredDifferenceSum\nTime elapsed:", endTime - startTime)


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

    print(function)

    result = leave_one_out(function, pathToDiseaseGeneFile, ppiGraph, param)

    with open(outputFile, "w") as of:
        of.write("Leave-one-out Validation Results:\n\nAlgorithm:\t\t{0}\nPPI Graph:\t\t{1}\nDisease Genes:\t\t{2}\nPercentage Correctly Found Genes:\t\t{3}%\n\n".format(algorithm, pathToPPINetworkFile, pathToDiseaseGeneFile, result*100))













    # # # load the full PPI Network
    # totalStartTime = time.time()
    # print("Starting leave one out validation.")
    # print("Loading graph")
    # PPI_Network = load_PPI_Network('../Data/9606.protein.links.v11.0.txt')
    # print("Loaded graph")
    # file_paths = ['../Data/endometriosis-proteins.diseasegenes.tsv','../Data/lymphoma-proteins.diseasegenes.tsv', '../Data/ischaemic-stroke-proteins.diseasegenes.tsv']
    # prior_paths = ['../Data/endometriosis-proteins-priors.diseasegenes.tsv','../Data/lymphoma-proteins-priors.diseasegenes.tsv', '../Data/ischaemic-stroke-proteins-priors.diseasegenes.tsv']
    # rwr_pr_params = [0.2, 0.4, 0.6, 0.8]
    # dk_params = [0.4, 0.6, 0.8, 1.0]

    # for file_index in range(3):
    #     for param_index in range(4):
    #         print("--------------------DISEASE GENE FILE: ", file_paths[file_index], "PARAMETER: ", rwr_pr_params[param_index], "---------------------------------")
    #         print("---------------RWR----------------------")
    #         result_rwr = leaveOneOut(rwr.random_walk, file_paths[file_index], PPI_Network, rwr_pr_params[param_index])
    #         print("percentage of genes improperly predicted for RWR:", result_rwr)
    #         print("---------------PR----------------------")
    #         result_pr = leaveOneOut(pr.page_rank, file_paths[file_index], PPI_Network, rwr_pr_params[param_index], prior_paths[file_index])
    #         print("percentage of genes improperly predicted for PR:", result_pr)
    #         print("---------------DK----------------------")
    #         result_dk = leaveOneOut(dk.diffusion_kernel, file_paths[file_index], PPI_Network, dk_params[param_index], prior_paths[file_index])
    #         print("percentage of genes improperly predicted for DK:", result_dk)

    # # result = leaveOneOut(pr.PageRank, 'diseaseGeneFile', PPI_Network, priors_file_path)
    # # print("Mean squared difference for PageRank:", result)

    # totalEndTime = time.time()
    # print("Finished leave one out validation!\nTotal time in hours:", (totalEndTime - totalStartTime)/3600)




if __name__ == '__main__':
    main()
