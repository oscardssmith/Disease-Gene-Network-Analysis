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




def leave_one_out(function, diseaseGeneFilePath, PPI_Network, param):
    print("Starting leaveOneOut function")

    # building list of disease genes
    diseaseGeneFile = open(diseaseGeneFilePath, 'r')
    allDiseaseGenes = diseaseGeneFile.read().splitlines()
    diseaseGeneFile.close()
   # print(allDiseaseGenes)
    numDiseaseGenes = len(allDiseaseGenes)
    rankThreshhold = 150

    numGenesNotFound = 0

    degree_list = [] #remove after graph is made (kate)
    in_out_list = [] #remove after graph is made (kate)

    graph_nodes = list(PPI_Network.nodes())
    startVector = load_start_vector(diseaseGeneFilePath, PPI_Network)
    startVector = (numDiseaseGenes/(numDiseaseGenes - 1)) * startVector
    # skipping
    for index, skipGene in enumerate(allDiseaseGenes):

        # find the skip gene in the start vector, make it zero
        index = graph_nodes.index(skipGene)
        node_degree = PPI_Network.degree(skipGene) #remove after graph is made (kate)
        degree_list.append(node_degree) #remove after graph is made (kate)
        newStartVector = startVector.copy()
        newStartVector[index] = 0
     #   startVector[index] = 0 
        priors_vector = np.zeros(PPI_Network.number_of_nodes())
        if function == pr.page_rank:
            priors_file_path = find_priors_file(diseaseGeneFilePath)
            priors_vector = pr.load_priors(priors_file_path, PPI_Network)
            priors_vector[index] = 0

        #run algorithm using modified disease gene file
        startTime = time.time()
        output= []
        if function == pr.page_rank:
            output = function(PPI_Network, newStartVector, priors_vector, param)
        else:
            print("sum of start vector:", np.sum(startVector))
            output = function(PPI_Network, newStartVector, param)
        endTime = time.time()
        print("finished algorithm. Time elapsed:", endTime - startTime)

        #find the predicted probability of the omitted gene and add it to the current sum
        startTime = time.time()
        foundGene = False
        for i in range(rankThreshhold):
            if output[i][0] == skipGene:
                foundGene = True
                print("Found the gene: ", skipGene, "at rank: ", i)
                in_out_list.append(1) #remove after graph is made (kate)
                break
        if not foundGene:
            numGenesNotFound +=1
            in_out_list.append(-1) #remove after graph is made (kate)

        endTime = time.time()
    
    # write the results of leave one out to a file
    disease_name = diseaseGeneFilePath.split(".")[0]
    output_name = "leave_one_out_1" + disease_name[5:]
    if function == pr.page_rank:
        output_name = output_name + "_pr.tsv"
    elif function == dk.diffusion_kernel:
        output_name = output_name + "_dk.tsv"
    elif function == rwr.random_walk:
        output_name = output_name + "_rwr.tsv"
    with open(output_name, "w") as output:
        for i in range(len(allDiseaseGenes)):
            output_string = allDiseaseGenes[i] + "\t" +str(degree_list[i]) + "\t" +str(in_out_list[i]) + "\n"
            output.write(output_string)
    
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
        f = "Data/" + f
        if 'priors' in f.split('.') and targetName in f.split('.'):
            return f
    print("No priors file exists for your specified disease gene set.\nCannot run leave-one-out with PageRank on this disease gene set.", targetName)
    sys.exit(0)


def main():

    algorithm = sys.argv[1]
    pathToPPINetworkFile = sys.argv[2]
    pathToDiseaseGeneFile = sys.argv[3]
    param = float(sys.argv[4])
    outputFile = sys.argv[5]

    print("loading data from files..")
    ppiGraph = compute_if_not_cached(load_graph, pathToPPINetworkFile, fileName=pathToPPINetworkFile)
    run_all = False
    if algorithm == "Algorithms/DiffusionKernel.py":
        function = dk.diffusion_kernel
    elif algorithm == "Algorithms/PageRank.py":
        function = pr.page_rank
    elif algorithm == "Algorithms/RandomWalk.py":
        function = rwr.random_walk
    elif algorithm == "All":
        run_all = True
    else:
        function = None
    
    if run_all:
        print("Starting leave one out for all three algorithms")
        leave_one_out(rwr.random_walk, pathToDiseaseGeneFile, ppiGraph, param)
        leave_one_out(pr.page_rank, pathToDiseaseGeneFile, ppiGraph, param)
        leave_one_out(dk.diffusion_kernel, pathToDiseaseGeneFile, ppiGraph, param)
        print("Finished leave one out for all three algorithms")
    else:
        result = leave_one_out(function, pathToDiseaseGeneFile, ppiGraph, param)


    print("Saving results to:", outputFile)
    with open(outputFile, "w") as of:
        of.write("Leave-one-out Validation Results:\n\nAlgorithm:\t\t{0}\nPPI Graph:\t\t{1}\nDisease Genes:\t\t{2}\nPercentage Correctly Found Genes:\t\t{3}%\n\n".format(algorithm, pathToPPINetworkFile, pathToDiseaseGeneFile, result*100))


if __name__ == '__main__':
    main()
