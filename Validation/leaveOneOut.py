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
sys.path.insert(1, '../RWR/')
sys.path.insert(1, '../DiffusionKernel/')
sys.path.insert(1, '../PageRank/')
sys.path.insert(1, '../Scripts/')
import RandomWalk as rwr
import DiffusionKernel as dk
import PageRank as pr
import loader
import time




def leaveOneOut(function, diseaseGeneFilePath, PPI_Network):
    print("Starting leaveOneOut function")

    diseaseGeneFile = open(diseaseGeneFilePath, 'r')
    allDiseaseGenes = diseaseGeneFile.read().splitlines()
    diseaseGeneFile.close()

    numDiseaseGenes = len(allDiseaseGenes)
    rankThreshhold = round(1.2 * numDiseaseGenes)

    numGenesNotFound = 0

    print("finished initialization, starting disease gene loop")

    for index, skipGene in enumerate(allDiseaseGenes):
        print("looping! skipping gene: ", skipGene)

        #create leave-one-out disease gene list
        # diseaseGeneList = []
        # for gene in allDiseaseGenes:
        #     if gene != skipGene:
        #         diseaseGeneList.append(gene)

        startVector = loader.load_start_vector(diseaseGeneFilePath, PPI_Network)
        startVector[index] = 0



        #run algorithm using modified disease gene file
        print("calling algorithm")
        startTime = time.time()
        output = function(PPI_Network, startVector)
        endTime = time.time()
        print("finished algorithm. Time elapsed:", endTime - startTime)

        #find the predicted probability of the omitted gene and add it to the current sum
        print("finding skipgene predicted probability")
        startTime = time.time()
        foundGene = False
        for i in range(rankThreshhold):
            if output[i][0] == skipGene:
                foundGene = True
                break
        if not foundGene:
            numGenesNotFound +=1

        endTime = time.time()
        print("added skipGene predictedProbability to squaredDifferenceSum\nTime elapsed:", endTime - startTime)


    print("------------------------\nFinished running algorithm with all disease genes left out\nCalculating mean squared difference")
    #Find average of all squared differences
    percentIncorrectlyRankedGenes = numGenesNotFound/numDiseaseGenes
    return percentIncorrectlyRankedGenes




def load_PPI_Network(filePath):
    return loader.load_graph(filePath)



def main():
    # # load the full PPI Network
    # totalStartTime = time.time()
    # print("Starting leave one out validation.")
    # print("Loading graph")
    # PPI_Network = load_PPI_Network('../Data/9606.protein.links.v11.0.txt')
    # print("Loaded graph")
    #
    # result = leaveOneOut(rwr.RandomWalk, '../Data/EndometriosisProteins.tsv', PPI_Network)
    # print("Mean squared difference for RWR:", result)
    #
    # #result = leaveOneOut(dk.DiffusionKernel, 'diseaseGeneFile', PPI_Network)
    # #print("Mean squared difference for Diffusion Kernel:", result)
    # 
    # result = leaveOneOut(pr.PageRank, 'diseaseGeneFile', PPI_Network)
    # print("Mean squared difference for PageRank:", result)
    #
    # totalEndTime = time.time()
    # print("Finished leave one out validation!\nTotal time in hours:", (totalEndTime - totalStartTime)/3600)

    # testing if RWR and PageRank are the same -- delete later!
    PPI_Network = load_PPI_Network('../Data/9606.protein.links.v11.0.txt') # load network
    startVector = loader.load_start_vector('../Data/EndometriosisProteins.tsv', PPI_Network)
    result_RWR = rwr.RandomWalk(PPI_Network, startVector)
    result_PR = pr.PageRank(PPI_Network, startVector)
    differenceScore =0
    for i in range(len(result_RWR)):
        if result_RWR[i][0] != result_PR[i][0]:
            differenceScore +=1
    print(differenceScore)



if __name__ == '__main__':
    main()
