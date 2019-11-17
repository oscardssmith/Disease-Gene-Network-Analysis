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
#import DiffusionKernel as dk
import PageRank as pr
import loader




def leaveOneOut(function, diseaseGeneFilePath, PPI_Network):
    print("Starting leaveOneOut function")

    diseaseGeneFile = open(diseaseGeneFilePath, 'r')
    allDiseaseGenes = diseaseGeneFile.read().splitlines()
    diseaseGeneFile.close()

    squaredDifferenceSum = 0

    print("finished initialization, starting disease gene loop")

    for skipGene in allDiseaseGenes:
        print("looping! skipping gene: ", skipGene)

        #create leave-one-out disease gene list
        diseaseGeneList = []
        for gene in allDiseaseGenes:
            if gene != skipGene:
                diseaseGeneList.append(gene)

        #run algorithm using modified disease gene file
        print("calling algorithm")
        startTime = time.time()
        output = function(diseaseGeneList, PPI_Network)
        endTime = time.time()
        print("finished algorithm. Time elapsed:", endTime - startTime)

        #find the predicted probability of the omitted gene and add it to the current sum
        print("finding skipgene predicted probability")
        startTime = time.time()
        for pair in output:
            if pair[0] == skipGene:
                predictedProbability = pair[1]
                squaredDifferenceSum += (1 - predictedProbability)**2
                break
        endTime = time.time()
        print("added skipGene predictedProbability to squaredDifferenceSum\nTime elapsed:", endTime - startTime)



    #Find average of all squared differences
    meanSquaredDifference = squaredDifferenceSum/len(allDiseaseGenes)
    return meanSquaredDifference




def load_PPI_Network(filePath):
    return loader.load_graph(filePath)



def main():
    # load the full PPI Network
    PPI_Network = load_PPI_Network('../Data/test-graph-data.tsv')
    print("Loaded graph")

    result = leaveOneOut(rwr.RandomWalk, '../Data/EndometriosisProteins.tsv', PPI_Network)
    print("Mean squared difference for RWR:", result)

    #result = leaveOneOut(dk.DiffusionKernel, 'diseaseGeneFile', PPI_Network)
    #print("Mean squared difference for Diffusion Kernel:", result)

    #result = leaveOneOut(pr.PageRank, 'diseaseGeneFile', PPI_Network)
    #print("Mean squared difference for PageRank:", result)




if __name__ == '__main__':
    main()
