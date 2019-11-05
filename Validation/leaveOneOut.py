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




def leaveOneOut(function, diseaseGeneFilePath, PPI_Network):
    diseaseGeneFile = open(diseaseGeneFilePath, 'r')
    allDiseaseGenes = diseaseGeneFile.read().splitlines()
    diseaseGeneFile.close()

    squaredDifferenceSum = 0

    for skipGene in allDiseaseGenes:
        #create leave-one-out disease gene list
        diseaseGeneList = []
        for gene in allDiseaseGenes:
            if gene != skipGene:
                diseaseGeneList.append(gene)

        #run algorithm using modified disease gene file
        output = function(diseaseGeneList, PPI_Network)

        #find the predicted probability of the omitted gene and add it to the current sum
        for pair in output:
            if pair[0] == skipGene:
                predictedProbability = pair[1]
                squaredDifferenceSum += (1 - predictedProbability)**2
                break


    #Find average of all squared differences
    meanSquaredDifference = squaredDifferenceSum/len(allDiseaseGenes)
    return meanSquaredDifference




def load_PPI_Network(filePath):
    return loader.load_graph(filePath)



def main():
    # load the full PPI Network
    PPI_Network = load_PPI_Network('../Data/test-graph-data.tsv')

    result = leaveOneOut(rwr.RandomWalk, '../Data/EndometriosisProteins.tsv', PPI_Network)
    print("Mean squared difference for RWR:", result)

    result = leaveOneOut(dk.DiffusionKernel, 'diseaseGeneFile', PPI_Network)
    print("Mean squared difference for Diffusion Kernel:", result)

    result = leaveOneOut(pr.PageRank, 'diseaseGeneFile', PPI_Network)
    print("Mean squared difference for PageRank:", result)




if __name__ == '__main__':
    main()
