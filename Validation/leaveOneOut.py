"""
This program assumes the following:
Each algorithm can be called with only two parameters:
    - the path to the TSV file containing the full PPI network
    - the path to the TSV file containing the known disease genes

Each algorithm produces results in a 2 dimensional array where each item is a tuple containing a gene name and it's ranking/probability,
 and the list is in order of probability, high to low, so that the index is equivalent to the ranking:
    - [[geneName, probability], [geneName, probability], [geneName, probability], [geneName, probability]]


"""
import sys
sys.path.insert(1, '../RWR/')
sys.path.insert(1, '../DiffusionKernel/')
sys.path.insert(1, '../PageRank/')
import RandomWalk as rwr
import DiffusionKernel as dk
import PageRank as pr




def leaveOneOut(function, diseaseGeneFilePath, PPINetworkFilePath):
    diseaseGeneFile = open(diseaseGeneFilePath, 'r')
    allDiseaseGenes = diseaseGeneFile.read().splitlines()
    diseaseGeneFile.close()

    squaredDifferenceSum = 0

    for skipGene in allDiseaseGenes:
        #create leave-one-out disease gene file
        tempGeneFile = open('tempGeneFile.tsv', 'w')
        for gene in allDiseaseGenes:
            if gene != skipGene:
                tempGeneFile.write(gene)

        #run algorithm using modified disease gene file
        output = function('tempGeneFile.tsv', PPINetworkFilePath)

        #find the predicted probability of the omitted gene and add it to the current sum
        for pair in output:
            if pair[0] == skipGene:
                predictedProbability = pair[1]
                break
        squaredDifferenceSum += (1 - predictedProbability)**2

    #Find average of all squared differences
    meanSquaredDifference = squaredDifferenceSum/len(allDiseaseGenes)
    return meanSquaredDifference









def main():
    #Set diseaseGeneFile and PPINetworkFile to the paths to the data you want to test

    result = leaveOneOut(rwr.RandomWalk, 'diseaseGeneFile', 'PPINetworkFile')
    print("Mean squared difference for RWR:", result)
    result = leaveOneOut(dk.DiffusionKernel, 'diseaseGeneFile', 'PPINetworkFile')
    print("Mean squared difference for Diffusion Kernel:", result)
    result = leaveOneOut(pr.PageRank, 'diseaseGeneFile', 'PPINetworkFile')
    print("Mean squared difference for PageRank:", result)




if __name__ == '__main__':
    main()
