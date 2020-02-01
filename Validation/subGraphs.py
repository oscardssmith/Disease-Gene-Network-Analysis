'''
Generating subgraphs for cytoscape
'''
import sys
sys.path.insert(1, '../Algorithms/')
sys.path.insert(1, 'Algorithms/')
sys.path.insert(1, '../Scripts/')
import loader

import networkx as nx

def main():
    PPI_Network = loader.load_graph('../Data/9606.protein.links.v11.0.txt')
    #ground_truth_files = ['../Data/MalaCard-protein-Endometriosis.diseasegenes.tsv', '../Data/MalaCard-protein-ischaemic-stroke.diseasegenes.tsv','../Data/MalaCard-protein-lymphoma.diseasegenes.tsv']
    file_paths = ['../Data/endometriosis-proteins.diseasegenes.tsv','../Data/lymphoma-proteins.diseasegenes.tsv', '../Data/ischaemic-stroke-proteins.diseasegenes.tsv']
    names = ["endometriosis.graphml", "lymphoma.graphml", "ischaemic_stroke.graphml"]
    for i in range(3):
        nodes_in_subgraph =[]
        with open(file_paths[i], 'r') as input:
            input = input.readlines()
            for line in input:
                protein = line.rstrip('\n')
                nodes_in_subgraph.append(protein)

                neighbors = list(PPI_Network.neighbors(protein))
                nodes_in_subgraph = nodes_in_subgraph + neighbors

        subgraph = PPI_Network.subgraph(nodes_in_subgraph)

        nx.write_graphml(subgraph, names[i])
        print("exported file named" + names[i])



if __name__ == '__main__':
    main()
