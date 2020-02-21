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
from leaveOneOut import leave_one_out

path_to_ppi = sys.argv[1]
all_dg_file_paths = ['Data/endometriosis-proteins.diseasegenes.tsv', 'Data/lymphoma-proteins.diseasegenes.tsv', 'Data/ischaemic-proteins.diseasegenes.tsv']
params = [0.2, 0.4, 0.6, 0.8, 1.0]
print("loading data from files..")
ppiGraph = compute_if_not_cached(load_graph, path_to_ppi, fileName=path_to_ppi)

for i in range(3):
	print("----------------------- DISEASE GENE FILE IS:", all_dg_file_paths[i], "---------------------------------")
	result_rwr = leave_one_out(rwr.random_walk, all_dg_file_paths[i], ppiGraph, 0.4)
	print("PERCENTAGE OF GENES FOUND FOR RANDOM WALK: ", result_rwr)
	result_pr = leave_one_out(pr.random_walk, all_dg_file_paths[i], ppiGraph, 0.4)
	print("PERCENTAGE OF GENES FOUND FOR PAGERANK: ", result_pr)
	result_dk = leave_one_out(dk.diffusion_kernel, all_dg_file_paths[i], ppiGraph, 0.4)
	print("PERCENTAGE OF GENES FOUND FOR DIFFUSION KERNEL:", result_dk)