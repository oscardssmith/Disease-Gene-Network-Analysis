import matplotlib.pyplot as plt
import numpy as np 


def unpack_output(file):
	predicted_in = []
	predicted_out = []
	with open(file, 'r') as input_file:
		for line in input_file.readlines():
			line = line.strip('\n').split('\t')
			if int(line[2]) == -1:
				predicted_out.append(line[1])
			else:
				predicted_in.append(line[1])
	return predicted_in, predicted_out

def main():
	# for each file
		#in, out = unpack_output(file)
		#fig, ax = plt.subplots()
		#ax.boxplot([in, out])
	pred_in = [5, 4, 5, 6, 2, 3]
	pred_out = [10, 11, 29, 2]
	fig, ax = plt.subplots()
	ax.boxplot([pred_in, pred_out])
	plt.show()
	return 0