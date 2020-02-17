import matplotlib.pyplot as plt
import numpy as np
import random


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
		#pred_in, pred_out = unpack_output(file)
	pred_in = [5, 4, 5, 6, 2, 3]
	pred_out = [10, 11, 29, 2]
	fig1, ax = plt.subplots()
	ax.boxplot([pred_in, pred_out])
	rand_x_in = []
	for i in range(len(pred_in)):
		rand_x_in.append(random.uniform(0.95, 1.05))
	ax.scatter(rand_x_in, pred_in, marker = ".", alpha = 0.7)
	rand_x_out = []
	for i in range(len(pred_out)):
		rand_x_out.append(random.uniform(1.95, 2.05))
	ax.scatter(rand_x_out, pred_out, marker = ".", alpha = 0.7)
	plt.ylabel("Degree")
	plt.xticks([1, 2], ("In", "Out"))
	plt.show()
	return 0
main()
