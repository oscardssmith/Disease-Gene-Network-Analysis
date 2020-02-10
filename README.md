# CompBioComps
Place for code for 2019 Comps

Disease Gene Prioritization with Algorithms

## PROJECT DESCRIPTION
 A gene is a sequence of DNA that codes for a protein that is associated with a particular biological function (you can think of genes and proteins as having a 1-1 relationship). A mutation, or aberration, that occurs in a gene and disrupts the normal function is called a disease gene.

 In our project, we run three graph algorithms: page rank, random walk with restart, and diffusion kernel. Each algorithm takes a set of known disease genes and uses this information to predict which other genes are associated with a given disease.

## GETTING STARTED
 In order to use our software package, first clone the GitHub repository onto your development environment. Open your terminal/command prompt and use the following command

 ```bash
 git clone https://github.com/oscardssmith/CompBioComps.git
 ```
Once cloned, use the command

```bash
cd CompBioComps
```
to make our repository your working directory. Now use the command
```bash
python3 run.py
```
to complete independencies, download appropriate files, and pull up an interface that allows you to interact with our code.

## INTERFACE
The interface will first ask for the task you would like to perform.

If you choose algorithm, you will be asked to choose which algorithm you would like to run (page rank, diffusion kernel, random walk with restart). Depending on the algorithm, you may also be asked to input a beta or R value. Suggested ranges exist within the interface. You will then be asked to select the PPI network you want to analyze. There will be one dataset by default unless you add a dataset of your own. You will then be asked to select the gene file you want to use. We have 9 gene files setup and ready to use. Finally, you will be asked to name your output file. More details on output are below.

If you choose validation, you will be asked which type of validation (ROC or leave one out). Once chosen, you will then go through the same process as above in terms of selecting an algorithm, parameters for the algorithm, dataset, gene file, and output file.

## OUTPUT
You will be asked to specify the name of your output file. This will appear in the results folder. There are three cases for the output files:
1. Running Algorithms- output file will be a csv file of protein names with probability, in descending order.
2. Area under the ROC curve- A picture of the ROC curve as a .png image
3. Leave one out- _________FILL ME IN__________________
