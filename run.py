# This script is intended to be the one-stop-shop for running any of our algorithms/tests, on any valid dataset.
# It presents the user with a clean command line menu interface for selecting the algorithm/program they want to run, as well as the datasets to use.


from signal import signal, SIGINT
import sys
import os
from termcolor import colored, cprint


# Interface Utility Functions

def resetScreen():
    os.system("clear")
    # Print some nice/standard header here so users stay oriented within our project?


def sigint_handler(signalReceived, frame):
    # Handle any cleanup here
    print('\nExiting..\n')
    exit(0)



# Selection/Execution Functions

def get_files_in_directory(path):
    return [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]

def get_ppi_data_files():
    ppiDataFiles = []
    for f in get_files_in_directory("Data/"):
        if 'ppi' in f.split('.'):
            ppiDataFiles.append(f)
    return ppiDataFiles

def get_disease_gene_files():
    diseaseGeneFiles = []
    for f in get_files_in_directory("Data/"):
        if 'diseasegenes' in f.split('.'):
            diseaseGeneFiles.append(f)
    return diseaseGeneFiles




def select_dataset():
    resetScreen()
    print("\n\nSelect the PPI network dataset you'd like to analyze:\n\n")
    
    datasets = {}
    for i, f in enumerate(get_ppi_data_files(), start=1):
        datasets[i] = f
        print(("\t- " + colored("{0}", "cyan") + ": {1}").format(i, f))

    print("\n\n\tNot seeing your data file? Make sure it is in the Data/ directory and has '.ppi' somewhere in its name.\n\n")

    choice = 0
    while choice == 0:
        try:
            choice = int(input("Select a dataset: >>"))
        except ValueError:
            cprint("please enter a number", "red")

    return datasets[choice]


def select_disease_gene_file():
    resetScreen()
    print("\n\nSelect the disease gene file you'd like to use:\n\n")

    diseaseGeneFiles = {}
    for i, f in enumerate(get_disease_gene_files(), start=1):
        diseaseGeneFiles[i] = f
        print(("\t- " + colored("{0}", "cyan") + ": {1}").format(i, f))


    print("\n\n\tNot seeing your data file? Make sure it is in the Data/ directory and has '.diseasegenes' somewhere in its name.\n\n")

    choice = 0
    while choice == 0:
        try:
            choice = int(input("Select a disease gene file: >>"))
        except ValueError:
            cprint("please enter a number", "red")

    return diseaseGeneFiles[choice]



def select_program():
    resetScreen()
    print("\n\nSelect the program you'd like to run:\n\n")

    cprint("---ALGORITHMS---\n", "green")
    print("\t- " + colored("1", "cyan") + ": Diffusion kernel")
    print("\t- " + colored("2", "cyan") + ": PageRank")
    print("\t- " + colored("3", "cyan") + ": Random walk with restart")
    cprint("\n---VALIDATION---\n", "green")
    print("\t- " + colored("4", "cyan") + ": Area under ROC curve")
    print("\t- " + colored("5", "cyan") + ": Leave one out cross validation")

    programs = [
        "" #zero index
        "DiffusionKernel/DiffusionKernel.py"
        "PageRank/PageRank.py"
        "RWR/Randomwalk.py"
        "Validation/AreaUnderROC.py"
        "LeaveOneOut.py"
    ]

    choice = 0
    while choice == 0:
        try:
            choice = int(input("Select a program: >>"))
        except ValueError:
            cprint("please enter a number", "red")
    print(len(programs))
    print(choice)

    return programs[choice]




def main():
    #Initialization
    signal(SIGINT, sigint_handler)
    resetScreen()
    print("Welcome to the script")
    input("Press enter to continue: >>")

    # Get user selections for what they want to run

    ppiDataset = select_dataset()

    

    diseaseGeneFile = select_disease_gene_file()

    program = select_program()


    # Run stuff using user selections

    print((colored("\nRunning ", "darkred") + "{0}" + colored("\n  on dataset ", "darkred") + "{1}," + colored("\n  using disease genes ", "darkred") + "{2}").format(program, ppiDataset, diseaseGeneFile))





if __name__ == "__main__":
    main()
