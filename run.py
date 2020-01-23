# This script is intended to be the one-stop-shop for running any of our algorithms/tests, on any valid dataset.
# It presents the user with a clean command line menu interface for selecting the algorithm/program they want to run, as well as the datasets to use.


from signal import signal, SIGINT
import sys
import os
from termcolor import colored, cprint



def resetScreen():
    os.system("clear")
    # Print some nice/standard header here so users stay oriented within our project?


def sigint_handler(signalReceived, frame):
    # Handle any cleanup here
    print('\nExiting..\n')
    exit(0)


def get_files_in_directory(path):
    return [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]

def get_ppi_data_files():
    ppiDataFiles = []
    for f in get_files_in_directory("Data/"):
        if 'ppi' in f.split('.'):
            ppiDataFiles.append(f)
    return ppiDataFiles




def select_dataset():
    resetScreen()
    print("\n\nSelect the PPI network dataset you'd like to analyze:\n\n")
    
    responses = {}
    for i, f in enumerate(get_ppi_data_files(), start=1):
        responses[i] = f
        print(("\t- " + colored("{0}", "cyan") + ": {1}").format(i, f))

    print("\n\n\tNot seeing your data file? Make sure it is in the Data/ directory and has '.ppi' somewhere in its name.\n\n")

    choice = 0
    while choice == 0:
        try:
            choice = int(input("Select a dataset: >>"))
        except ValueError:
            cprint("please enter a number", "red")



def select_disease_gene_file():
    pass


def select_program():
    pass 




def main():
    #Initialization
    signal(SIGINT, sigint_handler)
    resetScreen()
    print("Welcome to the script")
    input("Press enter to continue: >>")

    # Get user selections for what they want to run

    ppiDataset = select_dataset()

    #diseaseGeneFile = select_disease_gene_file()

    #program = select_program()


    # Run stuff using user selections





if __name__ == "__main__":
    main()
