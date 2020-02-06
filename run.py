# This script is intended to be the one-stop-shop for running any of our algorithms/tests, on any valid dataset.
# It presents the user with a clean command line menu interface for selecting the algorithm/validation method they want to run, as well as the datasets to use.


from signal import signal, SIGINT
import subprocess
import sys
import os
try:
    from termcolor import colored, cprint
    termcolorMissing = False
except ImportError:
    termcolorMissing = True
    def cprint(m, c):
        print(m)
    def colored(m, c):
        return m



# Interface Utility Functions

def resetScreen():
    os.system("clear")
    # Print some nice/standard header here so users stay oriented within our project?


def sigint_handler(signalReceived, frame):
    # Handle any cleanup here
    print('\nExiting..\n')
    exit(0)

def checkDependencies():
    print("\nChecking for required libraries...")

    if termcolorMissing:
        print("Module termcolor not found. Install termcolor python module for color.")

    error = False
    for lib in ['networkx', 'numpy', 'scipy', 'matplotlib', 'requests', 'tkinter']:
        error = checkPipLibrary(lib) or error
    
    if error:
        choice = input("\nWould you like to try to automatically install missing dependencies? (Y/n) >>")
        if choice.upper() != "Y":
            print("\nExiting..\n")
            sys.exit(0)
        else:
            print("\nInstalling dependencies..\n")
            subprocess.Popen(["python3", "setup.py"]).wait()
            print("\nExiting script.. please restart it.\n")
            sys.exit(0)
            
    else:
        cprint("Success.", "green")

    print("\nChecking for datasets..\n")

    dgf = get_disease_gene_files()
    ppif = get_ppi_data_files()

    if len(ppif) == 0:
        cprint("No PPI network found.", "red")
        choice = input("Would you like to automatically download the newest one from String-DB.org? (Y/n) >>")
        if choice.upper() != "Y":
            print("\nExiting..\n")
            sys.exit(0)
        else:
            print("Installing String dataset..")
            p = subprocess.Popen(["bash", "download-human-dataset.sh"])
            p.wait()
            if p.returncode != 0:
                cprint("Something went wrong. Exiting.", "red")
                sys.exit(1)
            else:
                print("\nString data successfully installed.\n")
    
    if len(dgf) == 0:
        print("\n There are no disease-gene files installed. Disease gene files must be in the Data/ directory and must include '.diseasegenes' in the file name.")
        print("Exiting..")
        sys.exit(0) 

    cprint("Success.", "green")   


def checkPipLibrary(lib):
    print("checking for {0}".format(lib), end=" ")
    error = False
    try:
        __import__(lib)
    except ImportError:
        print(colored("missing", "red"))
        error = True
    if not error:
        print(colored("yes", "green"))
    return error 





# Selection/Execution Functions

def get_files_in_directory(path):
    return [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]

def get_ppi_data_files():
    ppiDataFiles = []
    for f in get_files_in_directory("Data/"):
        if 'ppi' in f.split('.'):
            ppiDataFiles.append("Data/" + f)
    return ppiDataFiles

def get_disease_gene_files():
    diseaseGeneFiles = []
    for f in get_files_in_directory("Data/"):
        if 'diseasegenes' in f.split('.'):
            diseaseGeneFiles.append("Data/" + f)
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
    while choice == 0 or choice > len(datasets):
        try:
            choice = int(input("Select a dataset: >>"))
        except ValueError:
            cprint("please enter a number", "red")
        if choice > len(datasets):
            cprint("number must be between 1 and {0}".format(len(datasets)), "red")
            choice = 0

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
    while choice == 0 or choice > len(diseaseGeneFiles):
        try:
            choice = int(input("Select a disease gene file: >>"))
        except ValueError:
            cprint("please enter a number", "red")
        if choice > len(diseaseGeneFiles):
            cprint("number must be between 1 and {0}".format(len(diseaseGeneFiles)), "red")
            choice = 0

    return diseaseGeneFiles[choice]



def select_algorithm():
    resetScreen()
    print("\n\nSelect the algorithm you'd like to run:\n\n")

    cprint("---ALGORITHMS---\n", "green")
    print("\t- " + colored("1", "cyan") + ": Diffusion kernel")
    print("\t- " + colored("2", "cyan") + ": PageRank")
    print("\t- " + colored("3", "cyan") + ": Random walk with restart")

    print("\n\n")
    

    algorithms = {
        1:"Algorithms/DiffusionKernel.py",
        2:"Algorithms/PageRank.py",
        3:"Algorithms/RandomWalk.py"
    }

    choice = 0
    while choice == 0 or choice > len(algorithms):
        try:
            choice = int(input("Select an algorithm: >>"))
        except ValueError:
            cprint("please enter a number", "red")
        if choice > len(algorithms):
            cprint("number must be between 1 and {0}".format(len(algorithms)), "red")
            choice = 0

    return algorithms[choice]


def select_beta_value():
    resetScreen()
    print("\nDiffusion kernel allows you to specify a beta value that _____")
    print("\nPlease enter a beta value between " + colored("0", "cyan") + " and " + colored("2", "cyan") + ".")

    choice = float("inf")
    while choice < 0 or choice > 2:
        try:
            choice = float(input("Select beta: >>"))
        except ValueError:
            cprint("please enter a decimal number between 0 and 2", "red")
            continue
        if choice < 0 or choice > 2:
            cprint("please enter a decimal number between 0 and 2", "red")
            choice = float("inf")
    return choice


def select_pr_beta_value():
    resetScreen()
    print("\nPageRank allows you to specify a beta value that sets the probability of restarting from a known disease gene.\nA value of 0 means the algorithm will never 'restart', while a value of 1 means that the algorithm will only ever visit known disease genes. (always restart)\nMost researchers set a relatively low beta value, around .2")
    print("\nPlease enter a beta value between " + colored("0", "cyan") + " and " + colored("1", "cyan") + ".")

    choice = float("inf")
    while choice < 0 or choice > 1:
        try:
            choice = float(input("Select beta: >>"))
        except ValueError:
            cprint("please enter a decimal number between 0 and 1", "red")
            continue
        if choice < 0 or choice > 1:
            cprint("please enter a decimal number between 0 and 1", "red")
            choice = float("inf")
    return choice


def select_rwr_r_value():
    resetScreen()
    print("\nRandom Walk with Restart allows you to specify an R value that sets the probability of restarting from a known disease gene.\nA value of 0 means the algorithm will never 'restart', while a value of 1 means that the algorithm will only ever visit known disease genes. (always restart)\nMost researchers set a relatively low R value, around .2")
    print("\nPlease enter an R value between " + colored("0", "cyan") + " and " + colored("1", "cyan") + ".")

    choice = float("inf")
    while choice < 0 or choice > 1:
        try:
            choice = float(input("Select R: >>"))
        except ValueError:
            cprint("please enter a decimal number between 0 and 1", "red")
            continue
        if choice < 0 or choice > 1:
            cprint("please enter a decimal number between 0 and 1", "red")
            choice = float("inf")
    return choice


def select_validation():
    resetScreen()
    print("\n\nSelect the validation method you'd like to use:\n\n")

    cprint("---VALIDATION---\n", "green")
    print("\t- " + colored("1", "cyan") + ": None")
    print("\t- " + colored("2", "cyan") + ": Area under ROC curve")
    print("\t- " + colored("3", "cyan") + ": Leave one out cross validation")

    print("\n\n")

    validations = {
        1:"None",
        2:"Validation/AreaUnderROC.py",
        3:"Validation/LeaveOneOut.py"
    }

    choice = 0
    while choice == 0 or choice > len(validations):
        try:
            choice = int(input("Select a validation method: >>"))
        except ValueError:
            cprint("please enter a number", "red")
        if choice > len(validations):
            cprint("number must be between 1 and {0}".format(len(validations)), "red")
            choice = 0

    return validations[choice]





def main():
    #Initialization
    signal(SIGINT, sigint_handler)
    resetScreen()

    #Check for dependencies
    print("Welcome to the script")
    checkDependencies()

    # Get user selections for what they want to run

    ppiDataset = select_dataset()

    diseaseGeneFile = select_disease_gene_file()

    algorithm = select_algorithm()

    if algorithm == "Algorithms/DiffusionKernel.py":
        beta = select_beta_value()

    if algorithm == "Algorithms/PageRank.py":
        beta = select_pr_beta_value()

    if algorithm == "Algorithms/RandomWalk.py":
        R = select_rwr_r_value()



    validation = select_validation()


    # Confirm user selections
    resetScreen()
    print((colored("\nRunning:\t\t", "yellow") + "{0}" + colored("\n  on dataset:\t\t", "yellow") + "{1}" + colored("\n  using disease genes:\t", "yellow") + "{2}" + colored("\n\nValidating with:\t", "yellow") + "{3}").format(algorithm, ppiDataset, diseaseGeneFile, validation))
    input(colored("\nPress enter to continue (ctrl+c to cancel)..", "green"))

    # Run stuff
    if algorithm == "Algorithms/DiffusionKernel.py" or algorithm == "Algorithms/PageRank.py":
        numeric = beta
    else:
        numeric = R

    if validation == "None":
        cmd = "python3 {0} {1} {2} {3}".format(algorithm, ppiDataset, diseaseGeneFile, numeric)
        os.system(cmd)
    else:
        cmd = "python3 {0} {1} {2} {3} {4}".format(validation, algorithm, ppiDataset, diseaseGeneFile, numeric)
        os.system(cmd)



if __name__ == "__main__":
    main()
