import sys

ERROR_MESSAGE = "ERROR. File path needs to be specified in the command line."

"""
This script was used to create the test data file by taking the STRING data 
0for a specific organism and eliminating a few of the proteins to cut down 
the size of the network.
"""
def create_file(path):
    input_file = open(path, 'r')
    output_file = open(path + ".new", 'w+')
    for line in input_file:
        if line.strip() == "#":
            break
        data = line.strip().split(" ")
        if not (data[0] == "568816.Acin_0001" or data[1] == "568816.Acin_0001" or
                data[0] == "568816.Acin_0002" or data[1] == "568816.Acin_0002" or
                data[0] == "568816.Acin_0003" or data[1] == "568816.Acin_0003" or
                data[0] == "568816.Acin_0004" or data[1] == "568816.Acin_0004" or
                data[0] == "568816.Acin_0005" or data[1] == "568816.Acin_0005" or
                data[0] == "568816.Acin_0006" or data[1] == "568816.Acin_0006" or
                data[0] == "568816.Acin_0007" or data[1] == "568816.Acin_0007" or
                data[0] == "568816.Acin_0008" or data[1] == "568816.Acin_0008" or
                data[0] == "568816.Acin_0009" or data[1] == "568816.Acin_0009" or
                data[0] == "568816.Acin_0010" or data[1] == "568816.Acin_0010"):
            output_file.write(line)
    input_file.close()
    output_file.close()


def main():
    if len(sys.argv) != 2:
        print(ERROR_MESSAGE)
        sys.exit()
    create_file(sys.argv[1])


if __name__ == '__main__':
    main()
