"""
This script was used to create the test data file by taking the STRING data for a specific organism and eliminating a few of the proteins to cut down the size of the network.
"""




def create_file(path):
    inputFile = open(path, 'r')
    outputFile = open("test-graph-data.tsv", 'w')
    for line in inputFile:
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
            outputFile.write(line)

    inputFile.close()
    outputFile.close()





def main():
    create_file("test-graph-data-old.tsv")



if __name__ == '__main__':
    main()
