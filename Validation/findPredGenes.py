import sys
sys.path.insert(1, '../Imports/')
# Insert relative paths for calls from run.py
sys.path.insert(1, 'Imports/')
import StringNameConverter as snc

def find_pred_genes(outputFile, diseaseProteins):
    outputProts= []
    diseaseProts= []
    with open(outputFile, 'r') as input_file:
        for line in input_file:
            line = line.strip().split(',')
            outputProts.append(line[0].replace('"', ''))
        input_file.close()

    with open(diseaseProteins, 'r') as input_file2:
        for line in input_file2.readlines():
            line = line.strip('\n')
            diseaseProts.append(line)

    predictedProts = []
    threshold = 0
    for i in range(len(outputProts)):
        for j in range(len(diseaseProts)):
            if outputProts[i] != diseaseProts[j] and threshold<5:
                predictedProts.append(outputProts[i])
                threshold+=1
            break

    output = []
    table = snc.load_lookup_table()
    for item in predictedProts:
        output.append(snc.string_to_name(table, item))

    return output

def main():
    outputFile = 'Validation/ischaemic_dk_results.tsv'
    diseaseFile = 'Data/ischaemic-proteins.diseasegenes.tsv'
    print(find_pred_genes(outputFile, diseaseFile))

if __name__ == '__main__':
    main()
