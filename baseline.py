#Assuming no sequencing errors and each position has been sequenced with at least length 2
reads = raw_input("Reads File Name: ")
matrix = open(reads, 'r')
output = raw_input("Output File Name: ")
h0, h1, pos = "1", "0", 0
for read in matrix:
    if (pos == (len(read) - 2)):
        break
    if (read[pos] == '-' or read[pos+1] == '-'):
        continue
    if (read[pos] == h0[pos]):
        while pos < len(read) - 2 and read[pos+1] != '-':
            h0 += (read[pos+1])
            h1 += (str(1-int(read[pos+1])))
            pos+=1
    else:
        while pos < len(read) -2 and read[pos+1] != '-':
            h1 += (read[pos+1])
            h0 += (str(1-int(read[pos+1])))
            pos+=1
haplotypes = open(output, 'w')
haplotypes.write(h0 + '\n' + h1 + '\n')
haplotypes.close()
matrix.close()
if (raw_input("Check results? y/n: ") == 'y'):
    correct = open(raw_input("Haplotype file: "), 'r').readline()[:-1]
    if correct == h0 or correct == h1:
        print "CORRECT"
    else:
        print "WRONG"
