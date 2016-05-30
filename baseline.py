#Assuming no sequencing errors and each position has been sequenced with at least length 2
reads = raw_input("Reads File Name: ")
matrix = open(reads, 'r')
output = raw_input("Output File Name: ")
h0, h1, pos, currline = "1", "0", 0, 0
n = len(matrix.readline())
matrix.seek(0)
while True:
    matrix.seek((currline*n)+pos)
    read = matrix.read(2)
    if (pos == (n - 2)):
        break
    if (read[0] == '-' or read[1] == '-'):
        currline+=1
        continue
    if (read[0] == h0[pos]):
        hap = read[1]
        while pos < n - 2 and hap != '-':
            h0 += (hap)
            h1 += (str(1-int(hap)))
            pos+=1
            hap = matrix.read(1)
    else:
        hap = read[1]
        while pos < n - 2 and hap != '-':
            h1 += (hap)
            h0 += (str(1-int(hap)))
            pos+=1
            hap = matrix.read(1)
    currline+=1
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
