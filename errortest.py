#Given correct haplotypes, approximates the sequencing error rate from a read matrix. 
haps = open(raw_input("Haplotypes: "), 'r')
h0 = haps.readline()
h1 = haps.readline()
haps.close()
matrix = open(raw_input("Reads: "), 'r')
start,tot,err=0,0,0.0
numlines=0
for read in matrix:
    if (start >= len(read)):
        break
    numlines+=1
    if (read[start] == '-'):
        start+=1
    pos = start
    err0, err1 = 0, 0
    while (pos < len(read) and read[pos] != '-'):
        if read[pos] != h0[pos]:
            err0+=1
        elif read[pos] != h1[pos]:
            err1+=1
        tot+=1
        pos+=1
    err += float(min(err0,err1))
print "Error rate: ", err/tot
        

