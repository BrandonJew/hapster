#Generates test data for haplotype assembly. Be sure to run as "python datagen.py output" to specify output file (Default is output.txt)
import sys
import random
random.seed()
if len(sys.argv) == 2:
    out = str(sys.argv[1])
else:
    out = 'output'
n = int(raw_input("Number of SNPs: ")) #n = int(sys.argv[1])
errorrate = float(raw_input("Enter error rate:" )) #errorrate = float(sys.argv[3])
#Generate haplotypes
h0, h1 = "", ""
for _ in range(n):
    hap0 = random.randint(0,1)
    hap1 = 1 - hap0
    h0 += str(hap0)
    h1 += str(hap1)
output = open(out + "_haplotypes.txt", 'w')
output.write(h0 + '\n')
output.write(h1 + '\n')
output.close()
#Produce sequence reads
AVGREADS = 10
READDEV = 2
MINREADS = 7
AVGLEN = 5
LENDEV = 0.5
#MINLEN = 2
start = 0
err, tot = 0,0
reads = open(out + ".txt", 'w')
while start < n:
    numReads = max(int(round(random.gauss(AVGREADS, READDEV),0)), MINREADS)
    readLen = int(max(round(random.gauss(AVGLEN, LENDEV),0), AVGLEN))
    for _ in range(numReads):
        currHap = random.randint(0,1)
        reads.write('-'*start)
        for pos in range(start, min(start+readLen,n)):
            if currHap:
                if random.uniform(0,1) < errorrate:
                    err+=1
                    reads.write(str(1-int(h1[pos])))
                else:
                    reads.write(h1[pos])
            else:
                if random.uniform(0,1) < errorrate:
                    err+=1
                    reads.write(str(1-int(h0[pos])))
                else:
                    reads.write(h0[pos])
            tot+=1
        for i in range(min(start+readLen,n), n):
            reads.write('-')
        reads.write('\n')
    start += 1
reads.close()
print "Observed Error Rate: ", float(err)/tot
