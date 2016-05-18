import sys
import random
n = int(sys.argv[1])
#Generate haplotypes
h0, h1 = "", ""
for _ in range(n):
    hap0 = random.randint(0,1)
    hap1 = 1 - hap0
    h0 += str(hap0)
    h1 += str(hap1)
output = open(str(sys.argv[2]) + "_haplotypes.txt", 'w')
output.write(h0 + '\n')
output.write(h1 + '\n')
output.close()
#Produce sequence reads
AVGREADS = 15
READDEV = 4
MINREADS = 1
AVGLEN = 4
LENDEV = 2
MINLEN = 2
start = 0
reads = open(str(sys.argv[2]) + ".txt", 'w')
while start < n:
    numReads = max(int(round(random.gauss(AVGREADS, READDEV),0)), MINREADS)
    for _ in range(numReads):
        readLen = int(max(round(random.gauss(AVGLEN, LENDEV),0), MINLEN))
        currHap = random.randint(0,1)
        reads.write('-'*start)
        for pos in range(start, min(start+readLen,n)):
            if currHap:
                reads.write(h1[pos])
            else:
                reads.write(h0[pos])
        for i in range(min(start+readLen,n), n):
            reads.write('-')
        reads.write('\n')
    start += 1
reads.close()
