import array
import operator
import random

#returns proportion of correct SNPs assembled
def correctness(output):
	correct = open("output_haplotypes.txt", 'r').readline()
	tot,err0,err1 = len(correct)-1, 0,0
	for pos in range(0,tot):
		if correct[pos] != output[pos]:
		    err0 += 1
		if correct[pos] != str(1-int(output[pos])):
		    err1 += 1
	err = min(err0,err1)
	return 1 - (float(err)/tot)

def generatedata(n, errorrate):
	random.seed()
	out = 'output'
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
		                reads.write(str(1-int(h1[pos])))
		            else:
		                reads.write(h1[pos])
		        else:
		            if random.uniform(0,1) < errorrate:
		                reads.write(str(1-int(h0[pos])))
		            else:
		                reads.write(h0[pos])
		    for i in range(min(start+readLen,n), n):
		        reads.write('-')
		    reads.write('\n')
		start += 1
	reads.close()
	
def baseline():
	matrix = open("output.txt", 'r')
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
	matrix.close()
	return h0

def hapster():
	matrix = open("output.txt", 'r')
	h0, h1, start, backtrack = array.array('c',"1"), array.array('c',"0"), 0, 0
	opts = {'10': 0, '11': 0, '00' : 0, '01' : 0} #Opts contains possibilities for h0 at start
	begin,first = True, True
	while True:
		while True: #If backtrack goes too far, goes forward until usable read is reached
		    read = matrix.readline()
		    if (begin and (start < len(read) - 1) and (read[start] == '-' or read[start+1] == '-')):
		        continue
		    break
		if (start >= (len(read) - 2)):
		    break
		if (not begin) and (read[start] == '-' or read[start+1] == '-'):
		    best = max(opts.iteritems(), key=operator.itemgetter(1))[0] #Retrieves 2mer with highest hits
		    h0[start] = best[0]
		    h1[start] = str(1-int(best[0]))
		    h0.append(best[1])
		    h1.append(str(1-int(best[1])))
		    start += 1
		    opts = {'10': 0, '11': 0, '00' : 0, '01' : 0}
		    matrix.seek(max((matrix.tell())-(len(read)*backtrack), 0)) #Moving to new start position, so go back to earliest read that contain 2mers at this position
		    begin,first = True, False
		    continue
		if (read[start] == h0[start]):
		    opts[read[start:start+2]] += 1
		else:
		    opts[('0' + str(abs(int(read[start:start+2]) - 11)))[-2:]] += 1 #Inverts 2mer 
		if first: #Uses the number and length of reads at position 0 as the amount to backtrack when moving to a new start position
			pos = start
			while (read[pos] != '-'):
				backtrack += 1
				pos += 1
		begin = False
	matrix.close()
	return h0.tostring()
