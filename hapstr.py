#Assuming position has been sequenced with at least length 2 and correct haplotype is sequenced most frequently. 
import operator
import array
reads = raw_input("Reads File Name: ")
matrix = open(reads, 'r')
output = raw_input("Output File Name: ")
h0, h1, start, backtrack = array.array('c',"1"), array.array('c',"0"), 0, 0
opts = {'10': 0, '11': 0, '00' : 0, '01' : 0} #Opts contains possibilities for h0 at start
begin,first = True, True
currline = 0
n = len(matrix.readline())
matrix.seek(0)
while True:
    while True: #If backtrack goes too far, goes forward until usable read is reached
        matrix.seek((currline*n)+start)
        read = matrix.read(2)
        if (begin and (read[0] == '-' or read[1] == '-')):
            matrix.seek((n-2),1) #Goes forward one line
            currline+=1
            continue
        break
    if (start >= (n - 2)):
        break
    if (not begin) and (read[0] == '-' or read[1] == '-'):
        best = max(opts.iteritems(), key=operator.itemgetter(1))[0] #Retrieves 2mer with highest hits
        h0[start] = best[0]
        h1[start] = str(1-int(best[0]))
        h0.append(best[1])
        h1.append(str(1-int(best[1])))
        start += 1
        opts = {'10': 0, '11': 0, '00' : 0, '01' : 0}
        currline = max((currline-backtrack), 0) #Moving to new start position, so go back to earliest read that contain 2mers at this position
        begin,first = True, False
        continue
    if (read[0] == h0[start]):
        opts[read] += 1
    else:
        opts[('0' + str(abs(int(read) - 11)))[-2:]] += 1 #Inverts 2mer 
    if first: #Uses the number and length of reads at position 0 as the amount to backtrack when moving to a new start position
    	matrix.seek(currline*n)
        while (matrix.read(1) != '-'):
            backtrack += 1
    currline+=1
    begin = False
matrix.close()
haplotypes = open(output, 'w')
haplotypes.write(h0.tostring() + '\n' + h1.tostring() + '\n')
haplotypes.close()
if (raw_input("Check results? y/n: ") == 'y'):
    correct = open(raw_input("Haplotype file: "), 'r').readline()[:-1]
    if correct == h0.tostring() or correct == h1.tostring():
        print "CORRECT"
    else:
        print "WRONG"
