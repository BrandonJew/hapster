import testsrc as ts
import time
#TEST TIME COMPLEXITY
runs = 10
sizes = [100,1000,5000,10000]
output = open("timecomplexitydata.csv",'w')
output.write("snps,runs,baseline,hapster\n")
for size in sizes:
	output.write(str(size)+','+str(runs)+',')
	basetot,hapstertot = 0,0
	for run in range(runs):
		print "Starting size", str(size), "run #" + str(run)
		ts.generatedata(size,0)
		tstart = time.time()
		ts.baseline()
		basetot+=(time.time()-tstart)
		tstart = time.time()
		ts.hapster()
		hapstertot+=(time.time()-tstart)
	baseavg = str(float(basetot)/runs)
	hapsteravg = str(float(hapstertot)/runs)
	output.write(baseavg + ',' + hapsteravg + '\n')
output.close()
#TEST CORRECTNESS
testsize = 1000
errorinc = 0.05
errorrates = []
errorrate = 0
while errorrate < 1:
	errorrates.append(round(errorrate,2))
	errorrate+=errorinc
output = open("correctnessdata.csv",'w')
output.write("error_rate,runs,baseline,hapster\n")
for error in errorrates:
	output.write(str(error)+','+str(runs)+',')
	basetot, hapstertot = 0,0
	for run in range(runs):
		ts.generatedata(testsize,error)
		basetot += ts.correctness(ts.baseline())
		hapstertot += ts.correctness(ts.hapster())
	baseavg = str(float(basetot)/runs)
	hapsteravg = str(float(hapstertot)/runs)
	output.write(baseavg +',' + hapsteravg + '\n')
output.close()
