import math

sumLog=0

for i in range(1,200000):
    sumLog+=math.log10(i)
    print(i,"\t",sumLog)
