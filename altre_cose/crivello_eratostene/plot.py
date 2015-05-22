import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import sys

prime_file = sys.argv[1]
eratostene_file = sys.argv[2]
data1 = loadtxt("./" + prime_file)
data2 = loadtxt("./" + eratostene_file)

# prime
NUM1=data1[:,0]
TIME1=data1[:,1]   

# eratostene
NUM2=data2[:,0]
TIME2=data2[:,1]   

fig = plt.figure() # 
top    = fig.add_subplot(211) # 2 rows, 1 col, plot 1
bottom = fig.add_subplot(212) # 2 rows, 1 col, plot 2

top.set_title('ERATOSTENE-THREAD')
top.grid()
top.set_ylabel('time (s)')
top.plot(NUM1, TIME1, 'b-')

bottom.set_title('ERATOSTENE')
bottom.grid()
bottom.set_xlabel('n')
bottom.set_ylabel('time (s)')
bottom.plot(NUM2, TIME2, 'r-')

plt.subplots_adjust(hspace=0.5)

plt.savefig('prime-benchmark.pdf')
#plt.show()
