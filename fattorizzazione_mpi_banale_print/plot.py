import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import sys

file_name = sys.argv[1]
data1 = loadtxt("./" + file_name)
NUM=data1[:,0]  #  
TIME=data1[:,1]  # 

fig = plt.figure() # 
top = fig.add_subplot(111) # 1 riga, 1 colonna, figura 1

top.set_title('BRUTE FORCE')
top.grid()
top.set_xlabel('n')
top.set_ylabel('time')
top.plot(NUM, TIME)
#top.text(200,35,"Steps: "+str(int(STEPS))) 

plt.savefig('fatt-brute-force.pdf')
#plt.show()
