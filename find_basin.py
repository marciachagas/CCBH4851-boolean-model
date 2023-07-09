import sys, re
import numpy as np
import matplotlib.pyplot as plt

# check stationary genes in boolean evolution for last n time steps:
nsteps=1000

if len(sys.argv) < 2:
  print("Call this script parsing the file containing the trajectories calculated with ASSA-PBN")

file = sys.argv[1]

evolution = []
data = open(file, 'r')
for d in data:
    d = d.strip()
    state = re.split('\t', d)
    evolution.append(state)
data.close()
evolution.pop()
evolution.reverse()

if nsteps >= len(evolution):
    sys.exit("Number of steps is not valid, it is larger than the total time!")

genes = evolution.pop()
genes.pop()
ngenes = len(genes)

#Counter initialization
counter =  []
for j in range(ngenes): 
    counter.append(0)

# Counter incrementation
for t in range(nsteps): 
    for g in range(ngenes):
        if evolution[t][g] == evolution[t + 1][g]:
            counter[g] += 1
        else:
            counter[g] = 0 
    #print("#############################")
    #print(evolution[t], "   ", counter)
    #print(evolution[t+1])

# Normalization of the counter variable:
counter = [ c / nsteps for c in counter]

out = open("reached_basin.csv", "w")
for i in range(ngenes):        
    print(genes[i],"    ", counter[i], file = out)
out.close()

#generate a bar plot
y_pos = np.arange(len(genes))

# Create bars
plt.bar(y_pos, counter)

# Create names on the x-axis
plt.xticks(y_pos, genes)

# Show graphic
plt.show()


