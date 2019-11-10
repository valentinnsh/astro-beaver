import numpy as np
import math

f = open("data.dat", "w")

N = 13
f.write("#  " + str(N) + '\n')

N = 2**math.ceil(np.log2(N))

for i in range(0, N):
    tmp = np.random.rand(1, 2)*2-1
    f.write(str(tmp[0, 0]) + " " + str(tmp[0, 1]) + "\n")
