import matplotlib.pyplot as plt
import csv

x, y = [], []

#with open('result.dat') as f:
 #   lines = f.readlines()
  #  x = [line.split()[0] for line in lines]
   # y = [line.split()[1] for line in lines]

for line in open('result.dat', 'r'):
  values = [float(s) for s in line.split()]
  x.append(values[0])
  y.append(values[1])

# print(x, y)
plt.plot(x,y)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Interesting Graph\nCheck it out')
plt.show()


