#!/usr/bin/env python

from csvParse import *
#from numpy import zeros,tril,corrcoef,transpose,polyfit,mean
#from numpy.random import shuffle
#from scipy.stats import linregress
#from scipy.signal import residue
#import time

# parse our csv files, put into matrices
path = "C:/CEL/Mantel/"
xfile = "GDMatrix1000.csv"
yfile = "CDMatrix_Landscape1000.csv"
zfile = "CDMatrix_Distance1000.csv"
csv2zt(path,xfile)
csv2zt(path,yfile)
csv2zt(path,zfile)

'''
yfile = "CDMatrix_Landscape1000.csv"
zfile = "CDMatrix_Distance1000.csv"
xmat = csvParse(path, xfile)
ymat = csvParse(path, yfile)
zmat = csvParse(path, zfile)

# get lower tri of the matrices created above
lowerxmat = tril(xmat, -1)
lowerymat = tril(ymat, -1)
lowerzmat = tril(zmat, -1)

# get num rows and cols of matrix
matDimension = xmat.shape
matDimension = matDimension[0]  # [0] == [1] when using a square matrix, so we only need one value

# compute the number of values within the lower tri of the matrix
#x = 0
#for i in range(1,matDimension):
#  x = x + i

# algebraic form of the above summation
x = matDimension*(matDimension-1)/2

# create vectors of length computed above
a = zeros((x))
b = zeros((x))
c = zeros((x))

# reshape input matrices into a vector of all values
lowerxmat = lowerxmat.reshape(matDimension**2)
lowerymat = lowerymat.reshape(matDimension**2)
lowerzmat = lowerzmat.reshape(matDimension**2)

# grab all values of lower tri and put into vectors
index = 0
for i in range(1,matDimension):
  for j in range(0,i):
    k = (i*matDimension)+j
    a[index] = lowerxmat[k]
    b[index] = lowerymat[k]
    c[index] = lowerzmat[k]
    index = index+1

xresid = polyfit(c,a,1)
yresid = polyfit(c,b,1)
print "======== polyfit, degree 1 polynomial ========"
print "xresid: ",
print xresid
print "yresid: ",
print yresid

print ""

xresid = linregress(c,a)
yresid = linregress(c,b)
print "======== linregress ========"
print "xresid: ",
print xresid
print "yresid: ",
print yresid

print ""

#xsamplemean = mean(a)
#print "sample mean: ",
#print xsamplemean

slope = xresid[0]
intercept = xresid[1]

print "slope: ",
print slope
print "intercept: ",
print intercept

#actual value
measx = a[0]
#approx value
approx = (slope*measx)+intercept
#residual
residual = measx-approx
print "Measured: ", measx
print "approx: ", approx
print "residual: ", residual

print residue(c,a)

# compute normalized covariance of vectors a and b
mancor = corrcoef(a,b)[0][1]

# create a vector of zeros of length nperms... refactor this so nperms can be fed in elsewhere
nperms = 199
xperms = zeros((matDimension,matDimension))
mancorperms = zeros((nperms))

# create index vector
ind = zeros((matDimension))
for i in range(1,matDimension+1):
  ind[i-1] = i

t1 = time.time()
print t1

# perform multiple permutations of the mantel test
for i in range(0,nperms):
  shuffle(ind)
  for j in range(0,matDimension):
    for k in range(0,matDimension):
      xperms[j][k] = xmat[ind[j]-1][ind[k]-1]
  # get lower tri of xperms
  lowerxperms = tril(xperms, -1)
  lowerxperms = lowerxperms.reshape(matDimension**2)
  # grab all values of lower tri and put into a vector
  # the 'a' vector above is no longer used and can be repurposed
  index = 0
  for q in range(1,matDimension):
    for r in range(0,q):
      s = (q*matDimension)+r
      a[index] = lowerxperms[s]
      index = index+1
  # compute normalized covariance of vectors lowerxperms and b (lowery
  mancorperms[i] = corrcoef(a,b)[0][1]
  print "perm: ", i, " m: ", mancorperms[i]

#compute pvalue
pvalue = 0.
for i in range(0,nperms):
  if abs(mancor)<abs(mancorperms[i]):
    pvalue = pvalue + 1.

#what if we do pvalue/nperms ???
#is this +1 needed?
pvalue = (pvalue + 1.)/(nperms + 1.)

print "pvalue: ", pvalue
print "mancor: ", mancor

t2 = time.time()
print t2
print "Total time elapsed:",
print t2-t1
'''