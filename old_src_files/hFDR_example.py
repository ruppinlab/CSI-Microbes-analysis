import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import pandas as pd
#from rpy2.robjects import r, pandas2ri

# import the structSSI package
structSSI = importr("structSSI")

# robjects.r is the r instance

# creating rpy2 vectors
#res = robjects.StrVector(['abc', 'def'])
vec = robjects.FloatVector([0.01, 0.01, 0.5])
# to set the name - https://rpy.sourceforge.io/rpy2/doc-dev/html/vector.html#names
vec.names = robjects.StrVector(['family1', 'genus1', 'genus2'])

df = pd.DataFrame(data={"parent": ['family1', 'family1'], "child": ['genus1', 'genus2']})
# R matrices and arrays are just vectors with a dim attribute
l = df["parent"].tolist()
l = l.extend(df["child"].tolist())
v = robjects.StrVector(l)
#v = robjects.StrVector(['family1', 'family1', 'genus1', 'genus2'])
m = robjects.r['matrix'](v, nrow = 2)
# how to go from pandas dataframe to R matrix

# calling R functions
# rsum = robjects.r['sum']
# rsum(robjects.IntVector([1,2,3]))[0]

# out <- hFDR.adjust(vec, mat, alpha = 0.05)
hFDR = robjects.r["hFDR.adjust"]
out = hFDR(vec, m)
# how to subset s4class objects - https://rpy2.github.io/doc/latest/html/generated_rst/s4class.html#
print(out)
adj_pvals = out.slots['p.vals'][1]
