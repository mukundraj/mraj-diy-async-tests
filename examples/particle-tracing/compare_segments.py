'''
Script for comparing segments generated using iexchange and exchange.

'''


fname_ex = "./exchange.txt"
fname_iex = "./iexchange.txt"

with open(fname_ex) as f:
    content0 = f.readlines()
    content0 = [x.split() for x in content0]

sc0 = sorted(content0, key=lambda x: (float(x[0]), float(x[1]), float(x[2])))

with open(fname_iex) as f:
    content1 = f.readlines()
    content1 = [x.split() for x in content1]

sc1 = sorted(content1, key=lambda x: (float(x[0]), float(x[1]), float(x[2])))


if sc0 == sc1:
    print "Segments are identical"
else:
    print "Segments are not identical"


