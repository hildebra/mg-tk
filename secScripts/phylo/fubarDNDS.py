import sys, json, gzip
#from https://github.com/veg/hyphy/issues/1063

with gzip.open(sys.argv[1],"r") as fh:
	fubar1 = json.load (fh)

#estimate dnds
fubar = fubar1["grid"]
ds = sum (r[0] * r[2] for r in fubar)
dn = sum (r[1] * r[2] for r in fubar)



#print ("ds = %g dn = %g" % (ds, dn))
#print ("Weights")
#print ("\tdS > dN : %g" % sum (r[2] for r in fubar if r[0] > r[1]))
#print ("\tdS = dN : %g" % sum (r[2] for r in fubar if r[0] == r[1]))
#print ("\tdS < dN : %g" % sum (r[2] for r in fubar if r[0] < r[1]))

#estimate pos/neg selec
fb = fubar1["MLE"]["content"]["0"]
pos=0
neg=0
for r in fb:
	if r[3]>0.9:
		neg+=1
	if r[4]>0.9:
		pos+=1
ds2 = sum (r[0] for r in fb) / fubar1["input"]["number of sites"]
dn2 = sum (r[1] for r in fb) / fubar1["input"]["number of sites"]

#print ("pos %g neg %g selected sites" % (pos, neg))

print ("\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g" % 
	(fubar1["input"]["number of sequences"],fubar1["input"]["number of sites"] , neg, pos, dn, ds, dn2, ds2) )
