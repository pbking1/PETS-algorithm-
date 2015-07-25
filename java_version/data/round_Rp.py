import math
f = open("Rp.csv")
f1 = open("Rp1.csv", "w")
for i in f.readlines():
	a = float(i.strip("\n"))
	aa = round(a, 4)
	#print round(a, 4)
	if(len(str(aa)) == 4):
		print str(aa) + "00\n"
		f1.write(str(aa) + "00\n")
	elif(len(str(aa)) == 3):
		f1.write(str(aa) + "000\n")
	elif(len(str(aa)) == 5):
		print str(aa) + "0\n"
		f1.write(str(aa) + "0\n")
	else:
		print str(aa) + "\n"
		f1.write(str(aa) + "\n")
