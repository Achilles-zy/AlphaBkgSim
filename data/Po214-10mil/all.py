#!/usr/bin/env python
# -*- coding:utf-8 -*-
#/usr/bin/env python3
import os

def writeJson(path,maxid,dtype,srcid,energy,thstep,beam):
	for i in range(0,maxid):
		jsonname=str(dtype)+"_"+str(srcid)+"_"+str(energy).replace('.','')+"_"+str(i)
		dirname=str(dtype)+"_"+str(srcid)+"_"+str(energy).replace('.','')
		if os.path.exists(path+dirname)==False:
			os.mkdir(path+dirname)
		with open(path+dirname+"/"+jsonname+".json","w+",encoding="utf-8") as f:
			f.write("{\n")
			f.write('  "FileName": '+'"'+jsonname+'"'+","+"\n")
			f.write("\n")
			f.write('  "DetectorType": '+'"'+str(dtype)+'"'+","+"\n")
			f.write("\n")
			f.write('  "SourceID": '+str(srcid)+","+"\n")
			f.write("\n")
			f.write('  "SourceEnergy": '+str(energy)+","+"\n")
			f.write("\n")
			f.write('  "ID": '+str(i)+","+"\n")
			f.write("\n")
			f.write('  "DeadLayerThickness": '+str(i*thstep+0.000001)+","+"\n")
			f.write("\n")
			f.write('  "Beam": '+str(beam)+"\n")
			f.write("\n")
			f.write("}")
	return dirname

def pwriteJson(path,maxid,dtype,dgeo,srcid,energy,thstep,beam):
	for i in range(0,maxid):
		jsonname="p_"+str(dtype)+"_"+str(dgeo)+"_"+str(srcid)+"_"+str(energy).replace('.','')+"_"+str(i)
		dirname="p_"+str(dtype)+"_"+str(dgeo)+"_"+str(srcid)+"_"+str(energy).replace('.','')
		if os.path.exists(path+dirname)==False:
			os.mkdir(path+dirname)
		with open(path+dirname+"/"+jsonname+".json","w+",encoding="utf-8") as f:
			f.write("{\n")
			f.write('  "FileName": '+'"'+jsonname+'"'+","+"\n")
			f.write("\n")
			f.write('  "DetectorType": '+'"'+str(dtype)+'"'+","+"\n")
			f.write("\n")
			f.write('  "DetectorGeo": '+'"'+str(dgeo)+'"'+","+"\n")
			f.write("\n")
			f.write('  "SourceID": '+str(srcid)+","+"\n")
			f.write("\n")
			f.write('  "SourceEnergy": '+str(energy)+","+"\n")
			f.write("\n")
			f.write('  "ID": '+str(i)+","+"\n")
			f.write("\n")
			f.write('  "DeadLayerThickness": '+str(1)+","+"\n")
			f.write("\n")
			f.write('  "pLayerThickness": '+str(i*thstep)+","+"\n")
			f.write("\n")
			f.write('  "Beam": '+str(beam)+"\n")
			f.write("\n")
			f.write("}")
	return dirname


#path,maxid,dtype,srcid,energy,thstep,beam
#directory=writeJson("./jsonfiles/",8,"Natural",1,5.304,0.005,1000000)

#path,maxid,dtype,dgeo,srcid,energy,thstep,beam
pdirectory=pwriteJson("./jsonfiles/",11,"Enriched","flatBEGe",14,0,0.0001,10000000)

#Enriched/Natural
#MeV
#mm
# SrcID:
# 1: Vertical upper center
# 2: Vertical side center
# 3: Angular-even upper center
# 4: Angular-even side center
# 5: p+ surface
# 6: p+

# exe="./exampleB1 ./jsonfiles/"+directory+"/"

# inputfiles=os.listdir("./jsonfiles/"+directory+"/")
# print(inputfiles)
# if os.path.exists("./"+directory)==False:
# 	os.mkdir("./"+directory)

# for i in range(0,len(inputfiles)):
# 	cur_file=inputfiles[i]
# 	print(inputfiles[i])
# 	os.system(exe+cur_file)

exe="./exampleB1 ./jsonfiles/"+pdirectory+"/"

inputfiles=os.listdir("./jsonfiles/"+pdirectory+"/")
print(inputfiles)
if os.path.exists("./"+pdirectory)==False:
	os.mkdir("./"+pdirectory)

for i in range(0,len(inputfiles)):
	cur_file=inputfiles[i]
	print(inputfiles[i])
	os.system(exe+cur_file)
