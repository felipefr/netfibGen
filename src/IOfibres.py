import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import copy
import os, sys


strArray = lambda a : reduce(lambda x,y: str(x) + ' ' + str(y), a.flatten())

def writeInifileDefault(Ncoord,Ndof,op = 2):
	
	inifile=open("IniFile000.txt",'w')
	
	if(op=='Piola'): # Piola Style
		inifile.write("*Time\n0 0 0\n")
	
	elif(op=='SGP'):  #SGP Style
		inifile.write("*Initial Conditions\n")
	
	for j in range(Ncoord):
		for k in range(Ndof):
			inifile.write("0.0 ")
		inifile.write("\n")	

	inifile.close()
	
def writeParamDefault(Nfib,Nparam):

	paramfile=open("Param000.txt",'w')
	paramfile.write("*Parameter Groups\n")
	paramfile.write(str(Nfib+1) + "\n")
	
	paramfile.write("*Real Parameters\n")
	paramfile.write(str(0) + " ")
	for j in range(Nfib):
		paramfile.write(str(Nparam) + " ")
	paramfile.write("\n")

	for j in range(Nfib):
		for k in range(Nparam):
			paramfile.write("0.0 ")
		paramfile.write("\n")
	
	paramfile.write("*Integer Parameters\n")
	paramfile.write(str(0) + " ")
	for j in range(Nfib):
		paramfile.write(str(0) + " ")
	
	paramfile.close()


def writeParamCustom(Param):

	NparamGroups = len(Param)
	Nparam = len(Param[0])

	paramfile=open("Param000.txt",'w')
	paramfile.write("*Parameter Groups\n")
	paramfile.write(str(NparamGroups+1) + "\n")
	
	paramfile.write("*Real Parameters\n")
	paramfile.write(str(0) + " ")
	for j in range(NparamGroups):
		paramfile.write(str(Nparam) + " ")
	paramfile.write("\n")

	for j in range(NparamGroups):
		for k in range(Nparam):
			paramfile.write(str(Param[j,k]) + " ")
		paramfile.write("\n")
	
	paramfile.write("*Integer Parameters\n")
	paramfile.write(str(0) + " ")
	for j in range(NparamGroups):
		paramfile.write(str(0) + " ")
	
	paramfile.close()


def getColourPallette(net):

	nf = len(net.ElemFib)
	m = 0.8
	maxAf = np.max(net.Af)
	minAf = np.min(net.Af)
	
	print maxAf, minAf
	
	v = np.zeros(nf)
	
	colourPallette = []
	
	if(maxAf>minAf):
		v = (m/(maxAf - minAf))*(net.Af - minAf)
		
		for f in range(nf):
			colourPallette = colourPallette + [(v[f],v[f],v[f])]

	else:
		colourPallette = ['black']


	return colourPallette 
	
def writeFigNetwork(net,colourPallette=0,lw = 7.0,figNum = 1,filename = 'network.pdf'):
	
	if(colourPallette==0): 
		colourPallette = getColourPallette(net)
	elif(type(colourPallette) != type([]) ):
		colourPallette  = [colourPallette]
		
	if(net.ndim == 2):
		writeFigNetwork2D(net,colourPallette,lw,figNum,filename)
	elif(net.ndim == 3):	
		#~ writeFigNetwork3D(net,colourPallette,lw,figNum,filename)
		writeFigNetwork3D_mayavi(net)
		
	plt.tight_layout()
	plt.savefig(filename) 

def writeFigNetwork3D_mayavi(net):
	from mayavi import mlab
	
	mlab.figure(1, bgcolor=(0, 0, 0))
	mlab.clf()
	 
	Paux = net.P[net.flagNode>0,:]
	mlab.points3d(net.P[:, 0], net.P[:, 1], net.P[:, 2], net.flagNode , scale_factor=0.1, scale_mode='none', colormap='Blues', resolution=20)
	
	#~ mlab.points3d(Paux[:, 0], Paux[:, 1], Paux[:, 2], net.flagNode[net.flagNode>0] , scale_factor=0.1, scale_mode='none', colormap='Blues', resolution=20)
			
	for e in net.ElemFib:
		mlab.plot3d([net.P[e[0],0],net.P[e[1],0]],[net.P[e[0],1],net.P[e[1],1]],[net.P[e[0],2],net.P[e[1],2]])
		
	mlab.show()

def writeFigNetwork3D(net,colourPallette,lw,figNum,filename):

	k = 0
	Ncp = len(colourPallette)
	fig = plt.figure(figNum,(6,6))
	#ax = fig.add_subplot(111, projection='3d')
	ax = fig.gca(projection='3d')

	
	for f in net.ElemFib:
		X = net.P[f,:]
		ax.plot(X[:,0],X[:,1],X[:,2],linewidth = lw, color = colourPallette[k%Ncp])
		k = k + 1
	


def writeFigNetwork2D(net,colourPallette,lw,figNum,filename):

	k = 0
	Ncp = len(colourPallette)
	plt.figure(figNum,(6,6))
	for f in net.ElemFib:
		x = [net.P[f[0],0],net.P[f[1],0]]
		y = [net.P[f[0],1],net.P[f[1],1]]

		plt.plot(x,y, linewidth = lw, color = colourPallette[k%Ncp])
		k = k + 1
		

def writeParamFibres(net,Nparam,fignum = -1):
		
	print Nparam
	Param = np.zeros((net.nf,Nparam))
	
	# fortran to python convention
	Ipos_flag1 = 0
	Ipos_flag2 = 1 	
	Ipos_Lf = 2
	Ipos_Areaf = 3  
	Ipos_Vf = 4  
	Ipos_laf = 5 
	Ipos_af = 6
	
	if(net.ndim==2):	
		Ipos_yrel1 = 8
		Ipos_yrel2 = 10
		Ipos_Abar1 = 12
		Ipos_Abar2 = 13
		Ipos_normal1 = 14
		Ipos_normal2 = 16
		Ipos_nConnectFibre1 = 18
		Ipos_nConnectFibre2 = 19
		
		
		#~ Ipos_r0f = 20 # compatible with convention 0 of fibresMod
		Ipos_r0f = 24 # compatible with convention 1 of fibresMod

	elif(net.ndim==3):
		Ipos_yrel1 = 9
		Ipos_yrel2 = 12
		Ipos_Abar1 = 15
		Ipos_Abar2 = 16
		Ipos_normal1 = 17
		Ipos_normal2 = 20
		Ipos_nConnectFibre1 = 23
		Ipos_nConnectFibre2 = 24
		Ipos_r0f = 25
			
	# just to update
	net.updateDepedentVariables()
	
	ndim = net.ndim
	normal_bar = net.get_normal_bar()
	yG = net.get_yG()
	Bten = net.get_Bten()
	BtenInv = np.linalg.inv(Bten)
	
	BtenFile=open("Bten_yG.txt",'w')
	
	BtenFile.write( strArray(Bten) + '\n')
	BtenFile.write( strArray(BtenInv) + '\n')
	BtenFile.write( str(np.sum(net.Vf)) + '\n')
	BtenFile.write( strArray(yG) + '\n')


	for f in range(net.nf):
		p1 = net.ElemFib[f][0]
		p2 = net.ElemFib[f][1]
		
		Param[f,Ipos_flag1] = net.flagNode[p1]
		Param[f,Ipos_flag2] = net.flagNode[p2]
		Param[f,Ipos_Lf] = 	net.Lf[f]
		Param[f,Ipos_Areaf] = net.Af[f]
		Param[f,Ipos_Vf] = net.Vf[f]
		Param[f,Ipos_laf] = net.laf[f]
		Param[f,Ipos_af:Ipos_af+ndim] = net.af[f,:]
		Param[f,Ipos_yrel1:Ipos_yrel1 + ndim] = net.P[p1,:] - yG
		Param[f,Ipos_yrel2:Ipos_yrel2 + ndim] = net.P[p2,:] - yG 
		Param[f,Ipos_Abar1] = net.Abar[p1]
		Param[f,Ipos_Abar2] = net.Abar[p2] 
		Param[f,Ipos_normal1:Ipos_normal1 + ndim] = net.normal[p1,:] - normal_bar
		Param[f,Ipos_normal2:Ipos_normal2 + ndim] = net.normal[p2,:] - normal_bar
		Param[f,Ipos_nConnectFibre1] = len(net.InvConnect[p1])
		Param[f,Ipos_nConnectFibre2] = len(net.InvConnect[p2])		
		Param[f,Ipos_r0f] = net.r0f[f]		
		
		
	writeParamCustom(Param)
	
	
	if(fignum != -1): 
		os.system("cp Param000.txt  Param000_" + str(fignum) + ".txt ")



def writeMesh(X,Elem,auxElem,Ndof,Nsubstep):
	
	Ncoord=len(X)
	ElemTotal = Elem + auxElem
	NelemGroups= len(ElemTotal)
	Ndim = len(X[0,:])
	

	mesh=open("Mesh.txt",'w')
	mesh.write("*NODAL DOFs\n")
	mesh.write(str(Ndof) + "\n")
	mesh.write("*DIMEN\n")
	mesh.write("3\n")
	mesh.write("*COORDINATES\n")
	mesh.write(str(Ncoord) + "\n\n")
	
	if(Ndim == 2):
		for i in range(Ncoord):
			mesh.write(format(X[i,0],'.10e') + " " + format(X[i,1],'.10e') + " 0.0 \n")
	elif(Ndim == 3):
		for i in range(Ncoord):
			mesh.write(format(X[i,0],'.10e') + " " + format(X[i,1],'.10e') + " " + format(X[i,2],'.10e') + "\n")

	mesh.write("\n")
	mesh.write("*ELEMENT GROUPS\n")
	mesh.write(str(NelemGroups) + "\n")
	
	for i in range(NelemGroups):
		mesh.write(str(i + 1) + " " + str(len(ElemTotal[i])) + " Generic\n")

	mesh.write("\n")
	
	for elemGroup in ElemTotal:
		for e in elemGroup:
			mesh.write(str(len(e)) + "\n")

	mesh.write("\n")

	mesh.write("*INCIDENCE\n")
	
	for elemGroup in ElemTotal:
		for e in elemGroup:
			for v in e:
				mesh.write(str(v+1) + " ") # fortran convention
					
			mesh.write("\n")
	
	mesh.write("\n")
	mesh.write("*ELEMENT TYPE\n")

	k = 0
	for elemGroup in ElemTotal:
		k = k + 1
		for e in elemGroup:
			mesh.write(str(k) + "\n") 

	mesh.write("\n")
	mesh.write("*ELEMENT MAT\n")
	
	# just the elements in Elem have material associated
	k = 2
	for elemGroup in Elem:
		for e in elemGroup:
			mesh.write(str(k) + "\n")
			k = k + 1

	for elemGroup in auxElem:
		for e in elemGroup:
			mesh.write(str(1) + "\n")
	
	mesh.write("\n")
	mesh.write("*DIRICHLET CONDITIONS\n")
	
	for i in range(Nsubstep):
		for j in range(Ncoord):
			for k in range(Ndof):
				mesh.write("0 ")
			mesh.write("\n")
		mesh.write("\n")
		
			
	for i in range(Nsubstep):
		for j in range(Ncoord):
			for k in range(Ndof):
				mesh.write("0.0 ")
			mesh.write("\n")
		mesh.write("\n")
		
	mesh.close()
