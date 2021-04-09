import os, sys
import numpy as np
from fibresLib import get_afib_lf
#~ from numba import jit
#~ sys.path.append('/home/felipe/Dropbox/simuls18/fibresAuxLibrary')
from libSolverGP import *


def array2str(x):
	s = ''
	
	for y in x:
		s = s + str(y) + ' '
	
	return s


def getPKfromParam(paramFileName, pos_Vf = 4, pos_afib = 6, pos_sf = 13,NdimE = 2, engine = 'python'): # Fortran convention
	
	if(engine == 'python'):
		param = getFromParam(paramFileName)	
		if(NdimE == 2):
			return getPKfromParam2D_faster(param, pos_Vf, pos_afib, pos_sf)
		elif(NdimE == 3):
			return getPKfromParam3D(param, pos_Vf, pos_afib, pos_sf)
	elif(engine == 'fortran'):
		return getPKfromParam_fortranWrapper(paramFileName, pos_Vf, pos_afib, pos_sf, NdimE)
	else:
		print('engine not founded')

def getPKfromParam3D(param, pos_Vf, pos_afib, pos_sf): 
	
	NdimE = 3
	Vf = param[1,:,pos_Vf]
	sumVf = Vf.sum()
	af = param[1,:,pos_afib:pos_afib + NdimE]
	
	sf = param[:,:,pos_sf:pos_sf+NdimE]
	
	ntime = len(sf)
	nfib = len(af)
	
	PK = np.zeros((ntime,NdimE,NdimE)) # Fortran convention
	
	for i in range(ntime):
		for j in range(nfib):
			PK[i,:,:] = PK[i,:,:] + Vf[j]*np.outer(sf[i,j,:],af[j,:])
	
	PK = PK/sumVf
	
	return PK
	

def getPKfromParam2D(param, pos_Vf, pos_afib, pos_sf): 
	
	Vf = param[1,:,pos_Vf]
	sumVf = Vf.sum()
	af = param[1,:,pos_afib:pos_afib + 2]
	
	sf = param[:,:,pos_sf:pos_sf+2]
	
	ntime = len(sf)
	nfib = len(af)
	
	PK = np.zeros((ntime,4)) # Fortran convention
	
	for n in range(ntime):
		print('homogenising time :', n)
		for j in range(nfib):
			PK[n,0] = PK[n,0] + Vf[j]*af[j,0]*sf[n,j,0]
			PK[n,1] = PK[n,1] + Vf[j]*af[j,0]*sf[n,j,1]
			PK[n,2] = PK[n,2] + Vf[j]*af[j,1]*sf[n,j,0]
			PK[n,3] = PK[n,3] + Vf[j]*af[j,1]*sf[n,j,1]
	
	PK = PK/sumVf
	
	return PK
	
	
def getPKfromParam_fortranWrapper(paramFile, pos_Vf, pos_afib, pos_sf,NdimE):
	NTimeSteps,NParamCell, sumSize, nParamGroups, sizeGroups = getInfoParam(paramFile)
	
	CONVERTER = "/home/felipe/Dropbox/SOLVERS/softsAux/converters2018/converter.x"

	Nelem = nParamGroups - 1
	homogParam = [Nelem, NTimeSteps, NParamCell, NdimE, pos_Vf , pos_afib, pos_sf]
	
	
	inputConv = " 6 '" + paramFile + "' '" + array2str(homogParam) + "'"
	
	print(inputConv)
	os.system(CONVERTER + inputConv) 
	
	
	PK = np.loadtxt('Homog.txt').reshape((NTimeSteps,NdimE*NdimE))
	
	return PK


def getPKfromParam2D_faster(param, pos_Vf, pos_afib, pos_sf): # not so fast
	
	Vf = param[1,:,pos_Vf]
	sumVf = Vf.sum()
	af = param[1,:,pos_afib:pos_afib + 2]
	
	sf = param[:,:,pos_sf:pos_sf+2]
	
	nfib = len(af)
	ntime = len(sf)
		
	PK = np.einsum('lk,ilj',Vf.reshape(nfib,1)*af/sumVf,sf).reshape(ntime,4)
	
	return PK

	
def getVolumeFraction(paramFile, pos_Vf, VolOmegaMu):
	data = getFromParamSpecific(paramFile, [0] , [pos_Vf])
	sumVf = np.sum(data[0,:,0])
	return sumVf/VolOmegaMu
	

#~ def getPKfromParam(param, pos_Vf = 4, pos_afib = 6, pos_sf = 13):
	
	#~ Vf = param[1,:,pos_Vf]
	#~ sumVf = Vf.sum()
	#~ af = param[1,:,pos_afib:pos_afib + 2]
	
	#~ sf = param[:,:,pos_sf:pos_sf+2]
	
	#~ ntime = len(sf)
	#~ nfib = len(af)
	
	#~ PK = np.zeros((ntime,4))
	
	#~ PKtemp = np.zeros((2,2))
	#~ for t in range(ntime):
		#~ Pktemp = 0.0
		#~ for f in range(nfib):
			
			#~ for i in range(2):
				#~ for j in range(2):
					#~ PKtemp[i,j] = PKtemp[i,j] + Vf[f]*sf[t,f,i]*af[f,j]

	
		#~ PK[t,:] = PKtemp.flatten() # C convention
#		PK[t,:] = np.transpose(PKtemp).flatten() # Fortran convention
	
	#~ PK = PK/sumVf
	
	
	#~ return PK


def testRecruitmentActivation(param, pos_lfa = 5, pos_stretch = 23, tol = 0.0, op = 1 ):
	
	# if op == 0 , comparates strecht with lfa, if op == 1, the stretch is assumed to be a flag of activation
	
	# times x fibres x position
	lfa = param[0,:,pos_lfa]
	stretch = param[:,:,pos_stretch]

	print(stretch)

	ntime = len(stretch)
	nfib = len(stretch[0,:])
	
	percentRecruited = np.zeros(ntime)
	
	funcTest0 = lambda x,y : x > y
	funcTest1 = lambda x,y : x > 0.0
	
	funcTest = None
	if(op ==  0):
		funcTest = funcTest0
	elif(op == 1):
		funcTest = funcTest1
		
	
	for i in range(ntime):
		count = 0
		for j in range(nfib):
			if(funcTest(stretch[i,j], lfa[j] - tol)):
				count = count + 1

		percentRecruited[i] = float(count)/float(nfib)

	
	
	return percentRecruited


def getNormFluctuation(param, dataout, elemFib, addPts = 2 , Ipos_uf = 0, Ipos_Vf = 4,NdimE = 2):
	
	ntime = len(dataout)
	nfib = len(elemFib)
	
	normUf = np.zeros(ntime)
	
	Uf = dataout[:,:-addPts,Ipos_uf:Ipos_uf+NdimE]
	Vf = param[1,:,Ipos_Vf]
	sumVf = np.sum(Vf)
	
	for i in range(ntime):
		for j in range(nfib):
			p1 = elemFib[j,0] 
			p2 = elemFib[j,1]
			
			#~ normUf1_sq = np.dot(Uf[i,p1,:],Uf[i,p1,:])
			#~ normUf2_sq = np.dot(Uf[i,p2,:],Uf[i,p2,:])
			
			normUf1 = np.linalg.norm(Uf[i,p1,:])
			normUf2 = np.linalg.norm(Uf[i,p2,:])
			
			normUf[i] = normUf[i] + (Vf[j]/2.0)*( normUf1 + normUf2) 
		
		normUf[i] = normUf[i]/sumVf
		
	return normUf

def getLambda(X, param, dataout, elemFib , Ipos_uT, Ipos_uf, Ipos_Lf , addPts = 2):
	
	ntime = len(dataout)
	nfib = len(elemFib)
	
	lamb = np.zeros((ntime,nfib))
	lamb_bar = np.zeros((ntime,nfib))
	
	UT = dataout[:,:-addPts,Ipos_uT:Ipos_uT+2]
	Uf = dataout[:,:-addPts,Ipos_uf:Ipos_uf+2]
	Lf = param[1,:,Ipos_Lf] 
	
	X = X[:-addPts,0:2]

	xT = np.zeros(UT.shape)
	xTbar = np.zeros(UT.shape)

	for i in range(ntime):
		xT[i,:,:] = X + UT[i,:,:]
		xTbar[i,:,:] = X + UT[i,:,:] - Uf[i,:,:]
	
	for i in range(ntime):
		for j in range(nfib):
			p1 = elemFib[j,0] 
			p2 = elemFib[j,1]
			
			delta_xT = xT[i,p1,:] - xT[i,p2,:]
			delta_xTbar = xTbar[i,p1,:] - xTbar[i,p2,:]
			
			lamb[i,j] = np.linalg.norm(delta_xT)/Lf[j] 
			lamb_bar[i,j] = np.linalg.norm(delta_xTbar)/Lf[j] 
		
	return lamb, lamb_bar







#~ def getNormFluctuation(dataout, addPts = 2 , pos_uf = 0):
	
	#~ ntime = len(dataout)
	
	#~ normUf = np.zeros(ntime)
	
	#~ Uf = dataout[:,:-addPts,pos_uf:pos_uf+2]

	#~ nNodes = len(Uf[0,:,:])
	
	#~ for i in range(ntime):
		#~ for j in range(nNodes):
			#~ normUf[i] = normUf[i] + np.linalg.norm(Uf[i,j,:])
		
		
	#~ return normUf


	
	
def getMinRestrictionResidual(param, dataout, elemFib, addPts = 2 , Ipos_uf = 0, Ipos_flag1 = 0 , Ipos_flag2 = 1, Ipos_Vf = 4 , Ipos_afib = 6, Ipos_Areaf = 3,  
			Ipos_Abar1 = 12, Ipos_Abar2 = 13, Ipos_normal1 = 14, Ipos_normal2 = 16, Ipos_nConnectFibre1 = 18, Ipos_nConnectFibre2 = 19, op = 1):
	
	ntime = len(dataout)
	nfib = len(elemFib)
	
	resMR = np.zeros((ntime,2,2))
	resMRnorm = np.zeros(ntime)
	
	resVolAvg = np.zeros((ntime,2))
	resVolAvgNorm = np.zeros(ntime)
	
	Uf = dataout[:,:-addPts,Ipos_uf:Ipos_uf+2]
	
	Vf = param[1,:,Ipos_Vf]
	flag1 = param[1,:,Ipos_flag1]
	flag2 = param[1,:,Ipos_flag2]
	
	
	if(op == 0):
		Abar1 = param[1,:,Ipos_Areaf]
		Abar2 = param[1,:,Ipos_Areaf]
		normal1 = param[1,:, Ipos_afib : Ipos_afib + 2]
		normal2 = - param[1,:, Ipos_afib : Ipos_afib + 2]
		nConnect_frac1 = np.ones(nfib)
		nConnect_frac2 = np.ones(nfib)

	elif(op == 1):
		Abar1 = param[1,:,Ipos_Abar1]
		Abar2 = param[1,:,Ipos_Abar2]
		normal1 = param[1,:, Ipos_normal1 : Ipos_normal1 + 2]
		normal2 = param[1,:, Ipos_normal1 : Ipos_normal1 + 2]
		nConnect_frac1 = 1.0/param[1,:,Ipos_nConnectFibre1]
		nConnect_frac2 = 1.0/param[1,:,Ipos_nConnectFibre2]
		
	
	for i in range(ntime):
		for j in range(nfib):
			p1 = elemFib[j,0] 
			p2 = elemFib[j,1]
			
			if(flag1[j]>0):
				resMR[i,:,:] =  resMR[i,:,:] + Abar1[j]*nConnect_frac1[j]*np.outer(Uf[i,p1,:],normal1[j,:])
				
			if(flag2[j]>0):
				resMR[i,:,:] =  resMR[i,:,:] + Abar2[j]*nConnect_frac2[j]*np.outer(Uf[i,p2,:],normal2[j,:])
			
			
			resVolAvg[i,:] =  resVolAvg[i,:] + 0.5*Vf[j]*( Uf[i,p1,:] + Uf[i,p2,:] )
			
		resVolAvgNorm[i] =  np.linalg.norm(resVolAvg[i,:])
		resMRnorm[i] = np.linalg.norm(resMR[i,:,:].flatten())
			
			
		
	return resMR , resMRnorm, resVolAvg , resVolAvgNorm

unitVector2AngleSingle = lambda af : 180.0*np.arctan(af[1]/af[0])/np.pi
			
def unitVector2Angles(af):
	nf = len(af)
	angles = np.zeros(nf)
	
	for f in range(nf):
		angles[f] = unitVector2AngleSingle(af[f,:])

	return angles


def getAngleStats(dataout,regularAnglesFileName,X,elemFib,Ipos_u = 2, addPts = 2 ):
	
	ntime = len(dataout)
	nfib = len(elemFib)
	
	U = dataout[:,:-addPts,Ipos_u:Ipos_u+2]
	x = X[:-addPts] + U
	
	regularAngles = np.loadtxt(regularAnglesFileName)
	
	regularAnglesFlagP = regularAngles>0.0
	regularAnglesFlagM = np.logical_not(regularAnglesFlagP)
	
	
	anglesMeanP = np.zeros(ntime)
	anglesMeanM = np.zeros(ntime)
	anglesStdP = np.zeros(ntime)
	anglesStdM = np.zeros(ntime)
	
	
	anglesAux = np.zeros(nfib)
	
	for i in range(ntime):

		for j in range(nfib):
			anglesAux[j] = unitVector2AngleSingle(get_afib_lf(x[i,:,:],elemFib[j,:])[0])
		
		anglesMeanP[i] = np.mean(anglesAux[regularAnglesFlagP]) 
		anglesMeanM[i] = np.mean(anglesAux[regularAnglesFlagM])	
		anglesStdP[i] = np.std(anglesAux[regularAnglesFlagP]) 
		anglesStdM[i] = np.std(anglesAux[regularAnglesFlagM])
			
	
	return anglesMeanP , anglesMeanM , anglesStdP , anglesStdM 
	
	

def getAnglesDataframe(dataout,regularAnglesFileName,X,elemFib,idFile, Ipos_u = 2, addPts = 2 ):
	import pandas as pd	
	
	ntime = len(dataout)
	nfib = len(elemFib)
	
	U = dataout[:,:-addPts,Ipos_u:Ipos_u+2]
	x = X[:-addPts] + U
	
	regularAngles = np.loadtxt(regularAnglesFileName)
	
	regularAnglesFlag = regularAngles>0.0
	
	del regularAngles 

	time = np.zeros(nfib*ntime).astype('int')
	idFamily = np.zeros(nfib*ntime).astype('int')
	angles = np.zeros(nfib*ntime)
	
	k = 0
	for i in range(ntime):
		for j in range(nfib):
			angles[k] = unitVector2AngleSingle(get_afib_lf(x[i,:,:],elemFib[j,:])[0])
			time[k] = i
			
			if( regularAnglesFlag[j] ):
				idFamily[k] = 1
			else:
				idFamily[k] = 2
				angles[k] = - angles[k]
			
			k = k + 1

			
	d = {'angles' : angles , 'time' : time , 'idFamily' : idFamily , 'idFile' : (nfib*ntime)*[idFile]}
	
	data = pd.DataFrame(d)
	return data 
	
	

