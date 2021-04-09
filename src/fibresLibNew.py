#~ from generator import *
#~ from delaunay import *
from scipy.spatial import Delaunay
import numpy as np
import matplotlib.pyplot as plt
import copy
import os, sys
from IOfibres import *


class Network:
	
	def __init__(s): # default constructor 
		
		s.altenancyPointArea = [1,1]
		
		s.asymFac = 0 
		s.meanTheta = 0.25*np.pi
		s.nFibPrevision = 20
		
		s.ndim = 2
	
		s.P=[]
		s.ElemFib=[]
		s.ElemAngle = []
		s.ElemAngle2 = [] # second family of angles, it may be interesting
		s.ElemAngle3 = [] # third family of angles, it may be interesting
		
		s.delauTri = [] # it is optional, used to model the matrix, it improves convergence
		
		s.nx = 0
		s.ny = 0
		s.nf = 0
		s.np = 0
		s.flagNode = []
		s.flagFibres = []
		s.InvConnect = []
		s.InvConnectSign = []
		s.InvConnectPoint = []
		
		s.Lf = []
		s.Af = []
		s.Vf = []
		s.laf = []
		s.r0f = []
		s.af = []
		s.yrel = []
		s.yG = []
		s.normal_bar = []
		s.Bten = []
		s.BtenInvT = []
	
		s.normal = []
		s.Abar = []

	def allocFibreProperties(s):
		s.Af = np.zeros(s.nf)
		s.r0f = np.zeros(s.nf)
		s.laf = np.zeros(s.nf)
		
	def updateDepedentVariables(s):
		
		s.set_af_Lf_Vf()
		s.setNormalAndAbar()
		
		
	def setFlagNodes(s):
		s.flagNode = np.zeros(s.np,dtype = 'int')
		
		flagMatAux = np.array([[5,6,1],[8,7,3],[4,2,0]])
		
		for i in range(s.np):
			pintx = classify(s.P[i,0])
			pinty = classify(s.P[i,1])
			s.flagNode[i] = flagMatAux[pinty,pintx]
			


	def createAngleElements(s,maxAnglePerJoint = 6):
		
		#~ thetaCrit = 120.0 # more than this value desconsider the fibre pair to the angle
		
		for i in range(s.np):
			
			listAngle_i = []
			listTheta_i = []
			for j in range(len(s.InvConnect[i])):
				fj = s.InvConnect[i][j]
				
				for k in range(j+1,len(s.InvConnect[i])): # jump the same fibre and eliminate repetition in combination
					fk = s.InvConnect[i][k]
					
					sjk = s.InvConnectSign[i][j]*s.InvConnectSign[i][k] # changes the sign accordingly to get the right angle
					
					cosTheta = np.dot(sjk*s.af[fj],s.af[fk])
					# in case of some numerical errors
					if(cosTheta<-1.0):
						cosTheta = -1.0
					elif(cosTheta>1.0):
						cosTheta = 1.0
					
					theta = np.rad2deg(np.arccos(cosTheta)) 
					
					listAngle_i = listAngle_i + [[ i, s.InvConnectPoint[i][j], s.InvConnectPoint[i][k] ]]
					listTheta_i = listTheta_i + [theta]
				
			
			smallerTheta_indexes = np.argsort(listTheta_i)
			
					
			for ind in range(min(maxAnglePerJoint,len(listTheta_i))):	
				s.ElemAngle = s.ElemAngle + [listAngle_i[smallerTheta_indexes[ind]]]	
			
			
			for ind in range(maxAnglePerJoint,len(listTheta_i)):	
				s.ElemAngle2 = s.ElemAngle2 + [listAngle_i[smallerTheta_indexes[ind]]]



	def createAngleElements_threeGroups(s):
		
		# preparing 
		eps = 0.001
		s.removeFibres(combinedCriterium(VertHorizFibresCriterium,InSquareFibresCriterium),[[],[eps,eps,1.0-2*eps,1.0-2*eps]])

		s.setFlagsAndConnectivity()

		s.cleanDangledPoints()
		s.setFlagNodes()
		s.setFlagsAndConnectivity()

		s.allocFibreProperties()

		# dummy properties
		s.setProperty('A',1,[1.0])
		s.setProperty('la',1,[1.0])
		s.setProperty('r0',1,[5.0])

		s.updateDepedentVariables()

		s.setFlagNodes()
		
		# ----------------------- Angle creation itself ------------------------------------
		
		thetaCrit = 150.0*np.pi/180.0 # more than this value desconsider the fibre pair to the angle
		thetaCrit2 = 90.0*np.pi/180.0 # more than this value desconsider the fibre pair to the angle
		
		# GROUP 1 : internal angles ( thetaCrit = 150 degrees)
		# GROUP 2 : for angles of fibres on boundary
		# GROUP 3 : for angles of fibres on boundary
		
		for i in range(s.np):
			
			for j in range(len(s.InvConnect[i])):
				fj = s.InvConnect[i][j]
				
				for k in range(j+1,len(s.InvConnect[i])): # jump the same fibre and eliminate repetition in combination
					
					triple = [ i, s.InvConnectPoint[i][j], s.InvConnectPoint[i][k] ]					
					fk = s.InvConnect[i][k]
					sjk = s.InvConnectSign[i][j]*s.InvConnectSign[i][k] # changes the sign accordingly to get the right angle
					cosTheta = np.dot(sjk*s.af[fj],s.af[fk])
					theta = np.arccos(cosTheta) 
					
					
					if( s.flagNode[i] > 0 ) : # at boundary
						if(bool(s.flagNode[triple[1]] < 1) != bool(s.flagNode[triple[2]] < 1) and theta < thetaCrit2):
							s.ElemAngle3 = s.ElemAngle3 + [triple]
						elif(s.flagNode[triple[1]] < 1 and s.flagNode[triple[2]] and theta < thetaCrit):
							s.ElemAngle = s.ElemAngle + [triple]
					else:
						if(theta < thetaCrit):
							s.ElemAngle = s.ElemAngle + [triple]
						else:
							s.ElemAngle2 = s.ElemAngle2 + [triple]
		
		# --------------------------------------------------------------------------------------  
		
		s.updateDepedentVariables()
		
	

				
	# build inverse connectivity (and its sign (-1 or 1, starting or ending) ) and classify all points and lines as the follwings categories
	# for points : -1 = neighbour to boundary, 0 =  very interior, 1 = bottom, 2 = right, 3 = top, 4 = left
	# 5 = bottom/left , 6 = bottom/right, 7 = right/top, 8 = top/left
	# for lines : 0 = interior fibres, 1 = neighbourhood fibres.
	def setFlagsAndConnectivity(s):

		s.InvConnect = s.np*[[]]
		s.InvConnectSign = s.np*[[]]
		s.InvConnectPoint = s.np*[[]]
		
		for i in range(s.nf):
			f = s.ElemFib[i]
			s.InvConnect[f[0]] = s.InvConnect[f[0]] + [i]
			s.InvConnect[f[1]] = s.InvConnect[f[1]] + [i]
			s.InvConnectSign[f[0]] = s.InvConnectSign[f[0]] + [-1]
			s.InvConnectSign[f[1]] = s.InvConnectSign[f[1]] + [1]						

		s.flagFibres = np.zeros(s.nf,dtype = 'int')
	
		for i in range(s.np):
			for j in range(len(s.InvConnect[i])):
				fij = s.InvConnect[i][j]
				f = s.ElemFib[fij]
				s.flagFibres[fij] = 1
				sf =  s.InvConnectSign[i][j]
				if(sf == 1):
					s.InvConnectPoint[i] = s.InvConnectPoint[i] + [f[0]]
				elif(sf == -1):	
					s.InvConnectPoint[i] = s.InvConnectPoint[i] + [f[1]]


		# set flag for neighbourhood of boundary
		for i in range(s.np):		
			if(s.flagNode[i]>0):
				for j in s.InvConnectPoint[i]:
					if(s.flagNode[j] == 0):
						s.flagNode[j] = -1

				
	def setProperty(s,label,op,commonpar):
		prop = s.chooseProperty(label)
		par = [op] + commonpar		
		randomBox(prop,par)
			
	def set_af_Lf_Vf(s,updateFibres=[]):
		
		if(updateFibres == []):
			s.Lf = np.zeros(s.nf)
			s.af = np.zeros((s.nf,2))
			s.Vf = np.zeros(s.nf)
			updateFibres = np.arange(s.nf)
			
		for i in updateFibres:
			s.af[i,:], s.Lf[i] = get_afib_lf(s.P,s.ElemFib[i])
			s.Vf[i] = s.Af[i]*s.Lf[i]
		

	
	def setNormalAndAbar(s , updateNodes = []):
			
		vacc = np.zeros(2)
				
		if(updateNodes == []):
			s.normal = np.zeros((s.np,2))
			s.Abar = np.zeros(s.np)
			updateNodes = np.arange(s.np)[s.flagNode>0]
			
		for i in updateNodes:
			
			flag = s.flagNode[i]
			
			if(flag>0):
				
				normal_ik , kmax_i = getNormal_kmax(flag)
						
				vacc = 0.0
	
				for k in range(kmax_i):
					Aik = 0.0
					for f in s.InvConnect[i]:
						Aik = Aik + s.Af[f]*np.abs(np.dot(s.af[f],normal_ik[k,:]))
					
					vacc = vacc + Aik*normal_ik[k,:]
								
				s.Abar[i] = np.linalg.norm(vacc)
				s.normal[i,:] = (1.0/s.Abar[i])*vacc
		

	def createNetwork(s,op = 0): # op = 0 is for asymFac and 1 is by thetaEstimate

		#~ points = randomPointsWithBoundary(nf)
		if(op == 0):
			s.nx,s.ny = estimateGridPointsByAsymetry(s.nFibPrevision,s.asymFac)
		elif(op == 1):
			s.nx,s.ny = estimateGridPointsByTheta(s.nFibPrevision,s.meanTheta)
			
		# update estimatives
		s.asymFac = s.nx - s.ny
		s.meanTheta = np.arctan(float(s.nx-1)/float(s.ny-1))
			
		points = structuredPointsWithBoundary(s.nx,s.ny)
	 
		delau = Delaunay(points)
		
		s.P = delau.points # get points	
		s.delauTri = delau.simplices
		s.ElemFib = tri2Edges(s.delauTri) # transform triangles into edges
		
		s.nf = len(s.ElemFib)
		s.np = len(s.P)
		
	# Remove horizontal or vertical fibres. 
	def removeFibres(s,crit,param=[]):
				
		count = 0 
		k = 0 	
		while(k < s.nf):	
			if(crit(s,k,param)):
				del s.ElemFib[k]
				s.nf = s.nf - 1
			else:
				k = k + 1

	def modifyFibreProperty(s,label,foo,param,smoothFactor,binaryLevel = False):
		
		prop = s.chooseProperty(label)
		
		if(binaryLevel): # first parameter is the contrast
			for i in range(s.nf): 
				if(foo(s,i,param)):
					prop[i] = smoothFactor*prop[i]	
		else:
			for i in range(s.nf): 
				fac , inside = foo(s,i,param) 
				if(inside):
					prop[i] = smoothFactor*fac*prop[i]	
	

	def chooseProperty(s,label):
		if(label == 'r0'):
			prop = s.r0f
		elif(label == 'A'):
			prop = s.Af			
		elif(label == 'la'):
			prop = s.laf
		else:
			print 'property has not found'
			return 0
	
		return prop	

	def cleanDangledPoints(s):
		
		newInd = np.zeros(s.np,dtype = 'int')

		k = 0
		for i in range(s.np):
			if(len(s.InvConnect[i]) == 0): 
				newInd[i] = -1
				k = k + 1
			else:
				newInd[i] = i - k
		
		for f in s.ElemFib:
			f[0] = newInd[f[0]]
			f[1] = newInd[f[1]]
		
		for t in s.delauTri:
			t[0] = newInd[t[0]]
			t[1] = newInd[t[1]]
			t[2] = newInd[t[2]]
			
		k = 0
		s.delauTri = list(s.delauTri)
		nt = len(s.delauTri) 
		while(k<nt):
			if(any(s.delauTri[k] == -1)):
				del s.delauTri[k]
				nt = nt - 1
			else:
				k = k + 1
			
		s.P = s.P[newInd>-1,:]
		s.np = len(s.P)
		
	# update the points by adding a random pertubation. For periodic the the average. 
	def addPertubation(s,pertPoint,isPeriodic=False,facx = 1.0): 
		
		r = pertPoint*np.random.rand(s.np)
		theta = 2.0*np.pi*np.random.rand(s.np)
		
		rx = facx*r*np.cos(theta)
		ry = r*np.sin(theta)
				
		s.P[:,0] = s.P[:,0] + rx
		s.P[:,1] = s.P[:,1] + ry
		
		if(isPeriodic):
			for i in range(1,s.nx-1):
				j = 2*s.nx + s.ny - 3 - i # esoteric formula :)
				v = 0.5*(s.P[i,0] + s.P[j,0])
				s.P[i,0] = v
				s.P[j,0] = v
			
			for i in range(s.nx, s.nx + s.ny -2):
				j = 3*s.nx + 2*s.ny - 5 - i # esoteric formula :)
				v = 0.5*(s.P[i,1] + s.P[j,1])
				s.P[i,1] = v
				s.P[j,1] = v
			
			
			 
				
	def correctPoints(s,maxP,minP):
		
		s.P[s.P[:,0]<minP , 0] = minP 
		s.P[s.P[:,0]>maxP , 0] = maxP
		s.P[s.P[:,1]<minP , 1] = minP
		s.P[s.P[:,1]>maxP , 1] = maxP
		
		# 1 = bottom, 2 = right, 3 = top, 4 = left
		# 5 = bottom/left , 6 = bottom/right, 7 = right/top, 8 = top/left
		s.P[s.flagNode == 1, 1] = 0.0
		s.P[s.flagNode == 2, 0] = 1.0
		s.P[s.flagNode == 3, 1] = 1.0
		s.P[s.flagNode == 4, 0] = 0.0
		s.P[s.flagNode == 5, 0] = 0.0
		s.P[s.flagNode == 5, 1] = 0.0
		s.P[s.flagNode == 6, 0] = 1.0
		s.P[s.flagNode == 6, 1] = 0.0
		s.P[s.flagNode == 7, 0] = 1.0
		s.P[s.flagNode == 7, 1] = 1.0		
		s.P[s.flagNode == 8, 0] = 0.0
		s.P[s.flagNode == 8, 1] = 1.0
	
	# omega = 1 , total smooth, omega = 0 , no smooth
	def smooth(s,omega,isAll=True):
		
		listPoints = []

		if(isAll):
			listPoints = np.arange(s.np)
		else: 
			listPoints = np.arange(s.np)[s.flagNode < 1]
			
		np.random.shuffle(listPoints)
		
		
		for i in listPoints:
			pG = np.zeros(2)
			
			
			for j in range(len(s.InvConnect[i])):
				f = s.InvConnect[i][j]
				sf = s.InvConnectSign[i][j]
				if(sf == 1):
					pG = pG + s.P[s.ElemFib[f,0],:]
				elif(sf == -1):
					pG = pG + s.P[s.ElemFib[f,1],:]
			
			
			#~ print i, float(len(s.InvConnect[i])), s.P[i,:]
			pG = pG/float(len(s.InvConnect[i]))
			
			s.P[i,:] = omega*pG + (1.0 - omega)*s.P[i,:]
			
			
		
		

	# Create replication of RVEs where:
	# net0 is the Network of origin
	# Preg0 are points of the regular network of origin
	# dP0 are the fluctuations of the regular network of origin
	# pRep is the number o repetitions in each direction
	def replicateRVE(s,net0,Preg0,dP0,pRep): 
 	
		tol = 1.0e-5
		
		nP0 = len(Preg0)
		
		times = pRep + 1
		fac = 1.0/float(times)
		print 'multiplication factor ' , fac
		
		s.nx = (net0.nx-1)*times + 1
		s.ny = (net0.ny-1)*times + 1
		
		print s.nx, s.ny
		
		#~ s.nFibPrevision = 2*s.nx*s.ny - s.nx - s.ny + 1 
		s.nFibPrevision = 4*(s.nx-1)*(s.ny-1)
		s.asymFac = s.nx-s.ny
		s.meanTheta = np.arctan(float(s.nx-1)/float(s.ny-1))
		
		s.createNetwork()
		s.removeVertHoriFibers()
		
		s.setFlagsAndConnectivity()
		s.Af = np.ones(s.nf)
		s.laf = np.ones(s.nf)
		
		s.set_af_Lf_Vf()
		s.setNormalAndAbar()
		
		print s.nx, s.ny, s.nf , net0.nf
		
		dP = np.zeros((s.np,2))
			
		pClone = np.zeros(s.np).astype('int')
		fClone = np.zeros(s.nf).astype('int')
		
		
		# copying point pertubations
		for i in range(s.np):
			
			z = copy.deepcopy(s.P[i,:])
			
			Oz = choiceOriginQuadrant(z,times)
			z = (z - Oz)*times
			
			for j in range(nP0):
				z0 = Preg0[j,:]
				
				if(np.linalg.norm(z-z0)<tol):
					pClone[i] = j
					break
		
		
		for i in range(s.np):
			dP[i,:] = fac*dP0[pClone[i],:]
		
		s.P = s.P + dP
	
		
		# copying area, activation stretch, etc
		for i in range(s.nf):
			fi = copy.deepcopy(s.ElemFib[i,:])
			
			z = 0.5*(s.P[fi[0],:] + s.P[fi[1],:])

			Oz = choiceOriginQuadrant(z,times)
			z = (z - Oz)*times
			
			for j in range(net0.nf):
				fj = copy.deepcopy(net0.ElemFib[j,:])
				
				z0 = 0.5*(net0.P[fj[0],:] + net0.P[fj[1],:])

				if(np.linalg.norm(z-z0)<tol):
					fClone[i] = j
					break
		

		for i in range(s.nf):
			s.Af[i] = net0.Af[fClone[i]] ## it should be pre multplied by fac as well to guarantee the same volume
			s.laf[i] = net0.laf[fClone[i]]
			

		s.setNormalAndAbar()
			
		s.P = float(times)*s.P

		s.set_af_Lf_Vf()
		
	# Writing Methods
		
	def writeNetwork(s,Ndof,Nsubsteps,Nparam,fignum,opInifile = 'SGP', opAdditionalElements = 0, addAuxNodes = 2, addAuxNodeOnFibres = 1):
		
		Paux = copy.deepcopy(s.P)
		if(addAuxNodes > 0):
			Padd = np.array(addAuxNodes*[[999.9,999.9]])
			Paux = np.concatenate((Paux,Padd),axis=0) # just to be compatible with converter
		
		
		ElemFibAux = copy.deepcopy(s.ElemFib)
		
		if(addAuxNodeOnFibres == 1):
			AuxVertex = (len(Paux)-1)*np.ones(s.nf,dtype = 'int').reshape((s.nf,1))
			ElemFibAux = np.concatenate((ElemFibAux,AuxVertex), axis = 1 )
		
		writeInifileDefault(len(Paux),Ndof,opInifile)
		writeParamFibres(s,Nparam)
		
		
		if(opAdditionalElements == 0): # just global node
			writeMesh(Paux,[ElemFibAux],[[[len(Paux)-1]]],Ndof,Nsubsteps)
		elif(opAdditionalElements == 1): # triangles and global node
			writeMesh(Paux,[ElemFibAux],[s.delauTri,[[len(Paux)-1]]],Ndof,Nsubsteps)
		elif(opAdditionalElements == 2): # torsinal springs (angles) and global node
			writeMesh(Paux,[ElemFibAux],[s.ElemAngle,[[len(Paux)-1]]],Ndof,Nsubsteps)
		elif(opAdditionalElements == 3): # triangles, torsinal springs (angles) and global node
			writeMesh(Paux,[ElemFibAux],[s.delauTri,s.ElemAngle,[[len(Paux)-1]]],Ndof,Nsubsteps)
		elif(opAdditionalElements == 4): # torsinal springs (angles), torsional springs (oposite angles) and global node
			writeMesh(Paux,[ElemFibAux],[s.ElemAngle,s.ElemAngle2,[[len(Paux)-1]]],Ndof,Nsubsteps)
		elif(opAdditionalElements == 5): # torsinal springs (angles), torsional springs (oposite angles), boundary torsional springs and global node
			writeMesh(Paux,[ElemFibAux],[s.ElemAngle,s.ElemAngle2,s.ElemAngle3,[[len(Paux)-1]]],Ndof,Nsubsteps)
		
		
			
		os.system("cp Mesh.txt  Mesh" + str(fignum) + ".txt ")
		os.system("cp IniFile000.txt  IniFile000_" + str(fignum) + ".txt ")
		os.system("cp Param000.txt  Param000_" + str(fignum) + ".txt ")
		os.system("cp Bten.txt  Bten" + str(fignum) + ".txt ")


	def get_normal_bar(s):
		
		vacc = np.zeros(2)
		Acc = 0.0
		
		for i in range(s.np):
			if(s.flagNode[i]>0):
				vacc = vacc + s.Abar[i]*s.normal[i,:]
				Acc = Acc + s.Abar[i]
		
		normal_bar = vacc/Acc
		
		return normal_bar

	def getCentreFibre(s,k):
		f = s.ElemFib[k]	
		x0 = 0.5*(s.P[f[0],0] + s.P[f[1],0])
		y0 = 0.5*(s.P[f[0],1] + s.P[f[1],1])
		
		return x0,y0
	
	def get_yG(s):
		
		yG = np.zeros(2)
		
		for i in range(s.nf):
			f = s.ElemFib[i]
			yG = yG + 0.5*s.Vf[i]*(s.P[f[0],:] + s.P[f[1],:] )
			
		yG = yG/np.sum(s.Vf)
		
		return yG

	def get_Bten(s):
		
		Bten = np.zeros((2,2))
		
		for f in range(s.nf):
			Bten = Bten + s.Vf[f]*np.outer(s.af[f,:],s.af[f,:])
			
		Bten = Bten/np.sum(s.Vf)
		
		return Bten

# ================= auxiliary functions ==========================


# notice that nf = 4*(nx - 1)*(ny - 1) in regular grid excluding vertical and horizontal fibres

# computes actual number of fibres, given a prevision and asymetric factor (nx - ny = asymFac)
def estimateGridPointsByAsymetry(nfPrev,asymfac = 0):
	a = 1.0
	b = float(asymfac)
	c =  - 0.25*float(nfPrev)
	
	delta = np.sqrt(b*b - 4.0*a*c)
	y = (-b + delta)/(2.0*a)
	
	ny = y + 1.0
	nx = ny + float(asymfac)
	
	nx = np.round(nx).astype('int')
	ny = np.round(ny).astype('int')

	print 'estimated nx,ny=',nx, ny
	return nx , ny

# computes actual number of fibres, given a prevision and asymetric factor (nx - ny = asymFac)
def estimateGridPointsByTheta(nfPrev,theta = 0.25*np.pi):
	
	ny = 1.0 + np.sqrt(float(nfPrev)/(4.0*np.tan(theta)))
		
	nx = 1 + np.tan(theta)*(ny-1)
	
	nx = np.round(nx).astype('int')
	ny = np.round(ny).astype('int')

	print 'estimated nx,ny=',nx, ny
	return nx , ny

# create structured points, with boundary structured
# used in createMatrixAndNetwork
def structuredPointsWithBoundary(nx,ny): 

	nf = 2*ny*nx - ny - nx + 1
	npb = 2*(nx + ny - 2 ) 
	
	points = np.zeros((nf,2))
	points[0:npb,:] = borderPoints(0.0,0.0,1.0,1.0,nx,ny)
	
	tx = 1.0/(2.0*(nx-1.0))
	pxa = tx + np.linspace(0.0,1.0,nx)[0:nx-1]
	ty = 1.0/(2.0*(ny-1.0))
	pya = ty + np.linspace(0.0,1.0,ny)[0:ny-1]
	
	pxb = np.linspace(0.0,1.0,nx)[1:nx-1]
	pyb = np.linspace(0.0,1.0,ny)[1:ny-1]

	
	k = npb
	for px in pxa:
		for py in pya:
			points[k,0] = px  
			points[k,1] = py
			k = k + 1 

	for px in pxb:
		for py in pyb:
			points[k,0] = px
			points[k,1] = py
			k = k + 1
			

	return points

def tri2Edges(tri):
	
	e = []
	
	for t in tri:
		ei = []
		for i in range(3):
			if(i==0):
				ei = [t[0],t[1]]
			elif(i==1):
				ei = [t[1],t[2]]
			else:
				ei = [t[2],t[0]]
		
			if(ei[0]>ei[1]):
				ei = [ei[1],ei[0]]
		
			if(not(ei in e)):
				e = e + [ei]
	
	return e

# create structured boudary points
# used in structuredPointsWithBoundary
def borderPoints(a,b,c,d,nx,ny):
	
	x = np.linspace(a,c,nx)
	y = np.linspace(b,d,ny)
	
	N = 2*(nx + ny - 2)  
	
	p = np.zeros((N,2))
	
	k = 0
	for i in range(nx):
		p[k,0] = x[i] 
		p[k,1] = b 
		k = k + 1

	for j in range(1,ny):
		p[k,0] = c 
		p[k,1] = y[j]
		k = k + 1
		
	for i in range(1,nx):
		p[k,0] = x[nx-1-i] 
		p[k,1] = d
		k = k + 1		

	for j in range(1,ny-1):
		p[k,0] = a 
		p[k,1] = y[ny - 1 - j]
		k = k + 1
								
	return p
	 

# In a segment from 0 to 1, assigns 0 for the left extremity, 1 for the middle and 2 for the right extremity
# used in setFlagBound.
def classify(x):
	
	tol = 1.0e-8
	c = 0
	
	if(abs(x)<tol): # left extreme
		c = 0
	elif(abs(x-1.0)<tol): # right extreme
		c = 1
	else: # in the middle
		c = 2 
	
	return c
	

# returns a array [Ozx,Ozy] with the origin of the quadrant
# used in replicateRVE
def choiceOriginQuadrant(z,times):
	fac = 1.0/float(times)
	
	Ozx = 0.0
	Ozy = 0.0
	
	for i in range(times):
		Ozx = (times-i-1)*fac
		if(z[0]>Ozx): 
			break
	
	for j in range(times):
		Ozy = (times-j-1)*fac
		if(z[1]>Ozy): 
			break
	
	return np.array([Ozx,Ozy])


# get afib and lf given the list of positions and a edge
def get_afib_lf(P,e):
	dx = P[e[1],:] - P[e[0],:]
	
	lf = np.sqrt(np.dot(dx,dx))
	
	return dx/lf , lf



def getNormal_kmax(flag):
	
	# for points : -1 = neighbour to boundary, 0 =  very interior, 1 = bottom, 2 = right, 3 = top, 4 = left
	# 5 = bottom/left , 6 = bottom/right, 7 = right/top, 8 = top/left
	kmax = 2
	
	kmax_i = 0
	normal_ik = np.zeros((kmax,2)) 

	if(flag == 1):
		kmax_i = 1
		normal_ik[0,0] = 0.0 ; normal_ik[0,1] = -1.0 	
	elif(flag == 2):
		kmax_i = 1
		normal_ik[0,0] = 1.0 ; normal_ik[0,1] = 0.0
	elif(flag == 3):
		kmax_i = 1
		normal_ik[0,0] = 0.0 ; normal_ik[0,1] = 1.0
	elif(flag == 4):
		kmax_i = 1
		normal_ik[0,0] = -1.0 ; normal_ik[0,1] = 0.0
	elif(flag == 5):
		kmax_i = 2
		normal_ik[0,0] = 0.0 ; normal_ik[0,1] = -1.0
		normal_ik[1,0] = -1.0 ; normal_ik[1,1] = 0.0
	elif(flag == 6):
		kmax_i = 2
		normal_ik[0,0] = 0.0 ; normal_ik[0,1] = -1.0
		normal_ik[1,0] = 1.0 ; normal_ik[1,1] = 0.0
	elif(flag == 7):
		kmax_i = 2
		normal_ik[0,0] = 0.0 ; normal_ik[0,1] = 1.0
		normal_ik[1,0] = 1.0 ; normal_ik[1,1] = 0.0
	elif(flag == 8):	
		kmax_i = 2
		normal_ik[0,0] = 0.0 ; normal_ik[0,1] = 1.0
		normal_ik[1,0] = -1.0 ; normal_ik[1,1] = 0.0
	
	return normal_ik , kmax_i


def randomBox(V,par):
	
	n = len(V)
	#~ V = np.zeros(n)
	
	op = int(par[0])
	v = 0.0
	for i in range(n):
		if(op==1):
			v = par[1]
		elif(op==2): # around 1
			v = par[1] + par[2]*(np.random.rand() - 0.5)
		elif(op==3): # just greater than 1
			r = np.random.rand()
			if(r > par[3]): # commonpar[2] is the threshold
				v = par[1] + par[2]*(r - par[3])
			else:
				v = par[1] 
			
		elif(op==4): # random following normal distribution
			mu = par[1]
			sigma = par[2]
			threshold_min = par[3]
			threshold_max = par[4]
			v = sigma*np.random.randn() + mu
			
			v = max(threshold_min,v)
			v = min(threshold_max,v)
		elif(op==5): # random following gamma distribution
			mu = par[1]
			sigma = par[2]
			threshold_min = par[3]
			threshold_max = par[4]
			
			mu = mu - threshold_min
			
			theta = sigma*sigma/mu
			kappa = (mu/sigma)**2.0
			
			v = threshold_min + np.random.gamma(kappa,theta)
			while (v > threshold_max):
				v = threshold_min + np.random.gamma(kappa,theta)
			
		V[i] = v

	return V

def stripDecrease(net,k,param):
	
	a = param[0]
	b = param[1]
	c = param[2]
	delta = param[3]
	
	x0,y0 = net.getCentreFibre(k)
	
	d = np.abs( (c - a*x0 - b*y0)/np.sqrt(a*a + b*b) )
	
	gamma = 0.6 # later it will be amplified

	if(d<delta):
		smoothFac = gamma + (1.0 - gamma)*(d/delta)**2.0
		return smoothFac, True
		 
	return 1.0 , False 
	

def ballsDecrease(net,k,param):
	n = param[0] # number of balls
	
	smoothFac = 1.0
	
	x0,y0 = net.getCentreFibre(k)
	
	gamma = 0.5 # the default gamma is 0.5
	
	for i in range(n):
		cx = param[i*3 + 1]
		cy = param[i*3 + 2]
		radius = param[i*3 + 3]
		
		dist = np.sqrt((cx - x0)**2.0 + (cy - y0)**2.0) 
		
		if(dist < radius) :
			smoothFac = gamma + (1.0 - gamma)*(dist/radius)**2.0
			return smoothFac, True
		
	return 1.0, False
		
	
def InBallFibresCriterium(net,k,param):	
	x0 = np.zeros(2)
	x0[0] = param[0]
	x0[1] = param[1]
	R = param[2]
	
	px,py = net.getCentreFibre(k)
	
	d = np.sqrt( (px - x0[0])**2.0 + (py - x0[1])**2.0 )
	
	if(d<R):
		return True

	return False
	
def InSquareFibresCriterium(net,k,param):	
	x0 = param[0]
	y0 = param[1]
	Lx = param[2]
	Ly = param[3]
	
	px,py = net.getCentreFibre(k)
	
	if(px>x0 and px < x0 + Lx and py>y0 and py < y0 + Ly):
		return True
		
	return False


def VertHorizFibresCriterium(net,k,param=[]):
	tol = 0.1/float(net.nx*net.ny)
	
	f = net.ElemFib[k]
	dx = net.P[f[0],0] - net.P[f[1],0]
	dy = net.P[f[0],1] - net.P[f[1],1]
			
	if(np.abs(dx*dy) < tol):
		return True

	return False

def DangledFibresCriterium(net,k,param=[]):
	tol = 0.1/float(net.nx*net.ny)
	
	f = net.ElemFib[k]
	If0 = len(net.InvConnect[f[0]])
	If1 = len(net.InvConnect[f[1]])
	Ff0 = net.flagNode[f[0]]
	Ff1 = net.flagNode[f[1]]
	
	if( (If0 == 1  or If1 == 1) and (Ff0 <= 0 and Ff1 <= 0) ):
		return True

	return False
	
	
def RandomFibresCriterium(net,k, param):
	f = net.ElemFib[k]
	Ff0 = net.flagNode[f[0]]
	Ff1 = net.flagNode[f[1]]
	
	p = param[0] # Probability
	r = np.random.rand()
		
	if( Ff0 <= 0 and Ff1 <= 0 and r<p ): # only to internal fibres	
		return True

	return False


def inverseAtractor(P,param):
	tol = 1.0e-8
	x0 = np.zeros(2)
	x0[0] = param[0]
	x0[1] = param[1]
	R = param[2]
	Rmax = param[3]
	
	nP = len(P)
	
	for i in range(nP):
		d = np.linalg.norm(P[i,:] - x0)
		
		if(d<tol):
			P[i,0] = P[i,0] + R
		elif(d<Rmax):
			delta = (R - (R/Rmax)*np.linalg.norm(P[i,:] - x0))/np.linalg.norm(P[i,:] - x0)
			P[i,:] = P[i,:] + delta*(P[i,:] - x0)
	
	return P


def combinedCriterium(crit1,crit2):
	crit = lambda net,k,params: crit1(net,k,params[0]) and crit2(net,k,params[1]) 
	
	return crit
	
def convertToBinaryCriterium(critStar):
	crit = lambda net,k,param: critStar(net,k,param)[1]  
	
	return crit
	
