#~ from generator import *
#~ from delaunay import *
from scipy.spatial import Delaunay
import numpy as np
import matplotlib.pyplot as plt
import copy
import os, sys


class Network:
	
	def __init__(s): # default constructor 
		
		s.altenancyPointArea = [1,1]
		
		s.asymFac = 0 
		s.meanTheta = 0.25*np.pi
		s.nFibPrevision = 20
	
		s.P=[]
		s.ElemFib=[]
		
		s.delauTri = [] # it is optional, used to model the matrix, it improves convergence
		
		s.nx = 0
		s.ny = 0
		s.nf = 0
		s.np = 0
		s.npb = 0
		s.flagNode = []
		s.flagFibres = []
		s.InvConnect = []
		s.InvConnectSign = []
		s.InvConnectPoint = []
			
		s.Lf = []
		s.Af = []
		s.Vf = []
		s.lfa = []
		s.af = []
		s.yrel = []
		s.yG = []
		s.normal_bar = []
		s.Bten = []
		s.BtenInvT = []
	
		s.normal = []
		s.Abar = []
		
	    # methods
		# setFlagBound() 
		# setAf(op, commonpar)
		# set_lfa(op, commonpar) 
		# set_af_Lf_Vf() 
		# createNetwork() : create structured network based on delaunay triangulation (at the moment)
		# removeVertHoriFibers()
		# addPertubation(delta)


	# build inverse connectivity (and its sign (-1 or 1, starting or ending) ) and classify all points and lines as the follwings categories
	# for points : -1 = neighbour to boundary, 0 =  very interior, 1 = bottom, 2 = right, 3 = top, 4 = left
	# 5 = bottom/left , 6 = bottom/right, 7 = right/top, 8 = top/left
	# for lines : 0 = interior fibres, 1 = neighbourhood fibres.
	def setFlagsAndConnectivity(s):

		s.InvConnect = s.np*[[]]
		s.InvConnectSign = s.np*[[]]
		s.InvConnectPoint = s.np*[[]]
		
		
		for i in range(s.nf):
			s.InvConnect[s.ElemFib[i,0]] = s.InvConnect[s.ElemFib[i,0]] + [i]
			s.InvConnect[s.ElemFib[i,1]] = s.InvConnect[s.ElemFib[i,1]] + [i]
			s.InvConnectSign[s.ElemFib[i,0]] = s.InvConnectSign[s.ElemFib[i,0]] + [-1]
			s.InvConnectSign[s.ElemFib[i,1]] = s.InvConnectSign[s.ElemFib[i,1]] + [1]						

		s.flagNode = np.zeros(s.np,dtype = 'int')
		s.flagFibres = np.zeros(s.nf,dtype = 'int')
		
		flagMatAux = np.array([[5,6,1],[8,7,3],[4,2,0]])
		
		for i in range(s.np):
			pintx = classify(s.P[i,0])
			pinty = classify(s.P[i,1])
			s.flagNode[i] = flagMatAux[pinty,pintx]

			for j in range(len(s.InvConnect[i])):
				f = s.InvConnect[i][j]
				s.flagFibres[f] = 1
				sf =  s.InvConnectSign[i][j]
				if(sf == 1):
					s.InvConnectPoint[i] = s.InvConnectPoint[i] + [s.ElemFib[f,0]]
				elif(sf == -1):	
					s.InvConnectPoint[i] = s.InvConnectPoint[i] + [s.ElemFib[f,1]]


		# set flag for neighbourhood of boundary
		for i in range(s.np):		
			if(s.flagNode[i]>0):
				for j in s.InvConnectPoint[i]:
					if(s.flagNode[j] == 0):
						s.flagNode[j] = -1
						
	def setAf(s,op,commonpar):
		par = [op] + commonpar		
		s.Af = randomBox(s.nf,par)
		
	def set_lfa(s,op,commonpar):
		par = [op] + commonpar		
		s.lfa = randomBox(s.nf,par)
			
			
	def set_af_Lf_Vf(s,updateFibres=[]):
		
		if(updateFibres == []):
			s.Lf = np.zeros(s.nf)
			s.af = np.zeros((s.nf,2))
			s.Vf = np.zeros(s.nf)
			updateFibres = np.arange(s.nf)
			
		for i in updateFibres:
			s.af[i,:], s.Lf[i] = get_afib_lf(s.P,s.ElemFib[i,:])
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

	
	def rotateNetwork(s,thetaRot):
		
		Qrot = np.zeros((2,2))
		Qrot[0,0] = np.cos(thetaRot) 
		Qrot[0,1] = -np.sin(thetaRot) 
		Qrot[1,0] = np.sin(thetaRot) 
		Qrot[1,1] = np.cos(thetaRot) 
		
		for i in range(s.np):
			s.P[i,:] = np.dot(Qrot,s.P[i,:])
			s.normal[i,:] = np.dot(Qrot,s.normal[i,:])
			
		for i in range(s.nf):
			s.af[i,:] = np.dot(Qrot,s.af[i,:])
		

	def adjustDeterministicBoundArea(s):
		
		boundNodes = np.arange(s.np)[s.flagNode>0]
		
		Scx = 0.0
		Scy = 0.0
		Si = np.zeros(4)
		
		for i in boundNodes:
			
			flag = s.flagNode[i]
			
			if(flag > 4):
				Scx = Scx + s.Abar[i]*s.normal[i,0]
				Scy = Scy + s.Abar[i]*s.normal[i,1]
				
			else:
				Si[flag-1] = Si[flag-1] + s.Abar[i]
		
		
		delta_x = -(Si[1] - Si[3] + Scx)/(Si[1] + Si[3])
		delta_y = -(Si[2] - Si[0] + Scy)/(Si[0] + Si[2])
		
		
	
		for i in boundNodes:
			
			flag = s.flagNode[i]
			
			for f in s.InvConnect[i]:
				if(flag == 1): 
					s.Af[f] = s.Af[f]*(1.0 - delta_y)					
				elif(flag == 2): 
					s.Af[f] = s.Af[f]*(1.0 + delta_x)
				elif(flag == 3):
					s.Af[f] = s.Af[f]*(1.0 + delta_y)		
				elif(flag == 4):
					s.Af[f] = s.Af[f]*(1.0 - delta_x)
		
		
		s.setNormalAndAbar()
		
		print(delta_x , delta_y, functionalNBCnormal(s))	
		

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
		s.npb = 2*(s.nx + s.ny - 2 )
		#~ 2*(s.nx + s.ny - 2 ) ## valid only for structured, == len s.PtsBound??
		
	# Remove horizontal or vertical fibres. 
	# used in the main program
	# todo: generalise it to remove more general patterns
	def removeVertHoriFibers(s):

		tol = 0.1/float(s.nx*s.ny)
		
		ElemFibNew = []
			
		for f in s.ElemFib:
			dx = s.P[f[0],0] - s.P[f[1],0]
			dy = s.P[f[0],1] - s.P[f[1],1]
			
			if(np.abs(dx*dy) > tol):
				ElemFibNew = ElemFibNew + [f]

		s.ElemFib = np.array(ElemFibNew)
		s.nf = len(s.ElemFib)
		
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
		print('multiplication factor ' + str(fac))
		
		s.nx = (net0.nx-1)*times + 1
		s.ny = (net0.ny-1)*times + 1
		
		print(s.nx, s.ny)
		
		#~ s.nFibPrevision = 2*s.nx*s.ny - s.nx - s.ny + 1 
		s.nFibPrevision = 4*(s.nx-1)*(s.ny-1)
		s.asymFac = s.nx-s.ny
		s.meanTheta = np.arctan(float(s.nx-1)/float(s.ny-1))
		
		s.createNetwork()
		s.removeVertHoriFibers()
		
		s.setFlagsAndConnectivity()
		s.Af = np.ones(s.nf)
		s.lfa = np.ones(s.nf)
		
		s.set_af_Lf_Vf()
		s.setNormalAndAbar()
		
		print(s.nx, s.ny, s.nf , net0.nf)
		
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
			s.lfa[i] = net0.lfa[fClone[i]]
			

		s.setNormalAndAbar()
			
		s.P = float(times)*s.P

		s.set_af_Lf_Vf()
		
	# Writing Methods
		
	def writeNetwork(s,Ndof,Nsubsteps,Nparam,fignum,opInifile = 'SGP', opIncludeTri = 0, addAuxNodes = 2, addAuxNodeOnFibres = 1):
		
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
		
		
		if(opIncludeTri == 0):
			writeMesh(Paux,[ElemFibAux],[[[len(Paux)-1]]],Ndof,Nsubsteps)
		elif(opIncludeTri == 1):
			writeMesh(Paux,[ElemFibAux],[s.delauTri,[[len(Paux)-1]]],Ndof,Nsubsteps)
			
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
		
	
	def get_yG(s):
		
		yG = np.zeros(2)
		
		for f in range(s.nf):
			yG = yG + 0.5*s.Vf[f]*(s.P[s.ElemFib[f,0],:] + s.P[s.ElemFib[f,1],:] )
			
		yG = yG/np.sum(s.Vf)
		
		return yG

	def get_Bten(s):
		
		Bten = np.zeros((2,2))
		
		for f in range(s.nf):
			Bten = Bten + s.Vf[f]*np.outer(s.af[f,:],s.af[f,:])
			
		Bten = Bten/np.sum(s.Vf)
		
		return Bten

	def isPointPertubation(s,k):
	
		total = s.altenancyPointArea[0] + s.altenancyPointArea[1]
		
		mod = k%total

		return ( mod <  s.altenancyPointArea[0] )

	# Optimise from the perspective of the generic NBC.
	# todo: Create a general algorithm to optimise independently of the cost functional.
	# ParamOpt arguments:
	# maxit => maximum number of nodes mouvements
	# maxnochange => maximum number of mouvements without minimization, restart amp to amp0
	# amp0 => initial maximum amplitude of mouvement
	# alpha => ratio of decreasing o amp, as amp = alpha*amp0
	# dJmin => minimum decrease in each step
	# Jtol => Tolerance of minimisation to quit the loop before reach maxit
	# maxP => maximum coordinate allowed to P (in x or in y)
	# minP => idem of maxP, but it is the minimum
	def optimize(s,paramOpt,funcCost):
	
		J0 = funcCost(s)
		
		# initialises parameters
		maxit, maxnochange, dJmin, Jtol, ampP0, alphaP, gammaP , ampA0, alphaA, maxA, minA, pertPoint, alphaPert, restartPert, omegaSmooth, timesSmooth, Jtol2 = paramOpt
		knochange = 0
		k = 0
		kpert = 0
		
		ampP = copy.deepcopy(ampP0)
		
		Afmean = np.mean(s.Af)
		ampA0 = ampA0*Afmean
		ampA = copy.deepcopy(ampA0)
		maxA = maxA*Afmean
		minA = minA*Afmean
		
		maxP = gammaP*1.0 + (1.0 - gammaP)*np.max(s.P[s.flagNode<0,:].flatten())
		minP = (1.0 - gammaP)*np.min(s.P[s.flagNode<0,:].flatten())
		
		pertPoint0 = pertPoint*np.mean(s.Lf)
		pertPoint = copy.deepcopy(pertPoint0) 
		
		neighBoundPts = np.arange(s.np)[s.flagNode == -1]
		nnbp = len(neighBoundPts)
		np.random.shuffle(neighBoundPts)
		
		neighBoundFibres = np.arange(s.nf)[s.flagFibres == 1]
		nnbf = len(neighBoundFibres)
		np.random.shuffle(neighBoundFibres)
		
		kp = 0
		kf = 0
		ksmooth = 0
		
		pointPertubation = False
		# main loop
		while(k<maxit and J0>Jtol):
		
			status = True
			
			pointPertubation = s.isPointPertubation(k)
			
			if(pointPertubation): 
				#~ print 'pointPertubation' , kp
				# selects randomly the points on the bonundary neighbourhood 
				ipp = neighBoundPts[kp%nnbp]				
				kp = kp + 1
					
				dp, J = perturbPoint(funcCost,s,ampP,maxP,minP,ipp)
				
				status, J0, ampP, knochange = updatePertubation(J0, ampP, knochange, J, ampP0, alphaP, maxnochange, dJmin)
				
				if(status):
					s.P[ipp,:] = s.P[ipp,:] + dp	
					s.set_af_Lf_Vf(s.InvConnect[ipp])
					s.setNormalAndAbar(s.InvConnectPoint[ipp])
					
			else:
				#~ print 'areaPertubation' , kf
				ipp = neighBoundFibres[kf%nnbp]
				kf = kf + 1
								
				dA, J = perturbArea(funcCost,s,ampA,maxA,minA,ipp)
						
				status, J0, ampA, knochange = updatePertubation(J0, ampA, knochange, J, ampA0, alphaA, maxnochange, dJmin)
				
				if(status):
					s.Af[ipp] = s.Af[ipp] + dA	
				
			if(J0<Jtol2 and status and ksmooth < timesSmooth):
				print('smooth') 
				s.smooth(omegaSmooth,True)
				s.correctPoints(maxP,minP)
				s.set_af_Lf_Vf()	
				s.setNormalAndAbar()
				ksmooth = ksmooth + 1			
			
			if(not status and knochange == 0):
				print('pertubation') 
				s.addPertubation(pertPoint)
				#~ s.smooth(omegaSmooth,timesSmooth,True)
				s.correctPoints(maxP,minP)
				s.set_af_Lf_Vf()
				s.setNormalAndAbar()
				J0 = funcCost(s)
				if(kpert < restartPert):
					pertPoint = alphaPert*pertPoint
					kpert = kpert + 1
				else:
					pertPoint = copy.deepcopy(pertPoint0)
					kpert = 0
					
			k = k + 1
			
			print(J0)
			#~ print s.Abar
			
		#~ print J0



# ================= auxiliary functions ==========================


def updatePertubation(J0, amp, knochange, J, amp0, alpha, maxnochange, dJmin): 
	# test if the minimum is worthy
	if(J < J0):
		status = True
		
		if((J0-J)/J0 < dJmin):
			print(" (J0-J)/J0 < dJmin ")
			amp = amp0
			
		J0 = J
		knochange = 0
		
	else:
		status = False
		amp = alpha*amp
		
		knochange = knochange + 1
		if(knochange > maxnochange):
			print(" knochange > maxnochange ")
			amp = amp0
			knochange = 0
		
	
	return status, J0, amp, knochange


# computes functional for old BNC (\sum_{i\in setNodeGamma}  \sum_f [f,i] A_f a_f ) 
# used in the optimization algorithm
def functionalNBCfib(net): 
	
	vacc = np.zeros(2) 

	for i in range(net.npb):
		for j in range(len(net.InvConnect[i])):
			f = net.InvConnect[i][j]
			sfi = net.InvConnectSign[i][j]
			
			vacc = vacc + sfi*net.Af[f]*net.af[f]
			
	return np.sqrt(np.dot(vacc,vacc))
	
# computes functional for new BNC ( (\sum_{i\in setNodeGamma}  \bar{A}_i n_i ) ), todo
# used in the optimization algorithm
def functionalNBCnormal(net): 	
	vacc = np.zeros(2) 

	for i in range(net.np):
		if(net.flagNode[i]>0):
			vacc = vacc + net.Abar[i]*net.normal[i,:]
	
	print(vacc)
			
	return np.sqrt(np.dot(vacc,vacc))


def functionalNBCBoth(net): 	
	return np.sqrt(functionalNBCnormal(net)**2.0 + functionalNBCfib(net)**2.0) 




def perturbPoint(funcCost,net,amp,maxP,minP,ipp):
		
		# selects the random pertubation
	
		ampEfect = amp*np.random.rand(1)
		thet = 2.0*np.pi*np.random.rand(1)
		r = ampEfect*np.array([np.cos(thet),np.sin(thet)])
		
		# get pertubation in all possible directions
		p = np.zeros((8,2))
		p[0,0] = r[0] ; p[0,1] = r[1]; 
		p[1,0] = -r[0] ; p[1,1] = -r[1];
		p[2,0] = -r[0] ; p[2,1] = r[1];
		p[3,0] = r[0] ; p[3,1] = -r[1];
		p[4,0] = r[1] ; p[4,1] = r[0];
		p[5,0] = -r[1] ; p[5,1] = -r[0];
		p[6,0] = -r[1] ; p[6,1] = r[0];
		p[7,0] = r[1] ; p[7,1] = -r[0];
		
		# computes the corresponding cost functionals
		J = np.zeros(8)
		for j in range(8):
	
			net.P[ipp,:] = net.P[ipp,:] + p[j,:]
			net.set_af_Lf_Vf(net.InvConnect[ipp])
			net.setNormalAndAbar(net.InvConnectPoint[ipp])
			
			
			if( max(net.P[ipp,0],net.P[ipp,1]) > maxP or min(net.P[ipp,0],net.P[ipp,1]) < minP ): # penalises pertubations outside a range 
				J[j] = 999.9
			else:
				J[j] = funcCost(net)
			
			net.P[ipp,:] = net.P[ipp,:] - p[j,:]
			net.set_af_Lf_Vf(net.InvConnect[ipp])
			net.setNormalAndAbar(net.InvConnectPoint[ipp])
		
		
		# select the minimum
		ipmin = np.argmin(J)
		J = np.min(J)
		
		return p[ipmin,:] , J

def perturbArea(funcCost,net,amp,maxA,minA,ipp):
		
		# selects the random pertubation
		dA = amp*np.random.randn(1)
		
		Jp = 0.0
		Jm = 0.0
		
		#~ print 'delta area = ', dA
		
		
		net.Af[ipp] = net.Af[ipp] + dA
		
		net.setNormalAndAbar(net.ElemFib[ipp])
		
		if( net.Af[ipp] > maxA ):
			Jp = 999.9
		else: 
			Jp = funcCost(net)
		

		net.Af[ipp] = net.Af[ipp] - 2.0*dA
		net.setNormalAndAbar(net.ElemFib[ipp])
		if( net.Af[ipp] < minA ):
			Jm = 999.9
		else: 
			Jm = funcCost(net)
		
		net.Af[ipp] = net.Af[ipp] + dA
		net.setNormalAndAbar(net.ElemFib[ipp])
			
		if(Jp<Jm):
			return dA, Jp 
		else:
			return -dA, Jm


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

	print('estimated nx,ny=',nx, ny)
	return nx , ny

# computes actual number of fibres, given a prevision and asymetric factor (nx - ny = asymFac)
def estimateGridPointsByTheta(nfPrev,theta = 0.25*np.pi):
	
	ny = 1.0 + np.sqrt(float(nfPrev)/(4.0*np.tan(theta)))
		
	nx = 1 + np.tan(theta)*(ny-1)
	
	nx = np.round(nx).astype('int')
	ny = np.round(ny).astype('int')

	print('estimated nx,ny=',nx, ny)
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
	 
# get points and edges and convert into points on boundary, interior, fibres on boundary and interior, add as well dummy points
# used in createMatrixAndNetwork
def convertGeneric(p,e): 
	
	n = p.shape[0]
		
	xmin = np.min(p[:,0])
	ymin = np.min(p[:,1])
	xmax = np.max(p[:,0])
	ymax = np.max(p[:,1])
	
	p[:,0] = (p[:,0] - xmin)/(xmax - xmin)
	p[:,1] = (p[:,1] - ymin)/(ymax - ymin)

	ElemFib = []

	for ei in e:
		ElemFib = ElemFib + [ei]


	ElemFib = np.array(ElemFib)
	
	return p,ElemFib

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


def writeFigNetwork(net,c=0,lw = 7.0,figNum = 1,filename = 'network.pdf'):

	colorPallette = []
	
	if(c==0):
		nf = len(net.ElemFib)
		m = 0.8
		maxAf = np.max(net.Af)
		minAf = np.min(net.Af)
		
		v = np.zeros(nf)
		
		if(maxAf>minAf):
			v = (m/(maxAf - minAf))*(net.Af - minAf)
			
		for f in range(nf):
			colorPallette = colorPallette + [(v[f],v[f],v[f])]

	k = 0
	plt.figure(figNum,(6,6))
	for f in net.ElemFib:
		x = [net.P[f[0],0],net.P[f[1],0]]
		y = [net.P[f[0],1],net.P[f[1],1]]
	
		if(c == 0):
			plt.plot(x,y, linewidth = lw, color = colorPallette[k])
			k = k + 1
		else:
			plt.plot(x,y, linewidth = lw, color = c)

	plt.tight_layout()
	
	plt.savefig(filename) 


def writeParamFibres(net,Nparam,fignum = -1):
		
	print(Nparam)
	Param = np.zeros((net.nf,Nparam))
	
	# fortran to python convention
	Ipos_flag1 = 0
	Ipos_flag2 = 1 	
	Ipos_Lf = 2
	Ipos_Areaf = 3  
	Ipos_Vf = 4  
	Ipos_lfa = 5 
	Ipos_af = 6
	Ipos_yrel1 = 8
	Ipos_yrel2 = 10
	Ipos_Abar1 = 12
	Ipos_Abar2 = 13
	Ipos_normal1 = 14
	Ipos_normal2 = 16
	Ipos_nConnectFibre1 = 18
	Ipos_nConnectFibre2 = 19
	#~ Ipos_r0 = 20 # new compability
	# Ipos_r0 = 24 # old compability
	 
	
	# just to update
	net.set_af_Lf_Vf()
	net.setNormalAndAbar()
	
	normal_bar = net.get_normal_bar()
	yG = net.get_yG()
	Bten = net.get_Bten()
	BtenInv = np.linalg.inv(Bten)
	
	BtenFile=open("Bten.txt",'w')
	
	BtenFile.write( str(Bten[0,0]) + ' ' + str(Bten[0,1]) + ' ' + str(Bten[1,0]) + ' ' + str(Bten[1,1]) + '\n')
	BtenFile.write( str(BtenInv[0,0]) + ' ' + str(BtenInv[0,1]) + ' ' + str(BtenInv[1,0]) + ' ' + str(BtenInv[1,1]) + '\n')
	BtenFile.write( str(np.sum(net.Vf)) )


	for f in range(net.nf):
		p1 = net.ElemFib[f,0]
		p2 = net.ElemFib[f,1]
		
		Param[f,Ipos_flag1] = net.flagNode[p1]
		Param[f,Ipos_flag2] = net.flagNode[p2]
		Param[f,Ipos_Lf] = 	net.Lf[f]
		Param[f,Ipos_Areaf] = net.Af[f]
		Param[f,Ipos_Vf] = net.Vf[f]
		Param[f,Ipos_lfa] = net.lfa[f]
		Param[f,Ipos_af:Ipos_af+2] = net.af[f,:]
		Param[f,Ipos_yrel1:Ipos_yrel1 + 2] = net.P[p1,:] - yG
		Param[f,Ipos_yrel2:Ipos_yrel2 + 2] = net.P[p2,:] - yG 
		Param[f,Ipos_Abar1] = net.Abar[p1]
		Param[f,Ipos_Abar2] = net.Abar[p2] 
		Param[f,Ipos_normal1:Ipos_normal1 + 2] = net.normal[p1,:] - normal_bar
		Param[f,Ipos_normal2:Ipos_normal2 + 2] = net.normal[p2,:] - normal_bar
		Param[f,Ipos_nConnectFibre1] = len(net.InvConnect[p1])
		Param[f,Ipos_nConnectFibre2] = len(net.InvConnect[p2])		
		# Param[f,Ipos_r0] = net.r0[f]		
		
		
	writeParamCustom(Param)
	
	
	if(fignum != -1): 
		os.system("cp Param000.txt  Param000_" + str(fignum) + ".txt ")



def writeMesh(X,Elem,auxElem,Ndof,Nsubstep):
	
	Ncoord=len(X)
	ElemTotal = Elem + auxElem
	NelemGroups= len(ElemTotal)
	

	mesh=open("Mesh.txt",'w')
	mesh.write("*NODAL DOFs\n")
	mesh.write(str(Ndof) + "\n")
	mesh.write("*DIMEN\n")
	mesh.write("3\n")
	mesh.write("*COORDINATES\n")
	mesh.write(str(Ncoord) + "\n\n")
	
	for i in range(Ncoord):
		mesh.write(format(X[i,0],'.10e') + " " + format(X[i,1],'.10e') + " 0.0 \n")
		
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


def randomBox(n,par):

	V = np.zeros(n)
	
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
	
def stripDecrease(x0,y0,param):
	
	a = param[0]
	b = param[1]
	c = param[2]
	gamma = param[3]
	delta = param[4]
	
	d = np.abs( (c - a*x0 - b*y0)/np.sqrt(a*a + b*b) )
	
	smoothFac = 1.0 # default value
	
	if(d<delta):
		if(gamma < 0.0): # case of fixed decrease
			smoothFac = - gamma
		else:
			smoothFac = gamma + (1.0 - gamma)*(d/delta)**2.0

	return smoothFac

def ballsDecrease(x0,y0,param):
	n = param[0] # number of balls
	
	smoothFac = 1.0
	
	for i in range(n):
		cx = param[i*4 + 1]
		cy = param[i*4 + 2]
		radius = param[i*4 + 3]
		gamma = param[i*4 + 4]
		
		dist = np.sqrt((cx - x0)**2.0 + (cy - y0)**2.0) 
		
		if(dist < radius) :
			if(gamma < 0.0): # case of fixed decrease
				smoothFac = - gamma
			else:
				smoothFac = gamma + (1.0 - gamma)*(dist/radius)**2.0
		
	return smoothFac
		
def decreaseProperty(prop,foo,net,param):
	
	Prop = copy.deepcopy(prop)
	
	for f in range(net.nf):
		p1 = net.ElemFib[f,0]
		p2 = net.ElemFib[f,1]
		
		x0 = 0.5*(net.P[p1,0] + net.P[p2,0])
		y0 = 0.5*(net.P[p1,1] + net.P[p2,1])
		
		Prop[f] = foo(x0,y0,param)*prop[f]

	return Prop

def deleteFibresInBall(ElemFib,P,param):	
	x0 = np.zeros(2)
	x0[0] = param[0]
	x0[1] = param[1]
	R = param[2]
	
	nf = len(ElemFib)
	
	ElemFibNew = [] 
	
	for f in range(nf):
		p1 = ElemFib[f,0]
		p2 = ElemFib[f,1]
		
		px = 0.5*(P[p1,0] + P[p2,0])
		py = 0.5*(P[p1,1] + P[p2,1])
		
		d = np.sqrt( (px - x0[0])**2.0 + (py - x0[1])**2.0 )
		
		if(d>R):
			ElemFibNew = ElemFibNew + [[p1,p2]] 
			
	
	ElemFibNew = np.array(ElemFibNew,dtype = 'int')
	nfNew = len(ElemFibNew)
	
	return ElemFibNew, nfNew
	

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
	
	
	
	
	
	
