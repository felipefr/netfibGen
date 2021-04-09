#~ from generator import *
#~ from delaunay import *
from scipy.spatial import Delaunay
import numpy as np
import matplotlib.pyplot as plt
import copy
import os, sys
from IOfibres import *
import itertools
from funcAux import *
from scipy.spatial import Voronoi
from pyvoro import *

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
		
		s.delauSimplices = [] # it is optional, used to model the matrix, it improves convergence
		
		s.Lx = 1.0
		s.Ly = 1.0
		s.Lz = 0.0
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
		if(s.ndim == 2):
			flagMatAux = np.array([[5,6,1],[8,7,3],[4,2,0]])		
			s.flagNode = map( lambda x: flagMatAux[classify(x[1],s.Ly),classify(x[0],s.Lx)], s.P ) # weird convention  
		elif(s.ndim==3):
			s.flagNode = map( lambda x: (2-classify(x[0],s.Lx)) + (2-classify(x[1],s.Ly))*3 + (2-classify(x[2],s.Lz))*9, s.P )  
		
		s.flagNode = np.array(s.flagNode,dtype = 'int')

	def createAngleElements(s):
		
		thetaCrit = 120.0 # more than this value desconsider the fibre pair to the angle
		
		for i in range(s.np):
			
			for j in range(len(s.InvConnect[i])):
				fj = s.InvConnect[i][j]
				
				for k in range(j+1,len(s.InvConnect[i])): # jump the same fibre and eliminate repetition in combination
					fk = s.InvConnect[i][k]
					
					sjk = s.InvConnectSign[i][j]*s.InvConnectSign[i][k] # changes the sign accordingly to get the right angle
					
					cosTheta = np.dot(sjk*s.af[fj],s.af[fk])
					# in case of some numerical error
					if(cosTheta<-1.0):
						cosTheta = -1.0
					elif(cosTheta>1.0):
						cosTheta = 1.0
					
					theta = np.rad2deg(np.arccos(cosTheta)) 
					
					if(theta<thetaCrit):
						s.ElemAngle = s.ElemAngle + [[ i, s.InvConnectPoint[i][j], s.InvConnectPoint[i][k] ]]	

	# build inverse connectivity (and its sign (-1 or 1, starting or ending) ) and modify classification of interior points points as
	# 0 = interior points, -1 = neighbourhood points.
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
			s.af = np.zeros((s.nf,s.ndim))
			s.Vf = np.zeros(s.nf)
			updateFibres = np.arange(s.nf)
		
		for i in updateFibres:
			s.af[i,:], s.Lf[i] = get_afib_lf(s.P,s.ElemFib[i])
			s.Vf[i] = s.Af[i]*s.Lf[i]

	def setNormalAndAbar(s , updateNodes = []):
		print 'called setNormal'
		
		if(s.ndim==2):
			foo = getNormal_kmax2D
		elif(s.ndim==3):
			foo = getNormal_kmax3D
		
		vacc = np.zeros(s.ndim)
				
		if(updateNodes == []):
			s.normal = np.zeros((s.np,s.ndim))
			s.Abar = np.zeros(s.np)
			updateNodes = np.arange(s.np)[s.flagNode>0]
		
		
		#~ print updateNodes
		
		#~ print 'number of Dangled = ', s.getNumberDangledPoints()
		
		for i in updateNodes:
			
			flag = s.flagNode[i]
			
			if(flag>0):
				
				normal_ik , kmax_i = foo(flag,s.ndim)
			
				#~ print 'a' , i, normal_ik, kmax_i
				vacc = 0.0
	
				for k in range(kmax_i):
					Aik = 0.0
					#~ print 'b', k, s.InvConnect[i]
					for f in s.InvConnect[i]:
						#~ print 'c', f, s.af[f], s.Af[f]
						Aik = Aik + s.Af[f]*np.abs(np.dot(s.af[f],normal_ik[k,:]))
					
					#~ print 'd', Aik
					vacc = vacc + Aik*normal_ik[k,:]
			
				s.Abar[i] = np.linalg.norm(vacc)
				s.normal[i,:] = (1.0/s.Abar[i])*vacc
		
				
		
	def createNetwork(s, param, op = 0): # op = 0 is for asymFac and 1 is by thetaEstimate
			
		if(op<0):
			s.P = np.loadtxt(param[0])
			s.ElemFib = np.loadtxt(param[1],dtype = 'int')
			s.ndim = len(s.P[0])
			s.np = len(s.P)
			s.nf = len(s.ElemFib)
			# supposing 3D
			s.Lx = np.max(s.P[:,0])
			s.Ly = np.max(s.P[:,1])
			s.Lz = np.max(s.P[:,2])
			
			s.setFlagNodes()
			s.setFlagsAndConnectivity()		
			
			s.cleanDangledPoints()
			
			
		elif(s.ndim == 2 and op<2):
			s.createNetwork2D(param,op)
		elif(s.ndim == 2 and op == 2):
			#~ s.createNetworkVoronoi2D(param)
			s.createNetworkVoronoi2D_alternative(param)
		elif(s.ndim == 3 and op == 0):
			s.createNetwork3D(param)
		elif(s.ndim == 3 and op > 0 ):
			#~ s.createNetworkVoronoi3D(param)
			s.createNetworkVoronoi3D_alternative2(param)
			
		s.nf = len(s.ElemFib)
		s.np = len(s.P)
		
		s.setFlagNodes()
		s.setFlagsAndConnectivity()

	def createNetworkVoronoi2D(s,param): 

		nPoints = param[0]
		minDist = param[1]
		maxDist = param[2]
		nx = param[3]
		ny = param[4]
		
		points = generatePointsMinDist(nPoints,minDist,s.ndim)
			
		vor = Voronoi(points)

		s.P = vor.vertices
		e = vor.ridge_vertices
	
		s.ElemFib = cleanEdgesToInfinity(e)
		
		s.nf = len(s.ElemFib)
		s.np = len(s.P)
		
		s.setFlagNodes()
		s.setFlagsAndConnectivity()	
				
		s.cleanVoronoiPoints(maxDist)
		
		ptBound = borderPointsVoronoi2D(0.0,0.0,1.0,1.0,nx,ny)
		
		newFibres = connectInteriorWithBoundaryNodes(ptBound,s.P)
		s.ElemFib = s.ElemFib + newFibres
		s.P = np.concatenate((s.P,ptBound))
		
		s.Lx = np.max(s.P[:,0])
		s.Ly = np.max(s.P[:,1])
		
	def createNetworkVoronoi2D_alternative(s,param): 

		nPoints = param[0]
		minDist = param[1]
		maxDist = param[2]
		nx = param[3]
		ny = param[4]
		
		points = generatePointsMinDist(nPoints,minDist,s.ndim)
		ptBound = borderPointsVoronoi2D(0.0,0.0,1.0,1.0,nx,ny)
		points = np.concatenate((points,ptBound))
		
		vor = Voronoi(points)

		s.P = vor.points
		e = list(vor.ridge_points)
		e = map(list,e)
		
		s.ElemFib = copy.deepcopy(e)
		
		s.nf = len(s.ElemFib)
		s.np = len(s.P)
		
		s.setFlagNodes()
		s.setFlagsAndConnectivity()	
		
		s.Lx = np.max(s.P[:,0])
		s.Ly = np.max(s.P[:,1])
		
		
	def createNetworkVoronoi3D(s,param): 

		nPoints = param[0]
		minDist = param[1]
		maxDist = param[2]
		nx = param[3]
		ny = param[4]
		nz = param[5]
		
		points = generatePointsMinDist(nPoints,minDist,s.ndim)
			
		vor = Voronoi(points)

		s.P = vor.vertices
		e = list(vor.ridge_vertices)
		#~ e = list(vor.regions)
		e = map(list,e)
		e = cleanEdgesToInfinity(e)
		
		e = simplices2Edges(e)
		
		s.ElemFib = copy.deepcopy(e)
		
		s.nf = len(s.ElemFib)
		s.np = len(s.P)
		
		s.setFlagNodes()
		s.setFlagsAndConnectivity()	
				
		s.cleanVoronoiPoints(maxDist)
		
		ptBound = borderPointsVoronoi3D(0.0,0.0,0.0,1.0,1.0,1.0,nx,ny,nz)
		
		newFibres = connectInteriorWithBoundaryNodes(ptBound,s.P)
		s.ElemFib = s.ElemFib + newFibres
		s.P = np.concatenate((s.P,ptBound))
		
		s.Lx = np.max(s.P[:,0])
		s.Ly = np.max(s.P[:,1])
		s.Lz = np.max(s.P[:,2])
		
	def createNetworkVoronoi3D_alternative(s,param): 
		
		nPoints = param[0]
		minDist = param[1]
		maxDist = param[2]
		nx = param[3]
		ny = param[4]
		nz = param[5]
		
		
		points = generatePointsMinDist(nPoints,minDist,s.ndim)
		ptBound = borderPointsVoronoi3D(0.0,0.0,0.0,1.0,1.0,1.0,nx,ny,nz)
		points = np.concatenate((points,ptBound))
		
		vor = Voronoi(points)

		s.P = vor.points
		e = list(vor.ridge_points)
		e = map(list,e)
		
		s.ElemFib = copy.deepcopy(e)
		
		s.nf = len(s.ElemFib)
		s.np = len(s.P)
		
		s.setFlagNodes()
		s.setFlagsAndConnectivity()	
		
		s.Lx = np.max(s.P[:,0])
		s.Ly = np.max(s.P[:,1])
		s.Lz = np.max(s.P[:,2])
		
	def createNetworkVoronoi3D_alternative2(s,param): 
		
		nPoints = param[0]
		minDist = param[1]
		maxDist = param[2]
		nx = param[3]
		ny = param[4]
		nz = param[5]
		
		points = generatePointsMinDist(nPoints,minDist,s.ndim)
		
		voro = compute_voronoi(points, [[0.0,1.0],[0.0,1.0],[0.0,1.0]], 0.5)
		
		s.P, s.ElemFib = convertVoroPP2fibres(voro)

		s.nf = len(s.ElemFib)
		#~ s.ElemFib = s.ElemFib + recombineEdges(s.ElemFib,s.nf/2)
		#~ s.ElemFib = uniqueEdges(s.ElemFib)
		
		#~ s.eliminateFibresOverBoundary()
		
		s.nf = len(s.ElemFib)
		s.np = len(s.P)
		
		s.setFlagNodes()
		s.setFlagsAndConnectivity()	
		
		s.Lx = np.max(s.P[:,0])
		s.Ly = np.max(s.P[:,1])
		s.Lz = np.max(s.P[:,2])		
		
	def createNetwork2D(s,param,op): # op = 0 is for asymFac and 1 is by thetaEstimate

		if(op == 0):
			nfibprevision = param[0]
			asymFac = param[1]
			s.nx,s.ny = estimateGridPointsByAsymetry(nfibprevision,asymFac)
		elif(op == 1):
			nfibprevision = param[0]
			meanTheta = param[1]
			s.nx,s.ny = estimateGridPointsByTheta(nFibPrevision,meanTheta)
			
		# update estimatives

		# actual values
		s.asymFac = s.nx - s.ny 
		s.meanTheta = np.arctan(float(s.nx-1)/float(s.ny-1))
			
		points = structuredPointsWithBoundary(s.nx,s.ny)
	 
		delau = Delaunay(points)
		
		s.P = delau.points # get points	
		s.delauSimplices = delau.simplices
		s.ElemFib = simplices2Edges(s.delauSimplices) # transform triangles into edges
		
		s.Lx = np.max(s.P[:,0])
		s.Ly = np.max(s.P[:,1])
		

	def createNetwork3D(s,param):
		
		s.nx = param[0]
		s.ny = param[1]
		s.nz = param[2]
		
		s.Lx = param[3]
		s.Ly = param[4]
		s.Lz = param[5]
		
		px = np.linspace(0.0,s.Lx,s.nx) 
		py = np.linspace(0.0,s.Ly,s.ny) 
		pz = np.linspace(0.0,s.Lz,s.nz) 
		
		Npoints = s.nx*s.ny*s.nz
		
		points = np.zeros((Npoints,s.ndim))
		
		count = 0
		for i in range(s.nx):
			for j in range(s.ny):
				for k in range(s.nz):
					points[count,0] = px[i]
					points[count,1] = py[j]
					points[count,2] = pz[k]
					count = count + 1			
		
		delau = Delaunay(points)
		
		s.P = delau.points
		s.delauSimplices = delau.simplices
		
		s.ElemFib = simplices2Edges(s.delauSimplices)

		

	def cleanVoronoiPoints(s,maxDist):
		
		newInd = np.zeros(s.np,dtype = 'int')

		k = 0
		for i in range(s.np):
			if(isOutUnitSquare(s.P[i,:],maxDist)): # cut exterior nodes not close enough 
				print 'cut node'
				newInd[i] = -1
				k = k + 1
			else:
				newInd[i] = i - k
		
		for f in s.ElemFib:
			f[0] = newInd[f[0]]
			f[1] = newInd[f[1]]
			
		s.ElemFib = cleanEdgesToInfinity(s.ElemFib)
			
		s.P = s.P[newInd>-1,:]
		s.np = len(s.P)

	def eliminateFibresOverBoundary(s):
		s.setFlagNodes()
		
		k = 0
		
		while k < len(s.ElemFib):
			#~ e = [0,1]
			e = s.ElemFib[k]
			
			if(s.flagNode[e[0]]>0 and s.flagNode[e[1]]>0):
				s.ElemFib.pop(k)
				k = k - 1
				
			k = k + 1
	


	def removeFibres(s,crit,param=[]):
				
		count = 0 
		k = 0 	
		print s.nf, s.getNumberDangledPoints()
		while(k < s.nf):	
			if(crit(s,k,param)):
				del s.ElemFib[k]
				s.nf = s.nf - 1
			else:
				k = k + 1

		print s.nf, s.getNumberDangledPoints()
		s.setFlagNodes()
		s.setFlagsAndConnectivity()
		print s.np, s.getNumberDangledPoints()
		s.cleanDangledPoints()
		print s.np, s.getNumberDangledPoints()
		s.setFlagNodes()
		s.setFlagsAndConnectivity()
		print s.getNumberDangledPoints()
		
	def modifyFibreProperty(s,label,foo,param,smoothFactor,binaryLevel = False):
		
		prop = s.chooseProperty(label) # points to the property vector
		
		if(binaryLevel): # first parameter is the contrast
			for i in range(s.nf): 
				if(foo(s,i,param)[1]):
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


		
	def getNumberDangledPoints(s):
		k = 0
		for i in range(s.np):
			if(len(s.InvConnect[i]) == 0): 
				k = k + 1
		
		return k
		
	# update the points by adding a random pertubation. For periodic the the average. 
	def addPertubation(s,pertPoint,isPeriodic=False,facx = 1.0,facy=1.0,facz = 1.0): 
		
		r = pertPoint*np.random.rand(s.np)
		theta = 2.0*np.pi*np.random.rand(s.np)
		if(s.ndim == 3):
			phi = np.pi*np.random.rand(s.np)
		elif(s.ndim == 2):
			phi = 0.5*np.pi*np.ones(s.np)
		
		dp = np.zeros((s.np,s.ndim)) 
		
		dp[:,0] = facx*r*np.sin(phi)*np.cos(theta)
		dp[:,1] = facy*r*np.sin(phi)*np.sin(theta)
		if(s.ndim == 3):
			dp[:,2] = facz*r*np.cos(phi)
		
		#~ s.P[s.flagNode<1,:] = s.P[s.flagNode<1,:] + dp[s.flagNode<1,:]
		
		s.P[:,0:s.ndim] = s.P[:,0:s.ndim] + dp[:,0:s.ndim]
		
		if(isPeriodic and s.ndim == 2): # doesn't to 3d yet
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
			
			
	def correctPoints(s,maxP,minP): # bugged function
		
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
			
		
	def writeNetwork(s,Ndof,Nsubsteps,Nparam,fignum,opInifile = 'SGP', opAdditionalElements = 0, addAuxNodes = 2, addAuxNodeOnFibres = 1):
		
		Paux = copy.deepcopy(s.P)
		if(addAuxNodes > 0):
			Padd = np.array(addAuxNodes*[s.ndim*[999.9]])
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
			writeMesh(Paux,[ElemFibAux],[s.delauSimplices,[[len(Paux)-1]]],Ndof,Nsubsteps)
		elif(opAdditionalElements == 2): # torsinal springs (angles) and global node
			writeMesh(Paux,[ElemFibAux],[s.ElemAngle,[[len(Paux)-1]]],Ndof,Nsubsteps)
		elif(opAdditionalElements == 3): # triangles, torsinal springs (angles) and global node
			writeMesh(Paux,[ElemFibAux],[s.delauSimplices,s.ElemAngle,[[len(Paux)-1]]],Ndof,Nsubsteps)
			
		os.system("cp Mesh.txt  Mesh" + str(fignum) + ".txt ")
		os.system("cp IniFile000.txt  IniFile000_" + str(fignum) + ".txt ")
		os.system("cp Param000.txt  Param000_" + str(fignum) + ".txt ")
		os.system("cp Bten_yG.txt  Bten_yG_" + str(fignum) + ".txt ")


	def get_normal_bar(s):
		
		vacc = np.zeros(s.ndim)
		Acc = 0.0
		
		for i in range(s.np):
			if(s.flagNode[i]>0):
				vacc = vacc + s.Abar[i]*s.normal[i,:]
				Acc = Acc + s.Abar[i]
		
		normal_bar = vacc/Acc
		
		return normal_bar

	def getCentreFibre(s,k):
		f = s.ElemFib[k]	
		x = 0.5*(s.P[f[0],:] + s.P[f[1],:])
		
		return x
	
	def get_yG(s):
		
		yG = np.zeros(s.ndim)
		
		for i in range(s.nf):
			yG = yG + s.Vf[i]*s.getCentreFibre(i)
						
		yG = yG/np.sum(s.Vf)
		
		return yG

	def get_Bten(s):
		
		Bten = np.zeros((s.ndim,s.ndim))
		
		for f in range(s.nf):
			Bten = Bten + s.Vf[f]*np.outer(s.af[f,:],s.af[f,:])
			
		Bten = Bten/np.sum(s.Vf)
		
		return Bten
	
	def cleanDangledPoints(s):
		
		s.np = len(s.P)
		newInd = np.zeros(s.np,dtype = 'int')
		
		print len(s.InvConnect)
		
		k = 0
		for i in range(s.np):
			if(len(s.InvConnect[i]) == 0): 
				newInd[i] = -1
				k = k + 1
			else:
				newInd[i] = i - k
		
		for e in s.ElemFib:
			e[0] = newInd[e[0]]
			e[1] = newInd[e[1]]
		
		s.P = s.P[newInd>-1,:]
		
		s.np = len(s.P)

