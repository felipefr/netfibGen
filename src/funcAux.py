# ================= auxiliary functions ==========================

#~ from generator import *
#~ from delaunay import *
from scipy.spatial import Delaunay
import numpy as np
import matplotlib.pyplot as plt
import copy
import os, sys
from IOfibres import *
import itertools

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
	
	
def uniqueEdges(e):
	e = map(tuple, e)
	
	# eliminates the order, assumes the default order of sets  
	e = map( lambda x: x if(x[0]<x[1]) else (x[1],x[0]) , e) 
	 # eliminates duplicates by taking the set, later convert set to list, and tuple to list (for each edge)
	e = map(list,list(set(e)))
	
	return e
	
def unique_rows(a):
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))	



def unique_rows2(a):
	tol = 9.0e-2
	na = len(a)
		
	flag = range(na)

	for i in range(na):
		ai = a[i,:]
				
		dij = np.linalg.norm(a[i+1: ,:] - ai,axis = 1)
		
		for j in range(i+1,na):
			if(dij[j-i-1]<tol and flag[j]>flag[j-1]):
				flag[j] = i

			
			#~ aj = a[j]
			#~ dij = np.linalg.norm(ai-aj)
			
	flag = list(set(flag))
	
	return a[flag]   

def simplices2Edges(simpList):
	
	e = []
	
	for t in simpList:		
		e_t = list(itertools.combinations(t,2))		
		e = e + e_t

	e = uniqueEdges(e) 
	
	return e
	
def ciclicCombination(l):
	n = len(l)
	if(n<2):
		return []
	
	e_l = [(l[0],l[n-1])]
	
	for i in range(n-1):
		e_l = el + [(l[i],l[i+1])] 
	
	return e_l
	
def voronoiRegions2Edges(regions):
	
	e = []
	
	for t in regions:		
		e_t = ciclicCombination(t)		
		e = e + e_t

	e = uniqueEdges(e) 
	
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
def classify(x,L=1.0):
	
	tol = 1.0e-5
	c = 0
	
	if(abs(x)<tol): # left extreme
		c = 0
	elif(abs(x-L)<tol): # right extreme
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
	lf = np.linalg.norm(dx)
	
	return dx/lf , lf

def getNormal_kmax2D(flag,ndim):
	
	# for points : -1 = neighbour to boundary, 0 =  very interior, 1 = bottom, 2 = right, 3 = top, 4 = left
	# 5 = bottom/left , 6 = bottom/right, 7 = right/top, 8 = top/left
	kmax = ndim
	
	kmax_i = 0
	normal_ik = np.zeros((kmax,ndim)) 

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
	
def flag2Ind(flag):
	aux = flag%9
	iz = (flag-aux)/9
	ix = aux%3
	iy = (aux - ix)/3 
	
	iz = 2 - iz
	iy = 2 - iy
	ix = 2 - ix
	
	ind = [ix,iy,iz]

	return ind 
	
def getNormal_kmax3D(flag,ndim):
	
	normal_ik = np.eye(ndim) 
	
	ind = flag2Ind(flag)
	indNormal = []
	
	# the number of effective normal vector is the number of indexes different than 2
	for p in range(ndim):
		if(ind[p] == 2):
			normal_ik[p,p] = 0.0
		else:
			if(ind[p] == 0):
				normal_ik[p,p] = -normal_ik[p,p]
			
			indNormal = indNormal + [p]
	
	kmax_i = len(indNormal)
	normal_ik = normal_ik[indNormal,:]

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
	
def convert2LineEq(x0,x1): # helps with the implementation of stripDecrease.
	x0 = np.array(x0)
	x1 = np.array(x1)

	p = x1 - x0

	normal = np.zeros(2)
	
	normal[0] = p[1]
	normal[1] = -p[0]
	
	normal = normal/np.linalg.norm(normal)
	
	a = normal[0]
	b = normal[1]
	
	c = np.dot(normal,x0)
	
	return [a,b,c]
	
def convert2PlaneEq(x0,x1,x2): # helps with the implementation of stripDecrease.
	x0 = np.array(x0)
	x1 = np.array(x1)
	x2 = np.array(x2)

	p1 = x1 - x0
	p2 = x2 - x0

	normal = np.cross(p1,p2)

	normal = normal/np.linalg.norm(normal)
	
	a = normal[0]
	b = normal[1]
	c = normal[2]
	
	d = np.dot(normal,x0)
	
	return [a,b,c,d]
	
def planeDecrease(net,k,param): # this should be generalised for 3D
	
	a = param[0]
	b = param[1]
	c = param[2]
	d = param[3]
	delta = param[4]
	
	x0,y0,z0 = net.getCentreFibre(k)
	
	dist = np.abs( (d - a*x0 - b*y0 - c*z0)/np.sqrt(a*a + b*b + c*c) )
	
	gamma = 0.5 # later it will be amplified
	
	if(dist<delta):
		smoothFac = gamma + (1.0 - gamma)*(dist/delta)**2.0
		return smoothFac, True
		 
	return 1.0 , False 

def stripDecrease(net,k,param): # this should be generalised for 3D
	
	a = param[0]
	b = param[1]
	c = param[2]
	delta = param[3]
	
	x0,y0 = net.getCentreFibre(k)
	
	d = np.abs( (c - a*x0 - b*y0)/np.sqrt(a*a + b*b) )
	
	gamma = 0.5 # later it will be amplified

	if(d<delta):
		smoothFac = gamma + (1.0 - gamma)*(d/delta)**2.0
		return smoothFac, True
		 
	return 1.0 , False 
	

def stripDecreasePoints(net,k,param): # p(t) = p1 + t*(p2 - p1)
	
	p1 = np.array(param[0])
	p2 = np.array(param[1])
	delta = param[2]
	
	dp21 = p2 - p1
	p0 = net.getCentreFibre(k) # point to test
	
	t = np.dot(p0-p1,dp21)/np.dot(dp21,dp21) # t that minimises || p(t) - p0|| 
	pt = p1 + t*dp21
	
	d = np.linalg.norm(pt - p0)
		
	gamma = 0.5 # later it will be amplified

	if(d<delta and t>0.0 and t<1.0):
		smoothFac = gamma + (1.0 - gamma)*(d/delta)**2.0
		return smoothFac, True
		 
	return 1.0 , False 


def ballDecrease(net,k,param):
	x0 = np.zeros(net.ndim)
	x0[:] = param[0:net.ndim]
	R = param[net.ndim]
	
	p = net.getCentreFibre(k)
	d = np.linalg.norm(p-x0) 
		
	smoothFac = 1.0
	gamma = 0.5 # the default gamma is 0.5 (outside is gonna be adjusted)

	if(d < R) :
		smoothFac = gamma + (1.0 - gamma)*(d/R)**2.0
		return smoothFac, True
		
	return 1.0, False
		

def InBallFibresCriterium(net,k,param):	
	x0 = np.zeros(net.ndim)
	x0[:] = param[0:net.ndim]
	R = param[net.ndim]
	
	p = net.getCentreFibre(k)
	d = np.linalg.norm(p-x0)
	
	if(d<R):
		return True

	return False

def VertHorizFibresCriterium(net,k,param=[]):
	if(len(param)>0):
		tol = param[0]
	else:
		tol = 1.0e-3
		
	f = net.ElemFib[k]
	dx = net.P[f[0],:] - net.P[f[1],:]
			
	if(np.abs(np.prod(dx)) < tol):
		return True

	return False

def DangledFibresCriterium(net,k,param=[]):	
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
	

def generatePointsMinDist(n,minDist,NdimE = 2):
	
	P = np.zeros((n,NdimE))

	P[0,:] =  np.random.rand(NdimE)
	k = 1
	while (k<n):	
	
		Ptest = np.random.rand(NdimE)
	
		dp = P[:k,:] - Ptest[:]
		
		dist = np.linalg.norm(dp,axis = 1)
		
		if(np.min(dist) > minDist):
			P[k,:] = Ptest
			k = k+1
	
	return P

def cleanEdgesToInfinity(Edges): # represented as -1
	EdgesNew = []
	for ei in Edges:
		isValid = True
		for v in ei:
			if(v == -1):
				isValid = False
				break
			
		if(isValid):
			EdgesNew = EdgesNew + [ei]
	
	return EdgesNew
	
def isOutUnitSquare(pt,delta=0.0):
	r = False
	
	for pti in pt:
		if(pti<delta or pti>1.0-delta):
			r = True

	return r 
	
def borderPointsVoronoi2D(a,b,c,d,nx,ny):
	
	x = np.linspace(a,c,nx)
	y = np.linspace(b,d,ny)
	
	N = 2*(nx + ny - 2)
	
	p = np.zeros((N,2))
	
	k = 0
	for i in range(nx):
		for j in range(ny):
			
			if(i*j == 0 or (nx-1-i)*(ny-1-j)==0):				
				p[k,0] = x[i]
				p[k,1] = y[j] 
				k = k + 1
		
	return p
	
def borderPointsVoronoi3D(a,b,c,d,e,f,nx,ny,nz):
	
	x = np.linspace(a,d,nx)
	y = np.linspace(b,e,ny)
	z = np.linspace(c,f,nz)
	N = 2*(nx*ny + ny*nz + nx*nz - 2*(nx+ny+nz) + 4)
	
	p = np.zeros((N,3))
	
	count = 0
	for i in range(nx):
		for j in range(ny):
			for k in range(nz):
				
				if(i==0 or i==nx-1 or j==0 or j==ny-1 or k==0 or k==nz-1):				
					p[count,0] = x[i]
					p[count,1] = y[j] 
					p[count,2] = z[k] 
					count = count + 1
	
	
	return p

def connectInteriorWithBoundaryNodes(ptBound,ptInt):
	
	newFibres = []
		
	npi = len(ptInt)
	npb = len(ptBound)

	connectivity = 3
	minDist = 0.5

	for i in range(npb):
		dij = np.linalg.norm(ptBound[i,:]-ptInt,axis = 1)
		arg = np.argsort(dij)
		minDist = minDist*np.mean(dij)
		
		listNeighbours = [arg[0]]
		k = 0
		while(len(listNeighbours)<connectivity):
			k = k + 1
			pk = ptInt[arg[k],:]
			
			isValid = True
			for j in listNeighbours:
				pj = ptInt[j]
				if(np.linalg.norm(pk-pj)<minDist):
					isValid = False
					break
			
			if(isValid):
				listNeighbours = listNeighbours + [arg[k]]
			
			
		for j in listNeighbours:
			newFibres = newFibres + [[i + npi,j]]
			
	return newFibres

#~ def connectInteriorWithBoundaryNodes(ptBound,ptInt):
	
	#~ newFibres = []
		
	#~ npi = len(ptInt)
	#~ npb = len(ptBound)

	#~ connectivity = 3

	#~ for i in range(npb):
		#~ dij = np.linalg.norm(ptBound[i,:]-ptInt,axis = 1)
		#~ arg = np.argsort(dij)
		#~ for j in range(connectivity):
			#~ newFibres = newFibres + [[i + npi,arg[j]]]
			
	#~ return newFibres


def convertVoroPP2fibres(voro,points=[]):
	npts = len(points)
	
	# creating global vector
	Vert = []
	for cell in voro:
		Vert = Vert + cell['vertices']

	nvert = len(Vert)
	print nvert
	Vert = unique_rows2(np.array(Vert))
	
	nvert = len(Vert)

	print nvert
	#~ print raw_input() 

	# creating cell to global correspondence
	mapCell2Global = []
	for cell in voro:
		Pcell = cell['vertices']
		
		l = []
		for p in Pcell:
			dij = np.linalg.norm(Vert - p,axis = 1)
			
			l = l + [dij.argmin()]
		
		mapCell2Global = mapCell2Global + [l]
	
	#~ print mapCell2Global

	# creating Edges
	Edges = []
	for i,cell in enumerate(voro):
		adj = cell['adjacency'] 
		m = mapCell2Global[i]
		
		for j, aj in enumerate(adj):
			for k in aj:
				e = (m[j], m[k])
				if(get_afib_lf(Vert,e)[1]>1.0e-5):	
					Edges = Edges + [e]
				
		if(npts>0):
			for j in range(len(m)):
				e = (nvert + i, m[j]) # points will be concatenated with vertices
				if(get_afib_lf(Vert,e)[1]>1.0e-5):	
					Edges = Edges + [e]

	Edges = uniqueEdges(Edges)	
	
	if(npts==0):
		return Vert , Edges		
	else:
		return np.concatenate((Vert,points)) , Edges
	
	
def recombineEdges(Edges,nNew):
	
	nf = len(Edges)
	
	newEdges = []
	
	for i in range(nNew):
	
		e1 = np.random.randint(nf)
		e2 = np.random.randint(nf)
	
		if(e1 == e2):
			e2 = e1 + 1
			if(e2==nf):
				e2 = 0 
				
		e1 = Edges[e1]
		e2 = Edges[e2]
	
		newEdges = newEdges + [ [e1[0],e2[1]] , [e2[0],e1[1]] ] 
	
	
	return newEdges
	
