import numpy as np

def line(p1, p2):
    A = (p1[1] - p2[1])
    B = (p2[0] - p1[0])
    C = (p1[0]*p2[1] - p2[0]*p1[1])
    return A, B, -C

def intersection(L1, L2):
    D  = L1[0] * L2[1] - L1[1] * L2[0]
    Dx = L1[2] * L2[1] - L1[1] * L2[2]
    Dy = L1[0] * L2[2] - L1[2] * L2[0]
    if D != 0:
        x = Dx / D
        y = Dy / D
        return x,y
    else:
        return False

def fullpath(start,unitvector,bounds,linesize=0.01):
	xmax=bounds['x']
	zmax=bounds['z']['max']
	zmin=bounds['z']['min']
	pathx,pathy,pathz=[],[],[]
	pathx[0]=start[0]
	pathy[0]=start[1]
	pathz[0]=start[2]
	goodtogo=True
	while goodtogo:
		pathx.append(pathx[-1]+unitvector[0]*linesize)
		pathy.append(pathx[-1]+unitvector[1]*linesize)
		pathz.append(pathx[-1]+unitvector[2]*linesize)
		if (pathx[-1]**2+pathy[-1]**2)**0.5<xmax or zmin<=pathz[-1]<=zmax:
			goodtogo=False
	return pathx,pathy,pathz

def reducedims(x,y):
    if not(hasattr(x,'__iter__')) or not(hasattr(y,'__iter__')):
        raise AttributeError('need iterable')
    xeq=[(i**2+j**2)**0.5 for i,j in zip(x,y)]
    return xeq

def thetaarray(x2D,z):
    if not(hasattr(x,'__iter__')) or not(hasattr(y,'__iter__')):
        raise AttributeError('need iterable')
    theta=[]
    for i in range(len(x2D)-1):
        theta.append(np.arctan2(x2D[i+1]-x2D[i],z[i+1],z[i])])
    return theta
    
def hitpoint(pathst,pathe,lp1,lp2):
    '''
    Finds a hitpoint between a ray (at a start point) and a finite line in space.
    
    Input is the starting position of the ray, the 3D direction given by theta, where theta is the angle of the ray (at this point) to the vertical (z direction)
    '''
    ray=line(pathst,pathen)
    trace=line(lp1,lp2)
    intsx=intersection(ray,trace)
    
    if intsx is not False:
        if np.amin(pathst[0],pathen[0])<=ints[0]<=np.amax(pathst[0],pathen[0]) and \
        np.amin(pathst[1],pathen[1])<=ints[1]<=np.amax(pathst[1],pathen[1]) and \
        np.amin(lp1[0],lp2[0])<=ints[0]<=np.amax(lp1[0],lp2[0]) and\
        np.amin(lp1[1],lp2[1])<=ints[1]<=np.amax(lp1[1],lp2[1]):
            return intsx
        else:
            return False
    else:
        return False