#import nova.neutronics.nCrossSection as ncs
import numpy as np
from shapely import geometry as geom
import random

def perp( a ) :
    b = np.empty_like(a)
    b[0] = -a[1]
    b[1] = a[0]
    return b

# line segment a given by endpoints a1, a2
# line segment b given by endpoints b1, b2
# return 
def seg_intersect(a1,a2, b1,b2) :
    da = a2-a1
    db = b2-b1
    dp = a1-b1
    dap = perp(da)
    denom = np.dot( dap, db)
    num = np.dot( dap, dp )
    return (num / denom.astype(float))*db + b1

def intersect(A,B,checktype='all',returntype='normal'):
	fwdind=[]
	bckind=[]
	x=[]
	y=[]
	for i in range(len(Ax[:-1])):
		breakcheck=False
		u=[Ax[i+1]-Ax[i],Ay[i+1]-Ay[i]]
		for j in range(len(Bx[:-1])):
			v=[Bx[j+1]-Bx[j],By[i+1]-By[i]]
			d=(By[j]-Ay[i]+u[0]*(Ax[i]-Bx[j])/u[1])/(v[1]-u[1]*v[0]/u[0])
			c=(Ax[i]+v[0]*d-Bx[j])/u[0]
			if c>=0 and c<=1 and d>=0 and d<=1:
				fwdind.append(i+c)
				bckind.append(j+d)
				x.append(Ax[i]+u[0]*c)
				y.append(Ay[i]+u[1]*c)
				breakcheck=True
				break
		if breakcheck and checktype=='all':
			break
	if type(x) is None:
		return False
	else:
		if returntype=='normal':
			return [x,y],fwdind
		elif returntype=='coord':
			return x,y
		elif returntype=='int':
			return int(fwdind)

class npath(object):
	def __init__(self,bl,coarse=True,convergence=0.1,seed=None):
		self.xpath,self.ypath=bl.xpath,bl.ypath
		self.BSSx,self.BSSy=bl.supx.tolist(),bl.supy.tolist()
		self.setHCPBmats()
		if coarse:
			self.blshapegen()
		else:
			self.blfineshapegen()
		
		if type(convergence) is int:
			self.samples=convergence
		else:
			self.cc=convergence/100
		random.seed(a=seed)
		
	def setHCPBmats(self,enrich=70):
		PaFr=0.63
		H=ncs.element(1,1,1)
		Be=ncs.element(4,9,1)
		Li_en=ncs.element(3,[6,7],[90,10])
		Si=ncs.element(14,[28,29,30],[92.2,4.7,3.1])
		O=ncs.element(8,16,1)
		He=ncs.element(2,4,1)
		self.He_B=ncs.compound([He,H],[99.9,0.1],HeDens(1.1e5,673.15))
		self.He_M=self.He_B
		self.He_C=ncs.compound(He,1,HeDens(80e5,673.15))
		self.Be_sol=ncs.compound(Be,1,1850.0)
		self.Li4SiO4=ncs.compound([Li_en,Si,O],[4,1,4],2400.0)
		self.BePebbleBed=ncs.mixture([self.Be_sol,self.He_B],[PaFr,1-PaFr],preload=True)
		self.LiPebbleBed=ncs.mixture([self.Li4SiO4,self.He_M],[PaFr,1-PaFr],preload=True)
		#defining eurofer
		Fe=ncs.element(26,[54,56,57],[5.85,91.75,2.12])
		Cr=ncs.element(24,[50,52,53,54],[4.3,83.8,9.5,2.4])
		W=ncs.element(74,[180,182,183,184,186],[0.1,26.5,14.3,30.6,28.4])
		V=ncs.element(23,51,1)
		FeCont=100-(9+1.2+0.2+0.14)
		self.Eurofer=ncs.compound([Fe,Cr,W,V],[FeCont,9,1.2,0.2,0.14],7798.0)
		self.FW=ncs.mixture([self.Eurofer,self.He_C],[600,99],preload=True)
		#rough area of He channel (99mm^2) compared to metal (600mm^2)
		self.mani=ncs.mixture([self.Eurofer,self.He_C],[455,169],preload=True)
		#rough area of He channel in FW is 169mm^2, compared to 455mm^2 of SS
#		self.breeder=ncs.mixture([self.Be_sol,self.He_C,self.Eurofer],[231,12.5,22.5],preload=True)
#		self.multiplier=ncs.mixture([self.Li4SiO4,self.He_C,self.Eurofer],[77,12.5,22.5],preload=True)
		self.homogenous=ncs.mixture([self.Be_sol,self.He_C,self.Eurofer,self.Li4SiO4],\
			[231,25,45,77],preload=True)
		
	def blshapegen(self):
		point=[]
		self.mat=[]
		for j,(xarr,yarr) in enumerate(zip(self.xpath,self.ypath)):
			point.append(())
			for x,y in xarr,yarr:
				point[j]+=(x,y) # adding all points around the blanket box
			point[j]+=(xarr[0],yarr[0]) #closing the loop
			self.mat.append(self.homogenous)
		point.append(())
		for x,y in zip(self.BSSx,self.BSSy):
			point+=(x,y)
		self.mat.append(self.mani)
		self.blregion=geom.MultiLineString(point)
		
	def pathintersect(self,SP,EP,fidelity=100):
		point3D=[]
		point2D=[]
		for i in range(fidelity):
			point3D.append(SP+(EP-SP)*(i+1)/fidelity)
			req=np.sqrt(point3D[i][0]**2+point3D[i][2]**2)
			point2D.append(geom.Point([req,point3D[1]]))
		path2D=LineString(point2D)
		intersectlist=path2D.intersection(self.blregion)
		r=[]
		z=[]
		for p in list(intersectionlist):
			r.append(p.x)
			z.append(p.y)
		l=[]
		for i in range(len(z)-1):
			l.append(np.sqrt((r[i+1]-r[i])**2+(z[i+1]-z[i])**2))
		return l

	def path(self,SP):
		dirrand=random.uniform(0,2*np.pi) #random theta angle measured from vertically upwards
		u=np.array([np.sin(dirrand),np.cos(dirrand),0]) #u is velocity vector, directions are r, z and n, respectively
		EP=SP+100*u # calculating line end point, taken to be 100 metres away, suitably large
		alive=True
		while alive is True:
			Larr=self.pathintersect(SP,EP)
			interrand=random.random()
			
		