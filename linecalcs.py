import numpy as np
from math import isinf, inf

class regions(object):
	def __init__(self,point_arr_tple,mat,debug=False):
		self.debug=debug
		self.shape={'bound':{'r':point_arr_tple[0],'z':point_arr_tple[1],'m':mat}}
		self.shape['internal']=[]
		
	def innerbound(self,point_arr_tple,mat):
		self.shape['internal'].append({'r':point_arr_tple[0],'z':point_arr_tple[1],'m':mat})
		
	def interrogateinternal(self,start_tuple,dir_tuple):
		xo,yo,zo=start_tuple
		u,v,w=dir_tuple
		startarr=np.array([xo,yo,zo])
		vec=np.array([u,v,w])
		t=[]
		mats=[]
		L=[0]
		for i in range(len(self.shape['bound']['r'])):
			ri,zi=self.shape['bound']['r'][i-1],self.shape['bound']['z'][i-1]
			re,ze=self.shape['bound']['r'][i],self.shape['bound']['z'][i]
			m,c=decomposeline(ri,zi,re,ze)
			tt,nsol=intsctn(xo,yo,zo,u,v,w,c,m,debug=self.debug)
			for ti in tt:
				point=startarr+vec*ti
				#don't append intersection point if beyond the ranges of the line segment or if ti is negative
				if not(any([(point[0]**2+point[1]**2)**0.5>max(ri,re),point[2]>max(zi,ze),
					(point[0]**2+point[1]**2)**0.5<min(ri,re),point[2]<min(zi,ze),ti<0])):
					t.append(ti)
					mats.append(self.shape['bound']['m'])
		if len(t)>0:
			for shape in self.shape['internal']:
				for i in range(len(shape['r'])):
					ri,zi=shape['r'][i-1],shape['z'][i-1]
					re,ze=shape['r'][i],shape['z'][i]
					m,c=decomposeline(ri,zi,re,ze)
					tt,nsol=intsctn(xo,yo,zo,u,v,w,c,m,debug=self.debug)
					for ti in tt:
						point=startarr+vec*ti
						#don't append intersection point if beyond the ranges of the line segment or if ti is negative
						if not(any([(point[0]**2+point[1]**2)**0.5>max(ri,re),point[2]>max(zi,ze),
							(point[0]**2+point[1]**2)**0.5<min(ri,re),point[2]<min(zi,ze),ti<0])):
							t.append(ti)
							mats.append(shape['m'])
			sortedlist=sorted(zip(t,mats))
			t=[ti for ti,mi in sortedlist]
			mats=[mi for ti,mi in sortedlist]
			print(t,mats) if self.debug else ...
			for ti in t:
				print(u,v,w,ti) if self.debug else ...
				L.append((u**2+v**2+w**2)**0.5*ti-L[-1])
			return L,mats,t
	
class neutronstage(object):
	def __init__(self,start,vec,E,Rarr,l_rand):
		self.alive=True
		while self.alive:
			t,mat=self.findloc(start,vec,E,Rarr,l_rand)
			self.findreact(start,vec,E,mat,r_rand)
		
	def findloc(self,start,vec,E,Rarr,l_rand):
		L=[]
		mats=[]
		t=[]
		dens=1
		S=0
		for region in Rarr:
			Lt,mt,tt=region.interrogateinternal(start,vec)
			L.append(Lt)
			mats.append(mt)
			t.append(tt)
		t,mats,L=zip(sorted(zip(t,mats,L)))
		for tt,mt,Lt in zip(t,mats,L):
			for Li,mi in zip(Lt,mt):
				drop=np.exp(-Li*mt.mixXS(E))
				if dens*drop<l_rand:
					dS=mt.mixXS*np.log(dens/l_rand)
					Send=S+dS
					m_end=mi
				else:
					S+=Li
		try:
			tend=Send/(vec.dot(vec))**0.5
		except:
			tend=inf
			self.alive=False
		return tend,m_end
	
def intsctn(xo,yo,zo,u,v,w,ci,m,debug=False):
	if isinf(m): # if m is infinite, the line is purely horizontal
		#if so, there is one intersection
		if w!=0:
			t=[(ci/w)-zo]
		else:
			t=[ci]
		nsol=1
		return t,nsol
	#if not, the solution must be found via the quadratic solutions
	a=u**2+v**2-m**2*w**2
	b=2*(u*xo+v*yo-m*w*(ci+m*zo))
	c=xo**2+yo**2-ci**2-2*ci*m*zo-m**2*zo**2
	#If the lines are parallel, a=0
	if a==0:
		zeta0=-ci/m
		t0=-xo/u
		z0=zo+t0*w
		if z0==zeta0:
			#determine if the lines are coincident; if so there are infinite solutions
			nsol=3
			t=[inf,inf,inf]
		else:
			#if not, there are no solutions
			nsol=0
			t=[]
	else:
		if b**2<4*a*c:
			#other options are that the lines simply do not intersect and are not parallel
			t=[]
			nsol=0
		else:
			#Lastly, the most useful calculation; finding the intersection points
			det=(b**2-4*a*c)**0.5/(2*a)
			mid=-b/(2*a)
			if det==0:
				#There may be only one solution
				t=[mid]
				nsol=1
			else:
				#Or up to two unique solutions
				t=[mid+det,mid-det]
				nsol=2
	return t,nsol

def decomposeline(r1,z1,r2,z2):
	if z2!=z1:
		m=(r2-r1)/(z2-z1)
		c=r1-m*z1
	else:
		m=float('inf')
		c=z1 #if the line is horizontal, gradient is infinite and doesn't cross the r axis
		#c value is given as the z value only
	return m,c
