import ENDF6
import numpy as np
from matplotlib import pyplot as plt

amu=1.660539040e-27 # atomic mass unit, used to find number densities

class isotope:
	def __init__(self,Z,A):
		self.Z=Z
		self.A=A
		self.extrapE=[]
		self.extrapXS=[]
		self.m=A*amu
		f=open('Data\\{0:0=3d}_{1:0=3d}.endf'.format(Z,A))
		lines=f.readlines()
		#finding cross sections available
		dataav=ENDF6.list_content(lines)
		self.XStype=[]
		self.XSdata={}
		self.E={}
		for i in list(dataav):
			if i[1]==3:
				self.XStype.append(i[2])
				sec=ENDF6.find_section(lines,MF=3,MT=i[2])
				self.E[i[2]],self.XSdata[i[2]]=ENDF6.read_table(sec)
		self.XStotal=self.XSdata[1] # total cross section always determined by MT=1
		self.XSfrac={}
		for arr in self.XSdata:
			self.XSfrac[arr]=[]
			for E in self.E[arr]:
				try:
					self.XSfrac[arr].append(self.XSfind(E,MT=arr)/self.XSfind(E))
				except:
					self.XSfrac[arr].append(0)
		#finding titles for angle data
		self.angledata=[]
		for i in list(dataav):
			for j in self.XStype:
				if i[2]==j and i[1]==4:
					self.angledata.append(j)
		if len(self.extrapE)>0:
			if np.amin(self.extrapE)<14.1e6:
				print('Warning: extrapolated XS data at E<14.1MeV for isotope Z={}, A={}'.format(Z,A))
	
	def XSplot(self,MT=1):
		try:
			if type(MT) is list:
				for i in MT:
					plt.loglog(self.E[i],self.XSdata[i],label='{},{}:{}'.format(self.Z,self.A,i))
			else:
				plt.loglog(self.E[MT],self.XSdata[MT],label='{},{}:{}'.format(self.Z,self.A,MT))
		except:
			...
		plt.xlabel('Energy ($eV$)')
		plt.ylabel('Microscopic cross section (barn or $b$)')
		plt.legend(loc='best')
	
	def XSfind(self,E,MT=1):
		if E<np.amin(self.E[MT]) or E>np.amax(self.E[MT]):
			if MT!=1:
				XS=0
			elif MT==1:
				XS=np.interp(E,self.E[1],self.XSdata[1])
				self.extrapXS.append(XS)
				self.extrapE.append(E)
		else:
			XS=np.interp(E,self.E[MT],self.XSdata[MT])
		return XS
	
class element:
	def __init__(self,Z,Aarray,Acomp):
		'''
		Always returns microscopic cross sections of mixtures
		Multiply by bulk atomic density to determine macroscopic cross section of mixtures
		
		Inputs:
		Z is atomic number (proton number)
		Aarray is an integer or list of integers of the atomic masses of isotopes in the elemental mixture
		Acomp is of the same type as Aarray and defines the composition of each isotope, in percent
		
		'''
		if type(Acomp)==list:
			Acomp=[i/sum(Acomp) for i in Acomp]
		self.dataav=[]
		self.avm=0
		self.iso={}
		self.comp={}
		if type(Aarray) is list:
			for i,x in zip(Aarray,Acomp):
				self.iso[i]=isotope(Z,i)
				self.comp[i]=x/100
				self.avm+=self.iso[i].m*self.comp[i]
				for MT in self.iso[i].XStype:
					if MT not in self.dataav:
						self.dataav.append(MT)
		else:
			self.iso[Aarray]=isotope(Z,Aarray)
			self.comp=[1]
			self.avm=self.iso[Aarray].m
			self.dataav.extend(self.iso[Aarray].XStype)
		self.dataav.sort()

class compound:
	def __init__(self,elemarray,elemcomp,rho):
		'''
		Inputs:
		elemarray is an integer or list of element objects
			These element objects must be initialised before generating a compound
		elemcomp is the composition fraction of each item in elemarray,
			which can be either a percentage of total (adding up to 100%)
			or as in a chemical formula (i.e. H2O would be 2,1, U3O8 would be 3,8)
			the code takes the elemarray and normalises it for appropriate ratios
		rho is density in kg/m3
		'''
		if type(elemcomp)==list:
			tot=sum(elemcomp)
			elemcomp=[i/tot for i in elemcomp]		
		self.x=[]
		self.m=[]
		self.species=[]
		self.dataav=[]
		if type(elemarray) is list:
			for e,c in zip(elemarray,elemcomp):
				if type(e.iso) is dict:
					for ind in e.iso:
						self.x.append(c*e.comp[ind])
						self.m.append(e.iso[ind].m)
						self.species.append(e.iso[ind])
				else:
					self.x.append(c)
					self.m.append(e.iso.m)
					self.species.append(e.iso)
				for MT in list(e.dataav):
					if MT not in self.dataav:
						self.dataav.append(MT)
		else:
			self.x.append(1)
			for ind in elemarray.iso:
				self.m.append(elemarray.iso[ind].m)
				self.species.append(elemarray.iso[ind])
			self.dataav.extend(elemarray.dataav)
		
		self.dataav.sort()
		xsum=0
		for i,j in zip(self.x,self.m):
			xsum+=i*j
		self.N=rho/xsum
		#finding the isotope composition in a compound
		self.isocomp={}
		if type(elemarray) is list:
			for elem,comp in zip(elemarray,elemcomp):
				for iso in elem.iso:
					if iso in self.isocomp:
						self.isocomp[iso]+=self.N*comp*elem.comp[iso]*100
					else:
						self.isocomp[iso]=self.N*comp*elem.comp[iso]*100
		else:
			for iso in elemarray.iso:
				if iso in self.isocomp:
					self.isocomp[iso]+=self.N*elemcomp*elemarray.comp[iso]*100
				else:
					self.isocomp[iso]=self.N*elemcomp*elemarray.comp[iso]*100

	def totalXS(self,E,MT=1):
		if type(self.x) is list:
			totXS=0
			for x,isp in zip(self.x,self.species):
				try:
					totXS+=self.N*x*isp.XSfind(E,MT=MT)/10**28
				except:
					totXS+=0
		else:
			totXS=self.N*self.x*self.species(E,MT=MT)/10**28
		return totXS
	
class mixture:
	def __init__(self,comp_array,vol_frac_array,Eres=10**3,preload=False):
		'''
		The vol_array is a number or list of numbers giving the volume fraction for each compound
		For a pebble bed, you may describe the compounds as:
			Pebbles (Be)
			Fill gas (e.g. He)
			Structural support (Eurofer)
		Each of these are compounds, made up of elements
		The mixture is then generated as such:
		pebblebed=nCrossSection.mixture([pebbles,gas,SS],[75,15,10])
		'''
		self.N=[]
		self.comp_array=comp_array
		tot=sum(vol_frac_array)
		self.vol_frac_array=[i/tot for i in vol_frac_array] #normalising the volume fraction array
		self.dataav=[]
		if type(comp_array) is list:
			for comp,frac in zip(comp_array,vol_frac_array):
				self.N.append(comp.N*frac)
				for MT in comp.dataav:
					if MT not in self.dataav:
						self.dataav.append(MT)
		else:
			self.dataav.extend(comp_array.dataav)
		self.dataav.sort()
		if preload:
			lowlim=np.log10(0.01)
			uplim=np.log10(14.1e6)
			exparr=np.linspace(lowlim,uplim,int(Eres))
			self.Earr=10**exparr
			self.XSgen()
		
		self.isocomp={}
		if type(comp_array)==list:
			for comp,v_frac in zip(comp_array,self.vol_frac_array):
				for iso in comp.isocomp:
					if iso in self.isocomp:
						self.isocomp[iso]+=v_frac*comp.isocomp[iso]
					else:
						self.isocomp[iso]=v_frac*comp.isocomp[iso]
		else:
			for iso in comp_array.isocomp:
				if iso in self.isocomp:
					self.isocomp[iso]+=v_frac*comp_array.isocomp[iso]
				else:
					self.isocomp[iso]=v_frac*comp_array.isocomp[iso]
	
	def mixXS(self,E,MT=1):
		if type(self.comp_array) is list:
			avXS=0
			for comp,frac in zip(self.comp_array,self.vol_frac_array):
				avXS+=comp.totalXS(E,MT=MT)*frac
		else:
			avXS=self.comp_array.totalXS(E,MT=MT)*self.vol_frac_array
		return avXS
	
	def XSgen(self):
		self.XS={}
		for MT in self.dataav:
			self.XS[MT]=[]
			for E in self.Earr:
				self.XS[MT].append(self.mixXS(E,MT=MT))
				
	def XSplot(self,MT=1,name=None):
		try:
			if type(MT)==list:
				for r in MT:
					if type(name)==str:
						plt.loglog(self.Earr,self.XS[r],label='{}:MT{}'.format(name,r))
					else:
						plt.loglog(self.Earr,self.XS[r],label='MT{}'.format(r))
			else:
				if type(name)==str:
					plt.loglog(self.Earr,self.XS[MT],label='{}:MT{}'.format(name,MT))
				else:
					plt.loglog(self.Earr,self.XS[MT],label='MT{}'.format(MT))
			plt.xlabel('Neutron energy ($eV$)')
			plt.ylabel('Macroscopic Cross Section ($1/m$)')
			plt.legend(loc='best')
		except:
			...