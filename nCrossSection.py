import ENDF6
from findpath import pathfinder
import numpy as np
from math import isnan, isinf
import os, sys
from bisect import bisect_left as bisect
from matplotlib import pyplot as plt

amu=1.660539040e-27 # atomic mass unit, used to find number densities

class isotope:
	def __init__(self,Z,A):
		self.Z=Z
		self.A=A
		self.m=A*amu
		self.read()
			
	def read(self):
		self.extrapE=[]
		self.extrapXS=[]
		rootpath=pathfinder()
		data=os.path.join(rootpath,'Data')
		f=open(os.path.join(data,'{0:0=3d}_{1:0=3d}.endf'.format(self.Z,self.A)))
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
		#finding titles for angle data
		self.angledata=[]
		for i in list(dataav):
			for j in self.XStype:
				if i[2]==j and i[1]==4:
					self.angledata.append(j)
		print('Finished processing isotope ({},{})'.format(self.Z,self.A))
		if len(self.extrapE)>0:
			if np.amin(self.extrapE)<14.1e6:
				print('Warning: extrapolated XS data at E<14.1MeV for isotope Z={}, A={}'.format(Z,A))
		self.notloaded=False
	
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
		XStot=np.interp(E,self.E[1],self.XSdata[1])
		if MT==1:
			return XStot
		elif XStot==0:
			return 0.0
		else:
			if MT in self.XStype:
				if E in self.E[MT]:
					ind=self.E[MT].index(E)
					XS=self.XSdata[MT][ind]
					return XS
				elif E<min(self.E[MT]) or E>max(self.E[MT]):
					return 0.0
				else:
					i=bisect(self.E[MT],E)
					XSint_n=self.XSdata[MT][i]
					XSint_p=self.XSdata[MT][i-1]
					if XSint_p==0 and XSint_n==0:
						return 0.0
					if self.E[MT][i] in self.E[1]:
						totind_n=self.E[1].index(self.E[MT][i])
						XStot_n=self.XSdata[1][totind_n]
					else:
						XStot_n=np.interp(self.E[MT][i],self.E[1],self.XSdata[1])
					if XStot_n==0:
						r_n=0.0
					else:
						r_n=XSint_n/XStot_n
					if self.E[MT][i-1] in self.E[1]:
						totind_p=self.E[1].index(self.E[MT][i-1])
						XStot_p=self.XSdata[1][totind_p]
					else:
						XStot_p=np.interp(self.E[MT][i-1],self.E[1],self.XSdata[1])
					if XStot_p==0:
						r_p=0.0
					else:
						r_p=XSint_p/XStot_p
					E_c=E-self.E[MT][i-1]
					E_l=self.E[MT][i]-self.E[MT][i-1]
					r=(r_p+(r_n-r_p)*E_c/E_l)
					XS=r*XStot
					try:
						return XS[0]
					except:
						return XS
			else:
				return 0.0
	
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
				self.iso[(Z,i)]=isotope(Z,i)
				self.comp[(Z,i)]=x
				self.avm+=self.iso[(Z,i)].m*self.comp[(Z,i)]
				for MT in self.iso[(Z,i)].XStype:
					if MT not in self.dataav:
						self.dataav.append(MT)
		else:
			self.iso[(Z,Aarray)]=isotope(Z,Aarray)
			self.comp[(Z,Aarray)]=1
			self.avm=self.iso[(Z,Aarray)].m
			self.dataav.extend(self.iso[(Z,Aarray)].XStype)
		self.dataav.sort()

class compound:
	def __init__(self,elemarray,elemcomp,rho,mix=None):
		importconds=all([elemarray==None,elemcomp==None,rho==None,mix is not None])
		if not importconds:
			try:
				self.defaultinit(elemarray,elemcomp,rho)
			except:
				raise AttributeError('Need input of elements, composition array and density of compound, or a mixture!')
		elif importconds:
			self.importfrommix(mix)
		else:
			print('Something failed')
	
	def defaultinit(self,elemarray,elemcomp,rho):
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
		self.rho=rho
		self.species={}
		self.dataav=[]
		if type(elemarray) is list:
			for e,c in zip(elemarray,elemcomp):
				for i in e.iso:
					self.species[i]=e.iso[i]
				for MT in list(e.dataav):
					if MT not in self.dataav:
						self.dataav.append(MT)
		else:
			for i in elemarray.iso:
				self.species[i]=elemarray.iso[i]
			self.dataav.extend(elemarray.dataav)
		self.dataav.sort()
		#finding the isotope composition in a compound
		self.isocomp={}
		if type(elemarray) is list:
			for elem,comp in zip(elemarray,elemcomp):
				for iso in elem.iso:
					if iso in self.isocomp:
						self.isocomp[iso]+=comp*elem.comp[iso]
					else:
						self.isocomp[iso]=comp*elem.comp[iso]
		else:
			for iso in elemarray.iso:
				if iso in self.isocomp:
					self.isocomp[iso]+=elemarray.comp[iso]
				else:
					self.isocomp[iso]=elemarray.comp[iso]
		totalcomp=sum(self.isocomp.values())
		self.isofrac={}
		for iso in self.isocomp:
			self.isofrac[iso]=self.isocomp[iso]/totalcomp
		fracsum=0
		for (Z,A) in self.isocomp:
			fracsum+=A*self.isofrac[(Z,A)]
		self.N=self.rho/(amu*fracsum)
		for iso in self.isocomp:
			self.isocomp[iso]=self.isofrac[iso]*self.N

	def totalXS(self,E,MT=1):
		if type(self.x) is list:
			totXS=0
			for x,iso in zip(self.x,self.species):
				try:
					totXS+=self.N*x*iso.XSfind(E,MT=MT)/10**28
				except:
					totXS+=0
		else:
			totXS=self.N*self.x*self.species.XSfind(E,MT=MT)/10**28
		return totXS
	
	def importfrommix(self,mix):
		#import number density
		self.N=mix.N
		#import data available
		self.dataav=mix.dataav
		'''
		if type(mix.comp_array)==list:
			for frac,comp in zip(mix.vol_frac_array,mix.comp_array):
				for x,iso in zip(comp.x,comp.species):
					if iso not in self.species:
						self.species.append(iso)
						self.x.append(frac*x)
					else:
						i=self.species.index(iso)
						self.x[i]+=frac*x
		else:
			for x,iso in zip(mix.comp_array.x,mix.comp_array):
				if iso not in self.species:
					self.species.append(iso)
					self.x.append(x)
				else:
					i=self.species.index(iso)
					self.x[i]+=x
		'''
		self.species=mix.species
		self.rho=mix.rho
		self.isocomp=mix.isocomp
		self.isofrac=mix.isofrac

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
		self.N=0
		self.comp_array=comp_array
		try:
			tot=sum(vol_frac_array)
			self.vol_frac_array=[i/tot for i in vol_frac_array] #normalising the volume fraction array
		except:
			self.vol_frac_array=1
		self.dataav=[]
		self.rho=0
		self.isocomp={}
		self.species={}
		if type(comp_array) is list:
			for comp,frac in zip(comp_array,self.vol_frac_array):
				self.N+=comp.N*frac
				for MT in comp.dataav:
					if MT not in self.dataav:
						self.dataav.append(MT)
				self.rho+=comp.rho*frac
				for iso in comp.isocomp:
					if iso in self.isocomp:
						self.isocomp[iso]+=frac*comp.isocomp[iso]
					else:
						self.isocomp[iso]=frac*comp.isocomp[iso]
						self.species[iso]=comp.species[iso]
		else:
			self.dataav.extend(comp_array.dataav)
			self.N=comp_array.N
			self.rho=comp_array.rho
			for iso in comp_array.isocomp:
				if iso in self.isocomp:
					self.isocomp[iso]+=comp_array.isocomp[iso]
				else:
					self.isocomp[iso]=comp_array.isocomp[iso]
					self.species[iso]=comp_array.species[iso]
		self.dataav+=[101]
		self.dataav.sort()
		totalN=sum(self.isocomp.values())
		self.isofrac={}
		self.isomass={}
		for iso in self.isocomp:
			self.isofrac[iso]=self.isocomp[iso]/totalN
			self.isomass[iso]=self.isocomp[iso]*amu*iso[1]
		if preload:
			self.XSgen(Eres)
		else:
			self.gen=False
		
	def mixXS(self,E,MT=1):
		avXS=0
		if MT!=101:
			for iso in self.species:
				avXS+=self.isocomp[iso]*self.species[iso].XSfind(E,MT=MT)/10**28
		else:
			for iso in self.species:
				for MT in range(102,118):
					avXS+=self.isocomp[iso]*self.species[iso].XSfind(E,MT=MT)/10**28
		return avXS
	
	def XSgen(self,Eres):
		if not self.gen:
			lowlim=np.log10(0.01)
			uplim=np.log10(14.1e6)
			exparr=np.linspace(lowlim,uplim,int(Eres))
			self.Earr=10**exparr
			self.XS={}
			for MT in self.dataav:
				self.XS[MT]=[]
				for E in self.Earr:
					self.XS[MT].append(self.mixXS(E,MT=MT))
				print('Finished generating reaction {}'.format(MT))
		self.gen=True
				
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
	