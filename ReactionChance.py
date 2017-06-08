from nova.neutronics import nCrossSection as ncs
import datetime
import random
import numpy as np
e=np.exp(1)

class reactions(object):
	def __init__(self,seed=None):
		if seed='today':
			dateseed=int(datetime.datetime.now().strftime('%d%M%y'))
			self.randgen=random.seed(dateseed)
		elif type(seed)==int or type(seed)==float:
			self.randgen=random.seed(seed)
			
	def findreaction(self,Larr,ncsarr,E):
		RN3=self.randgen.random()
		#popfrac=1.0 #calculation works on ruling out percentiles, which are 1-popfrac
		Lsum=0
		Sigsum=0
		'''
		if first section is likely to react with 40% of neutrons and
			the neutron is in the 80th percentile, more sections must be checked until 80% is reached
		'''
		for L,ncs in zip(Larr,ncsarr):
			Sigsum+=L*E
			XStot=ncs.mixXS(E)
			popfrac=1-e**-Sigsum
			#popfrac*=e**-(XStot*L)
			if RN3<1-popfrac:
				dSigsum=-np.log(1-RN3)-Sigsum
				Ltot=Lsum+dSigsum/XStot
				RN4=self.randgen.random()
				for r in ncs.XS:
					XSpart.append(ncs.mixXS(E,MT=r)/XStot)
					if sum(XSpart)>=RN4:
						self.reaction(E,MT=r)
						r_end=r
						break
				return Ltot,r_end
			elif RN3>1-popfrac:
				Sigsum+=L*XStot
				Lsum+=L
				continue
	
	
	
