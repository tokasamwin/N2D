import linecalcs as n2, ReactionChance as reacha, nCrossSection as ncs

rarr=[1,11,11,1]
zarr=[0,0,10,10]

C=ncs.element(6,12,1)
graphite=ncs.compound(C,1,2000.0)
sol=ncs.mixture(graphite,1)

region=n2.regions((rarr,zarr),sol,debug=True)
r_in=[3,9,9,3]
z_in=[2,2,8,8]
diamond=ncs.compound(C,1,5000.0)
sol2=ncs.mixture(diamond,1)
region.innerbound((r_in,z_in),sol2)

xo,yo,zo=0,0,5
u,v,w=1,0,0
