#test to load melt data
import espressopp
import mpi4py.MPI as MPI
import sys
import os
import time
from espressopp import Real3D

num_chains = 1000
monomers_per_chain = 40
temperature1 = 1.0
temperature2 = 3.0
rc=pow(2.0,1.0/6.0)
skin=0.3
gamma = 1.5
dt = 0.001
bending = 0
tin = 0 #- initial time

if bending:
	kname = 'k1p5'
	ktheta = 1.5
else:
	kname = 'kzero'
	ktheta = 0.0

#inconfpath = '/data/pckr193/smrek/projects/phaseg/melt/data/configurations/N'+str(monomers_per_chain)+'/gamvar/'+kname+'/prc/'#
#inconffile = inconfpath+'prc_nc'+str(num_chains)+'N'+str(monomers_per_chain)+'dT'+str(temperature2)+'g'+str(gamma)+'dt0.005t'+str(tin)+'.dat'
inconffile = 'prc_nc1000N40dT3.0g1.5dt0.005t0.dat'
outconfpath = '/Users/niktre/simulation/50_lb/05_2tlb/simple_setup_test/output'#'./local_test/'
os.system('mkdir '+outconfpath)

pidf, typef, xposf, yposf, zposf, xvelf, yvelf, zvelf, Lxf, Lyf, Lzf = espressopp.tools.readxyz(inconffile)
scalefactor = int(Lxf)/Lxf
print scalefactor
box = (int(Lxf), int(Lyf), int(Lzf))

system         = espressopp.System()
system.rng     = espressopp.esutil.RNG()
system.bc      = espressopp.bc.OrthorhombicBC(system.rng, box)
system.skin    = skin
nodeGrid       = espressopp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
cellGrid       = espressopp.tools.decomp.cellGrid(box, nodeGrid, rc, skin)
system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)
  
integrator     = espressopp.integrator.VelocityVerlet(system)  
integrator.dt  = dt

#SIMPLE LANGEVIN
#thermostat = espressopp.integrator.LangevinThermostat(system)
#integrator.addExtension(thermostat)
#LANGEVIN DOUBLE THERMOSTAT
#thermostat             	= espressopp.integrator.LangevinThermostat2TId(system)
#thermostat.gamma1       = gamma
#thermostat.gamma2       = gamma
#thermostat.temperature1 = temperature1
#thermostat.temperature2 = temperature2
#thermostat.type1				= 0
#thermostat.type2				= 1
#thermostat.mpc = monomers_per_chain
#integrator.addExtension(thermostat)

fid = 'nc'+str(num_chains)+'N'+str(monomers_per_chain)+'dT'+str(temperature2)+'g'+str(gamma)+'LB'

#PRINTING PARAMETERS TO CHECK
print 'PARAMETERS OF THE RUN'
print '---------------------'
os.system('hostname')
print '{:22}{:25}'.format('infile:', inconffile)
print '{:22}{:25}'.format('num_chains:', num_chains)
print '{:22}{:25}'.format('monomers per chain:', monomers_per_chain)
print '{:22}{:25}'.format('temperature 2:', temperature2)
print '{:22}{:25}'.format('k_theta:', ktheta)
#print '{:22}{:25}{:30}'.format('thermostat gammas:', thermostat.gamma1, thermostat.gamma2)
print '{:22}{:25}'.format('time step:', dt)
print '{:22}{:25}'.format('fid:', fid)
print '{:22}{:25}'.format('outpath:', outconfpath)

#SETTING UP CONFORMATION
mass     = 1.0  
props    = ['id', 'type', 'mass', 'pos', 'v']
bondlist = espressopp.FixedPairList(system.storage)
if bending:
	anglelist = espressopp.FixedTripleList(system.storage)

for i in range(num_chains):
	chain = []
	bonds = []
	angles = []
	chtype = 0   #all chains of the same type - thermostat is using id to distinguish hot and cold
	for k in range(monomers_per_chain):
		idx  =  i * monomers_per_chain + k
		part = [ pidf[idx],  chtype, mass, espressopp.Real3D(xposf[idx]*scalefactor,yposf[idx]*scalefactor,zposf[idx]*scalefactor), espressopp.Real3D(xvelf[idx],yvelf[idx],zvelf[idx])]
		chain.append(part)
		if k>0:
			bonds.append((pidf[idx-1], pidf[idx]))
		if bending and k>1:
			angles.append((pidf[idx-2], pidf[idx-1], pidf[idx]))
	system.storage.addParticles(chain, *props)
	system.storage.decompose()
	bondlist.addBonds(bonds)
	if bending:
		anglelist.addTriples(angles)

#INTERACTIONS SETUP
#same LJ interaction between all particles
#LJ nonbonded
vl 						 = espressopp.VerletList(system, cutoff=rc)
interLJ00		     = espressopp.interaction.VerletListLennardJones(vl)
potLJ00			 		 = espressopp.interaction.LennardJones(epsilon = 1.0, sigma = 1.0, cutoff = rc,shift = 'auto')
interLJ00.setPotential(type1=0, type2=0, potential = potLJ00)
system.addInteraction(interLJ00)

# FENE bonds
potFENE   = espressopp.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
interFENE = espressopp.interaction.FixedPairListFENE(system, bondlist, potFENE)
system.addInteraction(interFENE)

# LATTICE BOLTZMANN SOLVENT
lb = espressopp.integrator.LatticeBoltzmann(system, nodeGrid)
integrator.addExtension(lb)

lb.nSteps = 10      # time scale contrast b/w LB and MD (= \delta t_LB / \delta t_MD)
lb.visc_b = 4.      # bulk viscosity (LJ-units)
lb.visc_s = 4.      # shear viscosity (LJ-units)
lb.lbTemp = temperature1 # cold temperature (LJ-units)
lb.fricCoeff = 20.  # coupling friction coefficient
lb.profStep = 1000  # profiling frequency (shows info on momentum conservation on the screen)

lb.highTemp = temperature2 # hot temperature
lb.chainLenMD = monomers_per_chain # monomers per chain to sort cold from hot

# set initial populations
initPop = espressopp.integrator.LBInitPopUniform(system,lb)
initPop.createDenVel(1.0, Real3D(0.,0.,0.))

# output files handling
lboutputScreen = espressopp.analysis.LBOutputScreen(system,lb)
OUT=espressopp.integrator.ExtAnalyze(lboutputScreen,lb.profStep)
integrator.addExtension(OUT)

#bending interation - only if bending switched on
# Cosine with FixedTriple list
if bending:
	potCosine = espressopp.interaction.Cosine(K=ktheta, theta0=0.0)
	interCosine = espressopp.interaction.FixedTripleListCosine(system, anglelist, potCosine)
	system.addInteraction(interCosine)


intime = time.time()
tsamp = 2000 #tau

t=tin
espressopp.tools.fastwritexyz(outconfpath+'prc_'+fid+'t'+str(int(t))+'.dat',system, velocities=True, unfolded = True)# mywrite writes the type as particleId/monomers_per_chain %2
nsteps = 2000#*int(thermostat.gamma1) #just a large number
tausteps = int(1.0/dt)

#HERE WE SHOULD ADD A SHORT WARMUP/EQUILIBRATION RUN DUE TO OTHER BOX SIZE
# it proly enough to run for 1600\tau
#integrator.run(monomers_per_chain*monomers_per_chain*tausteps)



for j in range(nsteps): 		#run nsteps x samp time
	espressopp.tools.analyse.info(system,integrator)
	integrator.run(tsamp*tausteps) 		#run for samp time * \tau
	t = integrator.step*integrator.dt + tin
	espressopp.tools.myfastwritexyz(outconfpath+'prc_'+fid+'t'+str(int(t))+'.dat',system, velocities=True, unfolded = True, mpc = monomers_per_chain)
        espressopp.tools.pdb.pdbwrite("traj.pdb", system, append=True)
        lb.saveLBConf() # dumps LB configuration for restart (for this set initial time step to the last one from previous simulation)

fintime = time.time()
espressopp.tools.analyse.final_info(system, integrator, vl, intime, fintime)

print "Done."
print "Running time: ", fintime-intime











