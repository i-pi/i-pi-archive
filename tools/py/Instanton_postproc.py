#!/usr/bin/env python2

""" ins_Q_xml.py

Reads all the information needed from the RESTART file and
compute partition functions of the reactant, TS or instanton.

Syntax:
            python ins_Q_xml.py <name_input> <mode>   <temperature (K)>  ( <nbeads>)
    
  Example: python inst_Q_xml.py   RESTART  instanton       300
           python inst_Q_xml.py   RESTART  reactant        300               50
           python inst_Q_xml.py   RESTART  TS              300

"""

#04Oct17

import numpy as np
import sys
import os
import math
import time

#I-PI path
ipi_path='/home/jeremy/projects/i-pi-mc'
#ipi_path='/home/litman/Yair/Instanton/I-PI-mc/online-repo/i-pi-mc'

if not (os.path.exists(ipi_path)):
   print 'We can not find ipi in %s' %ipi_path
   print 'Please correct the path'
   sys.exit()

sys.path.insert(0, ipi_path)
from ipi.engine.simulation import Simulation

def ch4_filter(pos,dynmat,natoms,m,asr):
    pos= pos[:,:-3]
    dynmat= dynmat[:-3,:-3]
    natoms = natoms-1
    m =m[:-1]
    asr='poly'
    return  pos,dynmat,natoms,m,asr 

def get_rp_freq(w0,nbeads,temp,asr=None):
    """ Compute the ring polymer frequencies for an 3D harmonic potential
	defined by the frequencies w0. """
    hbar=1.0
    kb=1
    betaP=1/(kb*nbeads*temp) 
    factor= (betaP*hbar)
    
    w=1.0
    if np.amin(w0)<0.0:
       print '@get_rp_freq: We have a negative frequency, something is going wrong.'
       sys.exit()

    if asr =='poly':
      nzero = 6        
    elif asr =='crystal':
      nzero = 3        
    else:
      nzero = 0

    for i in range(nzero):
      for k in range(1,nbeads):
       w*=factor*np.sqrt(4./(betaP*hbar)**2 * np.sin(np.absolute(k)*np.pi/nbeads)**2)
    #Yes, for each K w is nbeads
    
    for n in range(w0.size):
     for k in range(nbeads):
        w*=factor*np.sqrt(4./(betaP*hbar)**2 * np.sin(np.absolute(k)*np.pi/nbeads)**2+w0[n])
        #note the w0 is the eigenvalue ( the square of the frequency )
    return w
 
def spring_hessian(nbeads,natoms,omega2,m):
    """ Add spring terms to the extended hessian
        """
    ii = natoms * 3
    h=np.zeros([ii*nbeads,ii*nbeads])

    if nbeads == 1:
       return h

    m3 = get_m3(natoms,nbeads,m)
    h_sp = m3[0] * omega2

    # Diagonal
    diag2 = np.diag(2.0 * h_sp)

    for i in range(0, nbeads):
        h[i * ii:(i + 1) * ii, i * ii:(i + 1) * ii] += diag2

    # Non-Diagonal
    ndiag = np.diag(-h_sp)
    #quasi-band
    for i in range(0, nbeads - 1):
        h[i * ii:(i + 1) * ii, (i + 1) * ii:(i + 2) * ii] += ndiag
        h[(i + 1) * ii:(i + 2) * ii, i * ii:(i + 1) * ii] += ndiag

    # Corner
    h[ 0 :ii , (nbeads-1)*ii:(nbeads) * ii   ] += ndiag
    h[(nbeads - 1) * ii:(nbeads) * ii, 0:ii  ] += ndiag
    return h

def get_doble(q0,nbeads0,natoms,h0):
    """Takes nbeads, positions and hessian (only the 'physcal part') of the half polymer and 
       returns the equivalent for the full ringpolymer."""
    q=np.concatenate((q0, np.flipud(q0)), axis=0)
    nbeads = 2*nbeads0
    ii     = 3*natoms
    iii    = 3*natoms*nbeads0

    h = np.zeros((iii*2, iii*2))
    h[0:iii,0:iii]=h0

    #diagonal block
    for i in range(nbeads0):
        x = i*ii+iii
        y = ((nbeads0-1)-i)*ii
        h[x:x+ii,x:x+ii] = h0[y:y+ii,y:y+ii]

    return q,nbeads,h

def flo(st):
    #st=str(st)
    st=st.replace("+","")
    return  np.float(st)

def get_m3(natoms,nbeads,m):
    """Computes 'm3' mass vector, one value of mass for each degree of freedom"""
    m3 = np.zeros((nbeads,3*natoms),float)
    m3[:,0:3*natoms:3] = m
    m3[:,1:3*natoms:3] = m3[:,0:3*natoms:3]
    m3[:,2:3*natoms:3] = m3[:,0:3*natoms:3]
    return m3

def get_dynmat(h,q,natoms,nbeads,m):
    """Computes the dynmat"""
    m3 = get_m3(natoms,nbeads,m)  
    ii = natoms * nbeads
    mm = np.zeros((nbeads, natoms))
    for i in range(nbeads):
        mm[i] = m
    mm     = mm.reshape(ii)
    ism    = m3.reshape(ii * 3) ** (-0.5)
    ismm   = np.outer(ism, ism)
    dynmat = np.multiply(h, ismm)
    return dynmat

def red2comp(h, nbeads,natoms):
    """Takes the reduced physical hessian and construct the 'complete' one (all 0 included) """
    i = natoms * 3
    ii = 3*natoms*nbeads
    h0 = np.zeros((ii, ii), float)

    for j in range(nbeads):
        h0[j * i:(j + 1) * i, j * i:(j + 1) * i] = h[:, j * i:(j + 1) * i]
    return h0

def clean_dynmat(dynmat,q, natoms,nbeads,m,asr):
    """
        Removes the translations and rotations modes.
        IN  dynmat = dynamical matrix
            q      = positions
            natoms = number of atoms
            nbeads = number of beads
            m      = mass vector, one value for each atom
                     is returned or not. Defaults to False.
        OUT d      = non zero eigenvalues of the dynmatrix
            w      = dynmatrix with the external modes projected out"""

    #Set some util things
    ii = natoms * nbeads
    mm = np.zeros((nbeads, natoms))
    m3 = get_m3(natoms,nbeads,m)
    for i in range(nbeads):
        mm[i] = m
    mm     = mm.reshape(ii)
    ism    = m3.reshape(ii * 3) ** (-0.5)

    if asr =='none':
        hm=dynmat
    else:
       #Computes the centre of mass.
       com = np.dot(np.transpose(q.reshape((ii, 3))), mm) / mm.sum()
       qminuscom = q.reshape((ii, 3)) - com

       if asr =='poly':
           #Computes the moment of inertia tensor.
           moi  = np.zeros((3,3), float)
           for k in range(ii):
               moi-=np.dot(np.cross(qminuscom[k],np.identity(3)),np.cross(qminuscom[k],np.identity(3)))*mm[k]

           I,U=(np.linalg.eig(moi))
           R=np.dot(qminuscom,U)
           D=np.zeros((6,3*ii),float)

           #Computes the vectors along translations and rotations.
           #Translations
           D[0]=np.tile([1,0,0],ii)/ism
           D[1]=np.tile([0,1,0],ii)/ism
           D[2]=np.tile([0,0,1],ii)/ism
           #Rotations
           for i in range(3*ii):
               iatom=i/3
               idof=np.mod(i,3)
               D[3,i]=(R[iatom,1]*U[idof,2] - R[iatom,2]*U[idof,1])/ism[i]
               D[4,i]=(R[iatom,2]*U[idof,0] - R[iatom,0]*U[idof,2])/ism[i]
               D[5,i]=(R[iatom,0]*U[idof,1] - R[iatom,1]*U[idof,0])/ism[i]


           for k in range(6):
               D[k] = D[k] / np.linalg.norm(D[k])
            #Computes the transformation matrix.
           transfmatrix = np.eye(3*ii) - np.dot(D.T,D)
           hm = np.dot(transfmatrix.T,np.dot(dynmat,transfmatrix))

       elif asr == 'crystal':
           # Computes the vectors along translations.
           # Translations
           D = np.zeros((6, 3 * ii), float)
           D[0] = np.tile([1, 0, 0], ii) / ism
           D[1] = np.tile([0, 1, 0], ii) / ism
           D[2] = np.tile([0, 0, 1], ii) / ism

           for k in range(3):
               D[k] = D[k] / np.linalg.norm(D[k])
           # Computes the transformation matrix.
           transfmatrix = np.eye(3 * ii) - np.dot(D.T, D)
           hm = np.dot(transfmatrix.T, np.dot(dynmat, transfmatrix))

    ##Simmetrize to use linalg.eigh
    hmT = hm.T
    hm  = (hmT+hm)/2.0
    d, w = np.linalg.eigh(hm)


    #Count
    dd = np.sign(d) * np.absolute(d) ** 0.5 / (2 * np.pi * 3e10 * 2.4188843e-17) #convert to cm^-1
    #print dd[0:9]
    #Zeros
    condition = np.abs(dd) < 0.01 #Note that dd[] units are cm^1
    nzero= np.extract(condition, dd)

    if asr =='poly' and nzero.size != 6:
         print "Warning, we have %s symmetry but we have %d 'zero' frequencies" %(asr,nzero.size)
         print "Freq. in cm^-1:"
         print  nzero

    if asr =='crystal' and nzero.size != 3:
         print "Warning, we have %s symmetry but we have %d 'zero' frequencies" %(asr,nzero.size)
         print "Freq. in cm^-1:"
         print  nzero

    #Negatives
    condition = dd < -4.0 #Note that dd[] units are cm^1
    nneg = np.extract(condition, dd)
    print "We have %d 'neg' frequencies " % (nneg.size)

   # Now eliminate external degrees of freedom from the dynmatrix

    if nzero.size > 0 and (asr=='poly' or asr=='crystal'):
        if np.linalg.norm(nzero) > 0.01:
            print " Warning @Clean hessian: We have deleted %d 'zero' frequencies " % (nzero.size)
            print " but the norm is greater than 0.01 cm^-1.  This should not happen." % (nzero.size)

        d = np.delete(d,range(nneg.size,nneg.size+nzero.size))
        w = np.delete(w, range(nneg.size,nneg.size+nzero.size),axis=1)
        
        #dd = np.sign(d) * np.absolute(d) ** 0.5 / (2 * np.pi * 3e10 * 2.4188843e-17) #convert to cm^-1
        #print dd[0:9]


    if asr =='poly':
       return d,w,np.prod(I)
    else:
       return d,w, 1.0


#START
np.set_printoptions(precision=6, suppress=True, threshold=np.nan)


#I/O
inputt  = sys.argv[1]
case    = sys.argv[2]
temp    = float(sys.argv[3])/315774.66

if case not in list(['reactant','TS','instanton']):
   raise ValueError("We can not indentify the case. The valid cases are: 'reactant', 'TS' and 'instanton'")	

print ''
print 'We are ready to start'
print 'Reading %s ...' %inputt 
print '(This can take a while)'
time0=time.time()

#
simulation = Simulation.load_from_xml(inputt, custom_verbosity='low',request_banner=False)
#print simulation.__dict__
beads   = simulation.syslist[0].motion.beads.copy()
m       = simulation.syslist[0].motion.beads.m.copy()
nbeads  = simulation.syslist[0].motion.beads.nbeads
natoms  = simulation.syslist[0].motion.beads.natoms
if case == 'reactant':
    asr        = simulation.syslist[0].motion.asr
    nbeadsQ    = int(sys.argv[4])
else:
    action  = simulation.syslist[0].motion.action
    asr     = simulation.syslist[0].motion.hessian_asr
#cell    = simulation.syslist[0].motion.cell.copy()


if case != 'instanton' and nbeads >1:
   print 'Incompatibility between case and nbeads in %s.' %(inputt)
   print 'case %s' %case
   print 'beads %i' %nbeads
   sys.exit()


#Units. We use atomic units.
b2a  = 0.52917721
h2K  = 3.1668152e-06
kb   = 1.0
hbar = 1.0
amu  = 1822.8885
au2K = 315774.66
au2eV= 27.211383414215543

#Print information
time1=time.time()
print ''
print 'We have finished the reading in %f s.'%(time1-time0)
print 'We used %i beads in the calculation.' %nbeads
print 'We have %i atoms.' %natoms
print 'ASR mode is %s.' %asr
print ''

#Depending the case we read from the input file different things:
print ''
print 'We need to get/create the dynmat'
time1=time.time()
if case=='reactant':
    dynmat  = simulation.syslist[0].motion.dynmatrix.copy()
    pos     = beads.q
    #print 'Hey, we are using CH4_filter HERE'
    ###pos,dynmat,natoms,m,asr = ch4_filter(pos,dynmat,natoms,m,asr) 

elif case=='TS':
    pos     = beads.q
    hessian = simulation.syslist[0].motion.hessian.copy()
    dynmat  = get_dynmat(hessian,pos,natoms,nbeads,beads.m)
    temp2   = simulation.syslist[0].ensemble.temp

elif case=='instanton':
    hessian = simulation.syslist[0].motion.hessian.copy()
    mode    = simulation.syslist[0].motion.mode
    temp2   = simulation.syslist[0].ensemble.temp
    pots    = simulation.syslist[0].motion.old_u
    grads   = - simulation.syslist[0].motion.old_f
    
    if mode != 'rate':
       print 'This script only work for rates calculation'
       sys.exit()

    if np.absolute(temp-temp2)*au2K >2:
        print ' '
	print 'Mismatch between provided temperature and temperature in the calculation'
	sys.exit()
    print 'The instanton mode is %s' %mode
    print 'The temperature is %f K' %(temp*au2K)
    print 'The full ring polymer is made of %i' %(nbeads*2)

    h0      = red2comp(hessian,nbeads,natoms)  
    pos,nbeads,hessian2 = get_doble(beads.q,nbeads,natoms,h0)
    hessian = hessian2
      
    omega2  = (temp * nbeads * kb / hbar) ** 2
    spring  = spring_hessian(nbeads,natoms,omega2,m) 
    h       = np.add(hessian, spring)
    dynmat  = get_dynmat(h,pos,natoms,nbeads,m)


time2=time.time()
print 'We got the dynmat in %f s.'%(time2-time1)
print ' '

#We have read all the data so we can start
beta =1.0/(kb*temp)
betaP=1.0/(kb*(nbeads)*temp)
d,w,detI = clean_dynmat(dynmat,pos,natoms,nbeads,m,asr)               
print  "Final lowest ten frequencies (cm^-1)"
print   np.sign(d[0:10]) * np.absolute(d[0:10]) ** 0.5 / (2 * np.pi * 3e10 * 2.4188843e-17)  # convert to cm^-1


if case=='reactant':
    Qtras    = ( ( np.sum(m) ) / ( 2*np.pi*beta*hbar**2 ) )**1.5               #Compute Qtras
    #Qtras_rp = ( ( np.sum(m)*nbeadsQ ) / ( 2*np.pi*betaP*hbar**2 ) )**1.5      #Compute Qtras_rp (without freq factor) 

    if asr == 'poly':
       Qrot     = ( 8*np.pi*detI   / ( (hbar)**6 * (beta)**3 ))**0.5           #Compute Qrot
    #   Qrot_rp  = ( 8*np.pi*detI*nbeadsQ**3 / ( (hbar)**6 * (betaP)**3 ))**0.5 #Compute Qrot_rp (without freq factor)  
    else:
       Qrot    = 1.0
    #  Qrot_rp = 1.0

    ddd= (np.absolute((d)))  
    logQvib    = -np.sum( np.log( 2*np.sinh( (beta*hbar*np.sqrt(d)/2.0) )  ))   #Compute Qvib 
    w_rp       = get_rp_freq(d,nbeadsQ,temp)
    logQvib_rp = -np.log(w_rp)                                         #Compute Qvib_rp 
    #w_rp       = get_rp_freq(d,nbeadsQ,betaP,asr)
    #logQvib_rp_plus   = np.sum( np.log(w_rp) )                                #Compute Qvib_rp+extra_freq   

    print ''
    print ''
    print 'We are done'
    print 'Nbeads %i' %nbeadsQ
    print ''
    print 'Qtras: %f bohr^-3'      %(Qtras)
    print 'Qrot: %f'       %(Qrot)
    print 'logQvib: %f'    %(logQvib)
    print ''
    #print 'Qtras_rp: %f'   %(Qtras_rp/nbeadsQ**3)
    #print 'Qrot_rp: %f'    %(Qrot_rp/nbeadsQ**3)
    print 'logQvib_rp: %f' %(logQvib_rp)
    #print 'logQvib_rp+: %f' %(logQvib_rp_plus)

elif case=='TS':
    Qtras = ( ( np.sum(m) ) / ( 2*np.pi*beta*hbar**2 ) )**1.5    #Compute Qtras

    if asr == 'poly':
       Qrot = ( 8*np.pi*detI   / ( (hbar)**6 * (beta)**3 ))**0.5 #Compute Qrot
    else:
       Qrot = 1.0

    print 'Note: Deleted frequency for computing Qvib  %f cm^-1' % (np.sign(d[0]) * np.absolute(d[0]) ** 0.5 / (2 * np.pi * 3e10 * 2.4188843e-17)) 
    logQvib = -np.sum( np.log( 2*np.sinh( (beta*hbar*np.sqrt(np.delete(d,0))/2.0) )  ))   #Compute Qvib 


    print ''
    print ''
    print 'We are done'
    print ''
    print 'Qtras: %f bohr^-3'      %(Qtras)
    print 'Qrot: %f'       %(Qrot)
    print 'logQvib: %f'    %(logQvib)
    print temp2*au2K
    print 'Potential energy at TS:  %f ev, U/kbT %f'       %((action[0]*kb*temp2)*au2eV,action[0])  
    print 'TEST  %f ev'       %((action[0])*au2eV)  
    print ''
     
elif case=='instanton': 
    Qtras     = ( ( np.sum(m) ) / ( 2*np.pi*beta*hbar**2 ) )**1.5               #Compute Qtras
    #Qtras_rp = ( ( np.sum(m)*nbeadsQ ) / ( 2*np.pi*betaP*hbar**2 ) )**1.5      #Compute Qtras_rp (without freq factor) 

    if asr == 'poly':
       #Qrot_rp  = ( 8*np.pi*detI / ( (hbar)**6 * (betaP)**3 ))**0.5       #Compute Qrot_rp (without freq factor)   
       Qrot      = ( 8*np.pi*detI / ( (hbar)**6 * (betaP)**3 ))**0.5/(nbeads)**3    #Compute Qrot
    else:
       Qrot      = 1.0

    print ''    
    print 'Note: Deleted frequency for computing Qvib  %f cm^-1' % (np.sign(d[1]) * np.absolute(d[1]) ** 0.5 / (2 * np.pi * 3e10 * 2.4188843e-17)) 

    #Compute Qvib

    if asr !='poly':
        print 'Warning asr != poly'
        print 'First 10 eigenvalues'
        print  (np.sign(d[0:10]) * np.absolute(d[0:10]) ** 0.5 / (2 * np.pi * 3e10 * 2.4188843e-17))
        print 'You have to correct the Qvib, the given value is not OK'

    logQvib = -np.sum( np.log( betaP*hbar*np.sqrt(np.absolute(np.delete(d,1))) ))+6*np.log(nbeads)+np.log(nbeads) 
    print 'logQvib = -np.sum( np.log( betaP*hbar*np.sqrt(np.absolute(np.delete(d,1))) ))+6*np.log(nbeadsQ)+np.log(nbeadsQ) '
    BN     =  2*np.sum(beads.m3[1:,:]*(beads.q[1:,:] - beads.q[:-1,:])**2)
    BN_factor =np.sqrt(BN/(2*np.pi*betaP*hbar**2))
    

    print ''
    print ''
    print 'We are done'
    print 'Nbeads %i' %(nbeads)
    print '(diff only %i)' %(nbeads/2)
    print '1/(betaP*hbar): %f a.u.' %(1/(betaP*hbar))
    print 'BN %f     ' % BN
    print 'BN*N %f   ' % (BN*nbeads)
    print 'BN_factor %f     ' % BN_factor

    print 'Qtras: %f bohr^-3'      %(Qtras)
    print 'Qrot: %f'    %(Qrot)
    print 'log(Qvib*N): %f' %(logQvib)
    print 'S1/hbar %f'  %action[0] 
    print 'S2/hbar %f'  %action[1] 
    print 'S/hbar %f'   %(np.sum(action))

    import math
    Lambda = math.sqrt(2*math.pi*betaP)
    logprefactor = logQvib - np.log(nbeads) + 0.5*math.log(BN) - math.log(betaP*Lambda)
    print 'logprefactor = %g' % logprefactor
    print 'prefactor = %g' % math.exp(logprefactor)

    if 1: # new method
	 	N = nbeads
		V0 = 0 # or energy of reactants
		x = beads.q.reshape((N//2, -1, 3))
		V = pots
		f = x.shape[1] * 3
		from modules.trajectory import Trajectory
		mass = m.reshape((-1,1)) * np.ones(3)
		from modules.linalg import orthogonalize # note this is not scipy or np linalg, but Jeremy's

		# trajectory stuff
		print
		def cat(x):
			N0 = N//4
			return np.concatenate((x[N0::-1],x,x[-1:N0-1:-1]))
		V -= V0
		dVdx = grads
		d2Vdx2 = simulation.syslist[0].motion.hessian.copy().T
		d2Vdx2.shape = (N//2,f,f)

		print 'building trajectory...'
		traj = Trajectory(cat(x).reshape(N+1,f), cat(V), cat(dVdx).reshape(N+1,f), cat(d2Vdx2), betaP*hbar, m=mass.ravel(), massweight=True)
		print 'S', traj.S()
		d2Sdx2 = traj.d2Sdx2()
		d2Sdxdt = traj.d2Sdxdt()
		d2Sdt2 = traj.d2Sdt2()

		# rotate coordinates to q coords
		p = traj.dSdx()
		#print p
		p = (p[f:] - p[:f]) / 2
		#print p
		U = np.identity(f)
		U[:,0] = p
		U = orthogonalize(U)
		U = np.block([[U, np.zeros_like(U)], [np.zeros_like(U), U]])

		d2Sdq2 = np.dot(U.T, np.dot(d2Sdx2, U))
		d2Sdqdt = np.dot(d2Sdxdt, U)

		E = traj.dSdt()
		print 'E = %g' % E
		p = math.sqrt(2*(traj.V[0]-E)) # mass-weighted p
		print 'p = %g' % p, np.linalg.norm(traj.dSdx())
		print '1/beta = %e' % (1/beta)

		GY = traj.Gelfand_Yaglom()
		evals = np.sort(1/np.linalg.eigvals(betaP*GY))
		print 'Gelfand-Yaglom eigvals:', evals
		C = np.prod(evals.real)
		print 'Gelfand-Yaglom: C = %g' % C
		print 'Gelfand-Yaglom: C without zero-modes = %g' % (beta**6*C)

		prefactor = 1/math.sqrt(2*math.pi) * p * math.sqrt(-beta**6*C)
		d2SdQ2 = d2Sdq2[1:f,1:f] + d2Sdq2[1:f,f+1:] + d2Sdq2[f+1:,1:f] + d2Sdq2[f+1:,f+1:]
		print np.linalg.norm(d2SdQ2 - d2SdQ2.T)
		evals = np.linalg.eigvalsh(d2SdQ2)
		print 'eigvals d2SdQ2', evals
		prefactor /= math.sqrt(np.prod(evals[6:]))
		print 'prefactor from action derivs', prefactor

time2=time.time()
print ''
print 'Total time %f s.'%(time2-time1)

sys.exit()
