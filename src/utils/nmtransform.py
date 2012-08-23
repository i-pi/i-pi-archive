"""Contains functions for doing the inverse and forward normal mode transforms.

Classes:
   nm_trans: Uses matrix multiplication to do normal mode transformations.

Functions:
   FFT_nm_trans: Uses an FFT algorithm to do the normal mode transformation.
   FFT_inv_nm_trans: Uses an FFT algorithm to do the inverse normal mode
      transformation.
"""

__all__ = ['nm_trans', 'nm_rescale', 'FFT_nm_trans', 'FFT_inv_nm_trans' ]

import numpy as np

def mk_nm_matrix(nbeads):
   """ Gets the matrix that transforms normal modes to beads. q_nm = b2nm . q_beads. """

   b2nm = np.zeros((nbeads,nbeads))
   b2nm[0,:] = np.sqrt(1.0/nbeads)
   for i in range(1,nbeads/2+1):
      for j in range(nbeads):
         b2nm[i,j] = np.sqrt(2.0/nbeads)*np.cos(2*np.pi*j*i/float(nbeads))
   if (nbeads%2) == 0:
      b2nm[nbeads/2,0:nbeads:2] = np.sqrt(1.0/nbeads)
      b2nm[nbeads/2,1:nbeads:2] = -np.sqrt(1.0/nbeads)
   for i in range(nbeads/2+1, nbeads):
      for j in range(nbeads):
         b2nm[i,j] = np.sqrt(2.0/nbeads)*np.sin(2*np.pi*j*i/float(nbeads))
   return b2nm

def mk_rs_matrix(nb1, nb2):
   """ Gets the matrix that transforms a path with nb1 beads into one with nb2 beads. q_2 = T_12 q_1. """

   if (nb1>nb2):
      b1_nm=mk_nm_matrix(nb1)
      nm_b2=mk_nm_matrix(nb2).T

      # builds the "reduction" matrix that picks the normal modes we want to keep
      b1_b2=np.zeros((nb2, nb1), float)
      b1_b2[0,0]=1.0
      for i in range(1,nb2/2+1):
         b1_b2[i,i]=1.0; b1_b2[nb2-i,nb1-i]=1.0
      if (nb2%2==0):
         b1_b2[nb2/2,nb2/2]=0.5; b1_b2[nb2/2,nb2/2+1]=0.5;
      rs_b1_b2=np.dot(nm_b2,np.dot(b1_b2,b1_nm))
      return rs_b1_b2*np.sqrt(float(nb2)/float(nb1))
   else:
      return mk_rs_matrix(nb2, nb1).T*(float(nb2)/float(nb1))



class nm_trans:
   """Helper class to perform beads <--> normal modes transformation for vectors describing a path of nbeads beads. """
   def __init__(self, nbeads):
      self._b2nm=mk_nm_matrix(nbeads)
      self._nm2b=self._b2nm.T

   def b2nm(self, q):
      return np.dot(self._b2nm,q)

   def nm2b(self, q):
      return np.dot(self._nm2b,q)


class nm_rescale:
   """Helper class to rescale vectors describing a path of nbeads1 beads into those of a nbeads2 path -- and vice-versa. """
   def __init__(self, nbeads1, nbeads2):
      self._b1tob2=mk_rs_matrix(nbeads1,nbeads2)
      self._b2tob1=self._b1tob2.T*(float(nbeads1)/float(nbeads2))

   def b1tob2(self, q):
      return np.dot(self._b1tob2,q)

   def b2tob1(self, q):
      return np.dot(self._b2tob1,q)


#~ class nm_trans:
   #~ """Performs the normal mode transformation using a matrix multiplication.
#~
   #~ Holds helper routines for the normal mode transformations.
#~
   #~ Attributes:
      #~ cmatrices: A set of transformation matrices for the normal mode
         #~ transformations.
      #~ tmatrices: A set of transformation matrices for ring polymer contraction.
   #~ """
#~
   #~ def __init__(self, nbeads=None):
      #~ """Initialises nm_trans.
#~
      #~ Args:
         #~ beads: An optional number of beads to use to create the transformation
            #~ matrix.
      #~ """
#~
      #~ self.cmatrices = {}
      #~ self.tmatrices = {}
      #~ if nbeads is not None:
         #~ self.setup_transform(nbeads)
#~
   #~ def setup_transform(self, nbeads):
      #~ """Sets up matrices for normal-mode transformation.
#~
      #~ Args:
         #~ nbeads: The number of beads in the ring polymer.
      #~ """
#~
      #~ # Todo: optional Fourier transform?
#~
      #~ Cb2nm = np.zeros((nbeads,nbeads))
      #~ Cb2nm[0,:] = math.sqrt(1.0/nbeads)
      #~ for i in range(1,nbeads/2+1):
         #~ for j in range(nbeads):
            #~ Cb2nm[i,j] = math.sqrt(2.0/nbeads)*math.cos(2*math.pi*j*i/float(nbeads))
      #~ if (nbeads%2) == 0:
         #~ Cb2nm[nbeads/2,0:nbeads:2] = math.sqrt(1.0/nbeads)
         #~ Cb2nm[nbeads/2,1:nbeads:2] = -math.sqrt(1.0/nbeads)
      #~ for i in range(nbeads/2+1, nbeads):
         #~ for j in range(nbeads):
            #~ Cb2nm[i,j] = math.sqrt(2.0/nbeads)*math.sin(2*math.pi*j*i/float(nbeads))
#~
      #~ Cnm2b = Cb2nm.T.copy()
#~
      #~ self.cmatrices[nbeads] = (Cb2nm, Cnm2b)
#~
   #~ def setup_contract(self, nred, nbeads):
      #~ """Sets up matrices for ring polymer contraction.
#~
      #~ Args:
         #~ nred: The number of beads in the contracted ring polymer.
         #~ nbeads: The number of beads in the full ring polymer.
      #~ """
#~
      #~ try:
         #~ full = self.cmatrices[nbeads]
      #~ except KeyError:
         #~ self.setup_transform(nbeads)
         #~ full = self.cmatrices[nbeads]
      #~ try:
         #~ contracted = self.cmatrices[nred]
      #~ except KeyError:
         #~ self.setup_transform(nred)
         #~ contracted = self.cmatrices[nred]
#~
      #~ Tf2c = np.zeros((nred, nbeads))
      #~ for j in range(-nred/2+1, nred/2+1):
         #~ #taking the lower frequency normal modes only.
         #~ Tf2c += np.outer(contracted[1][:,j], full[0][j,:])
#~
      #~ Tc2f = Tf2c.T.copy()
      #~ #normalization to set contracted ring polymer to a different effecive
      #~ #temperature.
      #~ Tf2c *= math.sqrt(nred/float(nbeads))
      #~ Tc2f *= math.sqrt(nbeads/float(nred))
#~
      #~ self.tmatrices[nred] = (Tf2c, Tc2f)
#~
   #~ def b2nm(self, q):
      #~ """This performs the forward normal mode transformation.
#~
      #~ Args:
         #~ q: A 2 dimensional matrix in the bead representation. The first
            #~ dimension gives the different bead coordinates, and the second
            #~ the different degrees of freedom.
#~
      #~ Returns:
         #~ A matrix of the same shape as q, but in the normal mode representation.
      #~ """
#~
      #~ nbeads = len(q)
      #~ try:
         #~ return np.dot(self.cmatrices[nbeads][0], q)
      #~ except KeyError:
         #~ self.setup_transform(nbeads)
         #~ return np.dot(self.cmatrices[nbeads][0], q)
#~
   #~ def nm2b(self, qnm):
      #~ """This performs the forward normal mode transformation.
#~
      #~ Args:
         #~ qnm: A 2 dimensional matrix in the normal mode representation.
            #~ The first dimension gives the different bead coordinates, and
            #~ the second the different degrees of freedom.
#~
      #~ Returns:
         #~ A matrix of the same shape as qnm, but in the bead representation.
      #~ """
#~
      #~ nbeads = len(qnm)
      #~ try:
         #~ return np.dot(self.cmatrices[nbeads][1], qnm)
      #~ except KeyError:
         #~ self.setup_transform(nbeads)
         #~ return np.dot(self.cmatrices[nbeads][1], qnm)
#~
   #~ def contract(self, q, nred):
      #~ """This performs a ring polymer contraction.
#~
      #~ This is used in the ring polymer contraction scheme, to create an
      #~ appropriate smaller ring polymer from a larger one
#~
      #~ Args:
         #~ q: A 2 dimensional matrix in the bead representation. The first
            #~ dimension gives the full bead coordinates, and the second
            #~ the different degrees of freedom.
         #~ nred: The number of beads in the contracted ring polymer.
#~
      #~ Returns:
         #~ A matrix of the same number of atoms as q, but with a smaller
         #~ number of beads.
      #~ """
#~
      #~ nbeads = len(q)
      #~ if nred == nbeads:
         #~ return q
      #~ else:
         #~ try:
            #~ return np.dot(self.tmatrices[nred][0], q)
         #~ except KeyError:
            #~ self.setup_contract(nred, nbeads)
            #~ return np.dot(self.tmatrices[nred][0], q)
#~
   #~ def expand(self, q, nbeads):
      #~ """This performs the ring polymer expansion.
#~
      #~ Args:
         #~ q: A 2 dimensional matrix in the bead representation.
            #~ The first dimension gives the contracted bead coordinates,
            #~ and the second the different degrees of freedom.
         #~ nbeads: The number of beads in the full ring polymer.
#~
      #~ Returns:
         #~ A matrix with the same number of atoms as q, but with the
         #~ full number of beads.
      #~ """
#~
      #~ nred = len(q)
      #~ if nred == nbeads:
         #~ return q
      #~ else:
         #~ try:
            #~ return np.dot(self.tmatrices[nred][1], q)
         #~ except KeyError:
            #~ self.setup_contract(nred, nbeads)
            #~ return np.dot(self.tmatrices[nred][1], q)


def FFT_nm_trans(q):
   """Performs the normal mode transformation using FFT.

   Args:
      q: A 2 dimensional matrix in the bead representation. The first
         dimension gives the different bead coordinates, and the second
         the different degrees of freedom.

   Returns:
      A matrix of the same shape as q, but in the normal mode representation.
   """

   temp_mat = np.fft.rfft(q, axis=0)
   nbeads = len(q)
   if nbeads < 3:
      return temp_mat.real/math.sqrt(nbeads)

   nmodes = nbeads/2
   odd = nbeads - 2*nmodes  # 0 if even, 1 if odd

   temp_mat /= math.sqrt(nbeads)
   qnm = np.zeros(q.shape)
   qnm[0,:] = temp_mat[0,:].real

   if not odd:
      temp_mat[1:-1,:] *= math.sqrt(2)
      (qnm[1:nmodes,:], qnm[nbeads:nmodes:-1,:]) = (temp_mat[1:-1,:].real, temp_mat[1:-1,:].imag)
      qnm[nmodes,:] = temp_mat[nmodes,:].real
   else:
      temp_mat[1:,:] *= math.sqrt(2)
      (qnm[1:nmodes+1,:], qnm[nbeads:nmodes:-1,:]) = (temp_mat[1:,:].real, temp_mat[1:,:].imag)

   return qnm

def FFT_inv_nm_trans(qnm):
   """Performs the inverse normal mode transformation using FFT.

   Args:
      qnm: A 2 dimensional matrix in the normal mode representation. The first
         dimension gives the different normal mode coordinates, and the second
         the different degrees of freedom.

   Returns:
      A matrix of the same shape as qnm, but in the bead representation.
   """

   nbeads = len(qnm)
   if nbeads < 3:
      return np.fft.irfft(qnm*math.sqrt(nbeads), n=nbeads, axis=0)

   nmodes = nbeads/2
   odd = nbeads - 2*nmodes  # 0 if even, 1 if odd

   qnm_complex = np.zeros((nmodes+1, len(qnm[0,:])), complex)
   qnm_complex[0,:] = qnm[0,:]
   if not odd:
      (qnm_complex[1:-1,:].real, qnm_complex[1:-1,:].imag) = (qnm[1:nmodes,:], qnm[nbeads:nmodes:-1,:])
      qnm_complex[1:-1,:] /= math.sqrt(2)
      qnm_complex[nmodes,:] = qnm[nmodes,:]
   else:
      (qnm_complex[1:,:].real, qnm_complex[1:,:].imag) = (qnm[1:nmodes+1,:], qnm[nbeads:nmodes:-1,:])
      qnm_complex[1:,:] /= math.sqrt(2)

   qnm_complex *= math.sqrt(nbeads)

   return np.fft.irfft(qnm_complex, n=nbeads, axis=0)
