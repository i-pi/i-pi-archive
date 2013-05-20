"""Contains functions for doing the inverse and forward normal mode transforms.

Classes:
   nm_trans: Uses matrix multiplication to do normal mode transformations.
   nm_rescale: Uses matrix multiplication to do ring polymer contraction
      or expansion.

Functions:
   mk_nm_matrix: Makes a matrix to transform between the normal mode and bead
      representations.
   mk_rs_matrix: Makes a matrix to transform between one number of beads and
      another. Higher normal modes in the case of an expansion are set to zero.
"""

__all__ = ['nm_trans', 'nm_rescale', 'nm_fft', 'FFT_nm_trans', 'FFT_inv_nm_trans' ]

import numpy as np
import pyfftw

def mk_nm_matrix(nbeads):
   """Gets the matrix that transforms from the bead representation
   to the normal mode representation.

   If we return from this function a matrix C, then we transform between the
   bead and normal mode representation using q_nm = C . q_b, q_b = C.T . q_nm

   Args:
      nbeads: The number of beads.
   """

   b2nm = np.zeros((nbeads,nbeads))
   b2nm[0,:] = np.sqrt(1.0/nbeads)
   for i in range(1, nbeads/2+1):
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
   """Gets the matrix that transforms a path with nb1 beads into one with
   nb2 beads.

   If we return from this function a matrix T, then we transform between the
   bead and normal mode representation using q_2 = T . q_1

   Args:
      nb1: The initial number of beads.
      nb2: The final number of beads.
   """

   if (nb1 == nb2):
      return np.identity(nb1,float)
   elif (nb1 > nb2):
      b1_nm = mk_nm_matrix(nb1)
      nm_b2 = mk_nm_matrix(nb2).T

      #builds the "reduction" matrix that picks the normal modes we want to keep
      b1_b2 = np.zeros((nb2, nb1), float)
      b1_b2[0,0] = 1.0
      for i in range(1, nb2/2+1):
         b1_b2[i,i] = 1.0
         b1_b2[nb2-i, nb1-i] = 1.0
      if (nb2 % 2 == 0):
         #if the original number of beads is odd, and the new number is
         #even then the non-degenerate mode's contribution needs to be split
         #equally among the appropriate degenerate modes of the new ring polymer
         b1_b2[nb2/2, nb2/2] = 0.5
         b1_b2[nb2/2, nb1-nb2/2] = 0.5

      rs_b1_b2 = np.dot(nm_b2, np.dot(b1_b2, b1_nm))
      return rs_b1_b2*np.sqrt(float(nb2)/float(nb1))
   else:
      return mk_rs_matrix(nb2, nb1).T*(float(nb2)/float(nb1))


class nm_trans:
   """Helper class to perform beads <--> normal modes transformation.

   Attributes:
      _b2nm: The matrix to transform between the bead and normal mode
         representations.
      _nm2b: The matrix to transform between the normal mode and bead
         representations.
   """

   def __init__(self, nbeads):
      """Initializes nm_trans.

      Args:
         nbeads: The number of beads.
      """

      self._b2nm = mk_nm_matrix(nbeads)
      self._nm2b = self._b2nm.T

   def b2nm(self, q):
      """Transforms a matrix to the normal mode representation.

      Args:
         q: A matrix with nbeads rows, in the bead representation.
      """

      return np.dot(self._b2nm,q)

   def nm2b(self, q):
      """Transforms a matrix to the bead representation.

      Args:
         q: A matrix with nbeads rows, in the normal mode representation.
      """

      return np.dot(self._nm2b,q)


class nm_rescale:
   """Helper class to rescale a ring polymer between different number of beads.

   Attributes:
      _b1tob2: The matrix to transform between a ring polymer with 'nbeads1'
         beads and another with 'nbeads2' beads.
      _b2tob1: The matrix to transform between a ring polymer with 'nbeads2'
         beads and another with 'nbeads1' beads.
   """

   def __init__(self, nbeads1, nbeads2):
      """Initializes nm_rescale.

      Args:
         nbeads1: The initial number of beads.
         nbeads2: The rescaled number of beads.
      """

      self._b1tob2 = mk_rs_matrix(nbeads1,nbeads2)
      self._b2tob1 = self._b1tob2.T*(float(nbeads1)/float(nbeads2))

   def b1tob2(self, q):
      """Transforms a matrix from one value of beads to another.

      Args:
         q: A matrix with nbeads1 rows, in the bead representation.
      """

      return np.dot(self._b1tob2,q)

   def b2tob1(self, q):
      """Transforms a matrix from one value of beads to another.

      Args:
         q: A matrix with nbeads2 rows, in the bead representation.
      """

      return np.dot(self._b2tob1,q)



class nm_fft:
   """Helper class to perform beads <--> normal modes transformation
      using Fast Fourier transforms.

   Attributes:
      _b2nm: The matrix to transform between the bead and normal mode
         representations.
      _nm2b: The matrix to transform between the normal mode and bead
         representations.
   """

   def __init__(self, nbeads, natoms):
      """Initializes nm_trans.

      Args:
         nbeads: The number of beads.
      """

      self.a = pyfftw.n_byte_align_empty((nbeads, 3*natoms), 16, 'float32')
      self.b = pyfftw.n_byte_align_empty((nbeads//2+1, 3*natoms), 16, 'complex64')
      self.fft = pyfftw.FFTW(self.a, self.b, axes=(0,), direction='FFTW_FORWARD')
      self.ifft = pyfftw.FFTW(self.b, self.a, axes=(0,), direction='FFTW_BACKWARD')

      pass

   def b2nm(self, q):
      """Transforms a matrix to the normal mode representation.

      Args:
         q: A matrix with nbeads rows, in the bead representation.
      """

      self.a[:] = q
      self.fft()
      print self.b
      return FFT_nm_trans(q)

   def nm2b(self, nmq):
      """Transforms a matrix to the bead representation.

      Args:
         q: A matrix with nbeads rows, in the normal mode representation.
      """

      return FFT_inv_nm_trans(nmq)

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
      return temp_mat.real/np.sqrt(nbeads)

   nmodes = nbeads/2

   temp_mat /= np.sqrt(nbeads)
   qnm = np.zeros(q.shape)
   qnm[0,:] = temp_mat[0,:].real

   if nbeads % 2 == 0:
      temp_mat[1:-1,:] *= np.sqrt(2)
      (qnm[1:nmodes,:], qnm[nbeads:nmodes:-1,:]) = (temp_mat[1:-1,:].real, temp_mat[1:-1,:].imag)
      qnm[nmodes,:] = temp_mat[nmodes,:].real
   else:
      temp_mat[1:,:] *= np.sqrt(2)
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
      return np.fft.irfft(qnm*np.sqrt(nbeads), n=nbeads, axis=0)

   nmodes = nbeads/2
   odd = nbeads - 2*nmodes  # 0 if even, 1 if odd

   qnm_complex = np.zeros((nmodes+1, len(qnm[0,:])), complex)
   qnm_complex[0,:] = qnm[0,:]
   if not odd:
      (qnm_complex[1:-1,:].real, qnm_complex[1:-1,:].imag) = (qnm[1:nmodes,:], qnm[nbeads:nmodes:-1,:])
      qnm_complex[1:-1,:] /= np.sqrt(2)
      qnm_complex[nmodes,:] = qnm[nmodes,:]
   else:
      (qnm_complex[1:,:].real, qnm_complex[1:,:].imag) = (qnm[1:nmodes+1,:], qnm[nbeads:nmodes:-1,:])
      qnm_complex[1:,:] /= np.sqrt(2)

   qnm_complex *= np.sqrt(nbeads)

   return np.fft.irfft(qnm_complex, n=nbeads, axis=0)
