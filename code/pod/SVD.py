#from IPython import get_ipython
#get_ipython().magic('reset -sf')

import scipy
import scipy.linalg
import scipy.signal
from matplotlib import pyplot as plt
from pylab import genfromtxt
import os
# from termcolor import colored
import numpy as np
import numpy.fft as fft
import sys


# *************************************  Inputs  *************************************

nsnap = int(genfromtxt('./SVD_files/nsnap'))

# Load inputs file
with open('inputs.inp') as f:
    inputs = f.read()

# Split it into lines
inputs = inputs.splitlines()

# Find where POD INPUTS key word is
for i,lines in enumerate(inputs):
    if lines=='*POD INPUTS*': break
    
# Remove all unwanted content
inputs = [line.replace(' ','').split('!')[0] for line in inputs[i+1:] if line!= '']
        
# Assign values to variables        
nmodes = int(inputs[0])
if (nsnap < nmodes):
  nmodes = nsnap 
half_filter_list = [int(float(inputs[1])/2*nsnap)]

if (inputs[2]=='T' or inputs[2]=='t'): periodicity = 1
elif (inputs[2]=='F' or inputs[2]=='f'): periodicity = 0
else: raise ValueError('  INVALID PERIODICITY INPUT')  

if (inputs[3]!='G' and inputs[3]!='B'): raise ValueError ('   INVALID FILTER TYPE INPUT')
filter_list = inputs[3]

if (inputs[4]!='K' and inputs[4]!='P' and inputs[4]!='V'): raise ValueError ('   INVALID NORM INPUT')
norm_list = inputs[4]

dt = float(inputs[5])


# *************************************  Main code  *************************************

matrix = scipy.zeros((nsnap,nsnap),dtype=float)
os.system('mkdir -p figs/')
os.system('rm -rf figs/*')

# There is a loop here in case more than one filter is used
for filter_name in filter_list:

  for norm in norm_list:

    if norm == 'P': print  ('    Pressure norm')
    if norm == 'K': print  ('    Kinetic energy norm')
    os.system('mkdir -p figs/' + norm + '_norm/')
    
    # There is a loop here in case more than one filter size is used
    for Nf in half_filter_list:
      
      char_nf = '%03d' % Nf

      if norm == 'P': matrix[:,:] = genfromtxt('./SVD_files/corr_matrix_p.dat')
      if norm == 'K': matrix[:,:] = genfromtxt('./SVD_files/corr_matrix_k.dat')
      
      # Save full correlation matrix without any filtering
      plt.clf()   
      plt.contourf(matrix)
      plt.xlim([0,nsnap])
      plt.ylim([0,nsnap])
      plt.colorbar()
      label = norm + '_full_matrix_spod_Nf_000'
      plt.savefig('./figs/' + norm + '_norm/' + label + '.pdf')    


      #  Spectral POD stuff

      #~~~ Snapshot POD
      if Nf == 0:
        print ('    No filtering employed')

      #  Avoid shitwalking
      if 2*Nf+1 > nsnap:
        raise ValueError ('  Filter in spectral POD bigger than number of snapshots')
     
      #~~~ Filter set-up
      if Nf > 0:
        print ('     Spectral POD filter = ' + str( Nf))
        print ('     Total filter width = (2N+1) = ' +  str( 2*Nf+1))
        
        # Gaussian filter
        if filter_name == 'G':    
          print ('    -> Filter : Gaussian')
          sigma = Nf/4.0
          g = scipy.zeros(2*Nf+1)
          k = scipy.zeros(2*Nf+1)
          for i in range (0,2*Nf+1):
            k[i] = i-Nf
            g[i] = scipy.exp(-0.5*(k[i]/sigma)**2)
          
        # Boxcar filter    
        if filter_name == 'B':
          print ('    -> Filter : Boxcar')
          g = scipy.signal.boxcar(2*Nf+1,sym=True)
          
      # Filtering
      if Nf > 0:
        if periodicity == 0:
          print ('    -> No periodicity.. Stacking only zeros!')  
          mzeros = scipy.zeros((Nf,Nf))
          dummy1 = scipy.hstack((mzeros,mzeros,mzeros))
          dummy2 = scipy.hstack((mzeros,matrix,mzeros))
          dummy3 = scipy.hstack((mzeros,mzeros,mzeros))
          
          stack_matrix = scipy.vstack((dummy1,dummy2,dummy3))
        if periodicity == 1:
          print ('    -> Using periodicity.. Stacking the matrix itself!')
          dummy1 = matrix[-Nf:,-Nf:]
          dummy2 = matrix[-Nf:,:]
          dummy3 = matrix[-Nf:,:Nf]
          stack_matrix = np.hstack((dummy1,dummy2,dummy3))
          
          dummy1 = matrix[:,-Nf:]
          dummy3 = matrix[:,:Nf]
          
          stack_matrix = np.vstack( (stack_matrix,np.hstack( (dummy1,matrix,dummy3) ) ) )
          
          dummy1 = matrix[:Nf,-Nf:]
          dummy2 = matrix[:Nf,:]
          dummy3 = matrix[:Nf,:Nf]
          
          stack_matrix = np.vstack( (stack_matrix,np.hstack( (dummy1,dummy2,dummy3) ) ) )
          
        
        for i in range(Nf,Nf+nsnap):
          for j in range(Nf,Nf+nsnap):

            dummy = 0.0
            for l in range(0,2*Nf+1):
              k = l-Nf # Goes from -Nf to Nf+1
              dummy = dummy + stack_matrix[i+k,j+k]*g[l]

            matrix[i-Nf,j-Nf] = dummy/sum(g[:])

        plt.clf() 
        plt.contourf(matrix)
        plt.colorbar()
        label = norm + '_full_matrix_spod_Nf_' + str(char_nf) + '_filter_' + filter_name
        plt.savefig('./figs/' + norm + '_norm/' + label + '.pdf')    
        


# SVD 
V,S,_ = scipy.linalg.svd(matrix, full_matrices=1, compute_uv=1)
S = scipy.sqrt(S)
print ('    Saving temporal modes and eigenvalues ...')

file = open('./SVD_files/sigma.dat',"w")
file.write(norm  + ' \n')
file.write(str(nmodes)  + ' \n')
for x in S[:nmodes]:
  file.write('%.25E'  % x + ' \n') 
file.close()

file = open('./SVD_files/V_matrix.dat',"w")
file.write(norm  + ' \n')
file.write(str(nsnap) +  ' \n')
file.write(str(nmodes)  + ' \n')
file.write(str(dt)  + ' \n')
for x in V[:nsnap]:
  for y in x[:nmodes]:
    file.write('%.15E   '  % y )
    
  file.write(' \n') 
  
file.close()

# ****************************** PLOT SECTION ***************************************
try:
  import seaborn as sns
  sns.set_style("darkgrid")
except: 
  pass
print ('    Saving figures ...')

xf = np.fft.fftfreq(nsnap, d = dt)[1:]*2*np.pi
xf = [x for x in xf if (x<50) and (x>0)]
t  = [ dt*x for x in range(0,nsnap)]

for i in range (0,nmodes):	

  # Save temporal modes
  plt.clf() 
  plt.plot(t,V[:,i]*S[i])
  plt.grid(which = 'minor')
  char = '%02d' % (i+1)
  plt.savefig('./figs/eigenvec_' + char + '.png',dpi=450)
  
  # Save fft of temporal modes
  plt.clf()
  Vf = np.abs( fft.fft(V[:,i]) )[1:np.size(xf)+1]
  Vf = Vf / np.max(Vf)
  
  plt.semilogy(xf,Vf, label = 'Frequency of highest peak: %.4f' % (xf[np.argmax(Vf)]))
  plt.ylim(np.amax(Vf)*1e-4,np.amax(Vf)+0.1)
  plt.legend()
  plt.savefig('./figs/eigenvec_fft_' + char + '.png',dpi=450)

for i in range (0,nmodes-1):	
  # Save pairing of temporal modes
  plt.clf() 
  plt.plot(V[:,i],V[:,i+1])
  char1 = '%02d' % (i+1)
  char2 = '%02d' % (i+2)
  plt.savefig('./figs/modes_' + char1 + '_' + char2 + '.png',dpi=450)

# Distribution of singular values
S = S / np.sum(S)
SSum = 0.
SSum = [np.sum(S[:i+1]) for i,_ in enumerate(S)]
threshold = 0.80
i = 0
sumS = 0.
while (sumS < threshold):
  sumS = sumS + S[i]
  i = i + 1

print ('    First %.i modes (%.2f %%) represent %.2f %% of singular values spectrum' % (nmodes,nmodes/float(nsnap)*100.,np.sum(S[:nmodes]*100.)))
print ('    %.i modes (%.2f %%) are needed to represent %.2f %% of dataset' % (i,i/float(nsnap)*100.,sumS*100.))

xaxis = np.linspace(1,nsnap,nsnap)/nsnap

plt.figure(figsize=[15,10])

plt.subplot(2,3,1)
plt.title('Normalized singular values distribution')
plt.plot(xaxis[:nsnap-1],S[:nsnap-1],'ro-')
plt.xlabel('Percentage of files')

plt.subplot(2,3,2)
plt.title('Normalized singular values distribution (semilog)')
plt.semilogy(xaxis[:nsnap-1],S[:nsnap-1],'bo-')

plt.subplot(2,3,3)
plt.title('Sum of singular values')
plt.plot(xaxis,SSum,'go-')

plt.subplot(2,3,4)
plt.plot(xaxis[:nmodes],S[:nmodes],'ro-')

plt.subplot(2,3,5)
plt.semilogy(xaxis[:nmodes],S[:nmodes],'bo-')

plt.subplot(2,3,6)
plt.plot(xaxis[:nmodes],SSum[:nmodes],'go-')


plt.savefig('./figs/eigenval.png',dpi=750)


print ('    SVD complete!')
