#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from numpy import shape, dot, eye, set_printoptions, diag, sum
import numpy as np
from numpy.linalg import svd
import pdb
from spacepy.toolbox import feq
from matplotlib.mlab import find


"""
Data assimilation functions used for several projects

Authors
-------
Josef Koller, H. Godinez

jkoller@lanl.gov
Los Alamos National Laboratory

Copyright ©2010 Los Alamos National Security, LLC.
"""

# -----------------------------------------------
# Data Assimilation class: enKF
# -----------------------------------------------
class enKF(object):
	def __init__(self, ensembles = 50):
		self.Nens = ensembles # number of ensemble members

		return

	# -----------------------------------------------
	def __str__(self):
		return '<enKF instance, Nens='+self.Nens+'>'
	__repr__ = __str__

	# -----------------------------------------------
	def analysis(self, A, Psi, Inn, HAp):
		"""
		analysis subroutine after code example in Evensen 2003
		this will take the prepared matrices and calculate
		the analysis most efficiently, A will be returned

		Parameters
		==========
		A : numpy array
			What do I do?
		Psi :

		Inn :

		HAp :

		Returns
		=======
		out :

		Examples
		========

		"""
		# number of observations and ensemble size
		nobs, nens = HAp.shape
		nrmin = min(nobs,nens)
		NL = A.shape[0]

		# exit if only one observations
		if nobs is 1:
			print '\nERROR: nobs = 1 not supported this time\n'
			sys.exit(1)

		# compute HAp + Psi
		ES = HAp + Psi
		#print 'ES', ES
		# compute SVD of ES = HAp+Psi
		U,sig,Vt = svd(ES)
		sig = sig*sig # convert sig to eigenvalues

		# compute number of significant eigvals
		sigsum = sum(sig[0:nrmin-1])
		sigsum1 = 0.0
		for i in range(nrmin):
			if sigsum1/sigsum < 0.999:
				sigsum1 += sig[i]
				sig[i] = 1./sig[i]
			else:
				sig[i:] = 0.0

		print np.shape(sig), np.shape(U)
		# compute X1, X2, X3
		X1 = (sig*U).T
		X2 = dot(X1,Inn)
		X3 = dot(U,X2)

		# numerical check: sum of each col in X5 should be 1
		#X5 = (I+X4), X4 = HAp.T*X3
		#print HAp.shape, HAp.T.shape, X3.shape
		#X5 = (eye(nens) + dot(HAp.T,X3))
		#print diag(X5), diag(X5,k=-1)
		#savetxt('X5',X5)
		#print sum(X5,axis=0)

		# if nobs is small then compute representer first
		Reps = dot(A,HAp.T)
		B = dot(Reps,X3);
		C = A + B

		#C[np.where(C<0)] = 1e-99
		if 2*NL*nobs > nens*(nobs+NL):
			print '\nERROR: add routine when nobs is larger\n'
			#sys.exit(1)

		# check A for negative values
		#print np.log10(np.mean(A,axis=1))

		return C

	# -----------------------------------------------
        def analysis_Evensen(self, A, Psi, Inn, HAp):
		"""
                analysis subroutine after code example in
                Evensen 2003 this will take the prepared
                matrices and calculate the analysis most
                efficiently, A will be returned

		Parameters
		==========
		A :

		Psi :

		Inn :

		HAp :

		Returns
		=======
		out :

		Examples
		========
		"""
                # A   = ensemble of model states
                # Psi = E = observation perturbations
                # Inn = D = innovations
                # HAp = S = HA'
		nobs, nens = HAp.shape
		nrmin = min(nobs,nens)
		NL = A.shape[0]

                # compute HAp + Psi
                ES = HAp + Psi

                # Compute SVD of HAp + Psi  ->  U and sig,
                # using Eispack
                [U,sig,V] = np.linalg.svd(ES)
                U = U[:,0:nrmin]

                # convert eigenvalues
                sig = sig**2

                # Compute number of significant eigenvalues
                sigsum = np.sum(sig)
                sigsum1 = 0.0
                nrsigma = 0
                for i in range(nrmin):
                   if (sigsum1/sigsum < 0.999):
                      nrsigma = nrsigma+1
                      sigsum1=sigsum1+sig[i]
                      sig[i]=1.0/sig[i]
                   else:
                      sig[i:nrmin]=0.0
                      break
                   #endif
                #endfor

                # Compute X1
		X1 = (sig*U).T

                # Compute X2=X1*D
                X2 = np.dot(X1,Inn)

                # Compute X3=U*X2
                X3 = np.dot(U,X2)

                # Compute representers Reps=A'*S^T
                Reps =np.dot(A,HAp.transpose())

                # Compute A=A+Reps*X3
                C = A + np.dot(Reps,X3)

                return C

	# -----------------------------------------------
        def analysis_oneobs(self, A, Psi, Inn, HAp):
		"""
                analysis subroutine for a single observations
                with the EnKF. This is a special case.

		Parameters
		==========
		A :

		Psi :

		Inn :

		HAp :

		Returns
		=======
		out :

		Examples
		========
		"""

                # A   = ensemble of model states
                # Psi = E = observation perturbations
                # Inn = D = innovations
                # HAp = S = HA'

                # debugging command
                #pdb.set_trace()

		nobs, nens = HAp.shape
		nrmin = min(nobs,nens)
		NL = A.shape[0]

                # compute HPH^t
                HPHt = np.sum(HAp**2)

                # compute R
                R = np.sum(Psi**2)

                # get A'
		N1 = np.ones( (nens,nens) )/nens
		I = np.eye(nens)
		Ap = np.dot(A,I-N1)

                # compute PHt
                PHt = np.dot(Ap,HAp.T)

                # compute K
                K = 1/(HPHt + R)*PHt

                # Compute A=A+K*Inn
                C = A + np.outer(K,Inn)

                return C

	# -----------------------------------------------
	def getHAprime(self, HA):
		"""
		calculate ensemble perturbation of HA
		HA' = HA-HA_mean


		Parameters
		==========
		HA :

		Returns
		=======
		out :

		Examples
		========
		"""
		m, nens = np.shape(HA)
		N1 = np.ones( (nens,nens) )/nens
		I = np.eye(nens)
		S = np.dot(HA,I-N1)
		return S

	# -----------------------------------------------
	def getInnovation(self, y, Psi, HA):
		"""
		compute innovation ensemble D'

		Parameters
		==========
		y :

		Psi :

		HA :

		Returns
		=======
		out :

		Examples
		========
		"""
		m = np.size(y)
		nens = self.Nens
		# broadcoast y over ensembles to get Y
		Yens = np.ones( (m,nens) ) * y[:,np.newaxis]
		D = Yens + Psi - HA
		return D

	# -----------------------------------------------
	def getperturb(self, model, y):
		"""
		compute perturbations of observational vector

		Parameters
		==========
		model :

		y :

		Returns
		=======
		out :

		Examples
		========
		"""
		m = np.size(y)
		nens = self.Nens
		relstd = model.PSDerr
		Psi = np.random.standard_normal((m,nens)) * y[:,np.newaxis] * relstd
		#Psi = np.random.standard_normal((m,nens)) * np.mean(y)*relstd
		#Psi = np.zeros((m,nens))
		return Psi

	# -----------------------------------------------
	def getHA(self, model, Lobs, A):
		"""
		compute HA provided L vector of observations
		and ensemble matrix A

		Parameters
		==========
		model :

		Lobs :

		A :

		Returns
		=======
		out :

		Examples
		========
		"""
		# number of observations
		m = np.size(Lobs)
		# number of gridpoint and ensembles
		NL, nens = np.shape(A)
		# extract the grid
		Lgrid = model.Lgrid

		# dimension of HA will be (m x NL)(NL x nens) = m x nens
		HA = np.zeros( (m,nens) )
		i = 0
		for Lpos in Lobs:
                        # closest grid point
			index = np.argmin(np.abs(Lpos - Lgrid))
			if np.size(index) == 0:
				print '\nERROR in getHA\n'
				sys.exit(1)
			HA[i,:] = A[index,:]
			i += 1
		return HA

	# -----------------------------------------------
	def getHPH(self,Lobs,Pfxx):
		"""

		Parameters
		==========
		Lobs :

		Pfxx :

		Returns
		=======
		out :

		Examples
		========
		"""
		# dimension(L) = nobs(t_i)
		# dimension(Pfxx) = (NL,NL)

		Lgrid = dd['model']['Lgrid']
		NL = dd['model']['ngrid']

		HP = np.zeros((len(Lobs),NL))
		PH = np.zeros((NL,len(Lobs)))
		HPH = np.zeros((NL,NL))

		idx = np.zeros(len(Lobs), dtype=int)

		for Lval,iL in zip(Lobs, range(len(Lobs))):
			try:
				idx[iL] = np.argmin(np.abs(Lval - Lgrid))[0]
			except:
				print '\nERROR in getHPH\n'
				sys.exit(1)
			HP[iL,:] = Pfxx[:,idx[iL]]

		HPH = HP[:,idx]
		PH = Pfxx[:,idx]
		return PH, HPH

	# -----------------------------------------------
	def add_model_error(self, model, A, PSDdata):
		"""
		this routine will add a standard error to the ensemble states

		Parameters
		==========
		model :

		A :

		PSDdata :

		Returns
		=======
		out :

		Examples
		========
		"""
		Lgrid = model.Lgrid
		dL = Lgrid[1]-Lgrid[0]
		nens = int(self.Nens)
		#print y, L
                # in units of L (should not be smaller than dL)
		radius = dL*0.9
		radius = 1

		# create L vector where random noise will be added
		erridx = np.array([False]*len(Lgrid))
		for Lval in PSDdata['Lstar']:
			L1 = np.max((Lval-radius, Lgrid[0]))
			L2 = np.min((Lval+radius, Lgrid[-1]))
			Lidx = np.where( (Lgrid > L1) & (L2 > Lgrid) )[0]
			erridx[Lidx] = True

                # calculate the rel. average uncertainty based on (obs -
                # modelfcst)/obs
		fcst = np.mean(A, axis=1)
		rst = np.zeros( len(PSDdata['Lstar']) )
		for Lval, yval, i in zip(PSDdata['Lstar'], PSDdata['PSD'], range(len(PSDdata['Lstar'])) ):
			idx = np.argmin(np.abs(Lval - Lgrid))
			rst[i] = np.abs( yval - fcst[idx] )

		stdev = np.mean(rst)
		for iens in range(nens):
			f = A[erridx, iens]
			A[erridx, iens] = np.random.normal( f, scale = stdev )
		return A

	# -----------------------------------------------
	def add_model_error_obs(self, model, A, Lobs, y):
		"""
		this routine will add a standard error to the ensemble states

		Parameters
		==========
		model :

		A :

		Lobs :

		y :

		Returns
		=======
		out :

		Examples
		========
		"""
		Lgrid = model.Lgrid
		dL = Lgrid[1]-Lgrid[0]
		nens = int(self.Nens)
		#print y, L

                # in units of L (should not be smaller than dL)
		radius = dL*0.9
		radius = 1

		# create L vector where random noise will be added
		erridx = np.array([False]*len(Lgrid))
		for Lval in Lobs:
                    Lidx = np.where( feq(Lval,Lgrid) )[0]
                    erridx[Lidx] = True

                # calculate the rel. average uncertainty based on (obs -
                # modelfcst)/obs
		fcst = np.mean(A, axis=1)
		rst = np.zeros( len(Lobs) )
		for Lval, yval, i in zip(Lobs, y, range(len(Lobs)) ):
			idx = np.argmin(np.abs(Lval - Lgrid))
			rst[i] = np.abs( yval - fcst[idx] )

		stdev = np.mean(rst)
		for iens in range(nens):
                    A[erridx, iens] = np.random.normal( y, scale = stdev )
		return A

# -----------------------------------------------
# Data Assimilation class: enKF
# -----------------------------------------------
class direct(object):
	"""
	"""

	def __init__(self):
		return

	# -----------------------------------------------
	def __str__(self):
		"""
		"""
		return '<direct insertion instance>'
	__repr__ = __str__

	# -----------------------------------------------
	def analysis(self, model, PSDdata, HAp):
		"""
		analysis subroutine with simple direct insertion.

		Parameters
		==========
		model :

		PSDdata :

		HAp :

		Returns
		=======
		out :

		Examples
		========
		"""
		NL = model.NL
		#Nobs = len(PSDdata['PSD']

		# number of observations and ensemble size
		nobs, nens = HAp.shape
		nrmin = min(nobs,nens)
		NL = A.shape[0]

		# exit if only one observations
		if nobs is 1:
			print '\nERROR: nobs = 1 not supported this time\n'
			sys.exit(1)

		# compute HAp + Psi
		ES = HAp + Psi
		#print 'ES', ES
		# compute SVD of ES = HAp+Psi
		U,sig,Vt = svd(ES)
		sig = sig*sig # convert sig to eigenvalues

		# compute number of significant eigvals
		sigsum = sum(sig[0:nrmin-1])
		sigsum1 = 0.0
		for i in range(nrmin):
			if sigsum1/sigsum < 0.999:
				sigsum1 += sig[i]
				sig[i] = 1./sig[i]
			else:
				sig[i:] = 0.0

		print np.shape(sig), np.shape(U)
		# compute X1, X2, X3
		X1 = (sig*U).T
		X2 = dot(X1,Inn)
		X3 = dot(U,X2)

		# numerical check: sum of each col in X5 should be 1
		#X5 = (I+X4), X4 = HAp.T*X3
		#print HAp.shape, HAp.T.shape, X3.shape
		#X5 = (eye(nens) + dot(HAp.T,X3))
		#print diag(X5), diag(X5,k=-1)
		#savetxt('X5',X5)
		#print sum(X5,axis=0)

		# if nobs is small then compute representer first
		Reps = dot(A,HAp.T)
		B = dot(Reps,X3);
		C = A + B

		#C[np.where(C<0)] = 1e-99
		if 2*NL*nobs > nens*(nobs+NL):
			print '\nERROR: add routine when nobs is larger\n'
			#sys.exit(1)

		# check A for negative values
		#print np.log10(np.mean(A,axis=1))
		return C

# -----------------------------------------------
def average_window(PSDdata, Lgrid):
    """
    combine observations on same L shell in

    Parameters
    ==========
    model :

    PSDdata :

    HAp :

    Returns
    =======
    out :

    Examples
    ========
    """
    # sort observations first in L
    idx = PSDdata['Lstar'].argsort()
    Lobs = PSDdata['Lstar'][idx]
    y = PSDdata['PSD'][idx]

    # map Lobs onto closest grid point
    for i, obs in enumerate(y):
       Lobs[i] = Lgrid[ np.argmin(np.abs(Lobs[i]-Lgrid)) ]

    # identify unique grid-points of Lobs
    tmpLobs = np.unique(Lobs)

    # declare average observation array
    tmpy = np.zeros_like(tmpLobs)

    # run through all unique grid-points and compute average observation
    for i, iL in enumerate(tmpLobs):
        # identify idex of each unique grid-point
        idx = np.where(feq(iL,Lobs))
        # compute average observation for each unique grid-point
        tmpy[i] = np.average(y[idx])

    # walk through all grid points and find how many obs available
    #for i, iL in enumerate(Lgrid):
    #   idx = np.where( iL == Lobs[:])[0]

    # assign Lobs and y
    Lobs = tmpLobs
    y = tmpy
    return Lobs, y

# -----------------------------------------------
# -----------------------------------------------
# -----------------------------------------------
def getobs4window(dd, Tnow):
    """
    get observations in time window [Tnow - Twindow, Tnow]
    from all satellites lumped together into one y vector

    Parameters
    ==========
    model :

    PSDdata :

    HAp :

    Returns
    =======
    out :

    Examples
    ========
    """
    Tstart = Tnow - dd['model']['Twindow']
    Tend = Tnow

    # read all PSD observation and append into y,L
    y = []; L = []
    for satkey in dd['data']['satids']: # for every s/c
        PSD = dd['data'][satkey]['PSD']
        TAI = dd['data'][satkey]['TAI']
        Lvals = dd['data'][satkey]['L']
        indices = np.where( (Tstart < TAI) & (TAI <= Tend) )
        y = np.append(y, PSD[indices])
        L = np.append(L, Lvals[indices])

    # sort them and average if they fall into same L range
    ind = L.argsort() # return indices to a sorted y
    y = y[ind]
    L = L[ind]

    #check if observation available
    if len(y) == 0: return dd, np.array([]), np.array([])

    i = 0
    Lmean = []; ymean = []
    Lgrid = dd['model']['Lgrid']
    for iL in Lgrid[1:-1]:
        i += 1
        L_lower = np.mean([ Lgrid[i], Lgrid[i-1] ])
        L_upper = np.mean([ Lgrid[i+1], Lgrid[i] ])
        indices = np.where( (L_lower < L) & (L <= L_upper) )
        if np.size(indices) is 0: continue

        Lmean = np.append(Lmean, iL)
        ymean = np.append(ymean, np.mean(y[indices]))

    # add Lmean, ymean to dd
    try: # first time executing this
        dd['data']['dwindow'].append( (Tnow, Lmean, ymean) )
    except:
        dd['data']['dwindow'] = []
        dd['data']['dwindow'].append( (Tnow, Lmean, ymean) )

    # these line are only for averaging over L* as well!

    Lmean = np.mean(Lmean)
    # move to actual grid point
    print Lmean
    idx = find(Lgrid>Lmean)[0]

    Lmean = Lgrid[idx]
    ymean = np.mean(ymean)
    return dd, np.array([Lmean]), np.array([ymean])

# -----------------------------------------------
def output(init,result):
    """
    write results to file and be done

    Parameters
    ==========
    model :

    PSDdata :

    HAp :

    Returns
    =======
    out :

    Examples
    ========
    """
    return

# -----------------------------------------------
# ----- older stuff (not used) ------------------
# -----------------------------------------------
def assimilate_JK(dd):
    """
    this version is currently not working
    main function to assimilate all data provided in init

    Parameters
    ==========
    model :

    PSDdata :

    HAp :

    Returns
    =======
    out :

    Examples
    ========
    """
    np.random.seed(123)

    Lgrid = dd['model']['Lgrid']
    dd['model']['initPSD'] = np.array( 4e-7*np.exp( -(Lgrid-5)**2/2) )

    nens = dd['kalman']['nens']
    NL = dd['model']['ngrid']
    f0 = np.array(dd['model']['initPSD'])
    Tgrid = dd['model']['Tgrid']

    dd['results'] = {}
    dd['results']['fa'] = {}
    dd['results']['ff'] = {}

    # --- initialize some stuff
    # broadcast all f0 into A and randomize
    A = np.ones( (NL,nens) )*f0[:,np.newaxis]
    #A = 10**np.random.normal(np.log10(A),0.3)  # watch out for negative values

    dd['results']['fa']['Tgrid'] = np.zeros(len(Tgrid))
    dd['results']['ff']['Tgrid'] = np.zeros(len(Tgrid))
    dd['results']['fa']['PSD'] = np.zeros( (f0.shape[0],Tgrid.shape[0]) )
    dd['results']['fa']['PSD'][:,0] = f0
    dd['results']['fa']['Tgrid'][0] = Tgrid[0]
    dd['results']['ff']['PSD'] = np.zeros( (f0.shape[0],Tgrid.shape[0]) )
    dd['results']['ff']['PSD'][:,0] = f0
    dd['results']['ff']['Tgrid'][0] = Tgrid[1]

    minPSD = np.min(A)/100

    # time loop to assimilate everything
    for Tnow, tcnt in zip(Tgrid, range(len(Tgrid)-1)):

        # make forecast using all ensembles in A
        icnt = 0
        for f in A.T:
            fnew = forecast(dd, f, Tnow)
            A[:,icnt] = fnew[:]
            icnt += 1

        #print np.log10(np.mean(A,axis=1))
        Tnow = Tnow + dd['model']['Twindow']
        tcnt = tcnt + 1
        dd['results']['ff']['PSD'][:,tcnt] = np.mean(A, axis=1)
        dd['results']['ff']['Tgrid'][tcnt] = Tnow

        # get observations for time window ]Tnow-Twindow,Tnow]
        dd, L,y = getobs4window(dd, Tnow)

        # check if len(y) ==0
        if len(y) ==0:
            dd['results']['fa']['PSD'][:,tcnt] = np.mean(A, axis=1)
            dd['results']['fa']['Tgrid'][tcnt] = Tnow
            continue

        # perturb observations so that ensemble get sufficient spread
        HA = getHA(dd,L,A) # project ensemble states to obs. grid: HA(nobs,nens)
        e = np.zeros( (len(y),nens) )  # this is Psi in Evensen 2003
        D = np.zeros( (len(y),nens) )
        err_z = 0.3
        for yval, iobs in zip(y, range(len(y))):
            relstd = yval*err_z
            rnd = np.random.normal( yval*np.ones(nens), relstd)
            rnd[np.where(rnd<minPSD)] = minPSD
            e[iobs,:] = rnd - yval
            D[iobs,:] = yval + e[iobs,:]

        # add model error
        relstd = np.zeros(nens)
        for yval, iobs in zip(y, range(len(y))):
            idx = np.where( feq(L[iobs],Lgrid) )
            relstd[:] = 0.5*( yval - HA[iobs,:] )
            A[idx,:] = np.random.normal( HA[iobs,:], np.abs(relstd))
            idx2 = np.where(A[idx,:] < minPSD); A[idx,idx2] = minPSD # add min flux here

        # get residual or innovation in the classical sense
        Res = D - HA # dim = (nobs,nens)

        # error of covariant matrix
        ffmean = np.mean(A, axis=1)
        V = np.zeros((NL,nens))
        for iens in range(nens):
            V[:,iens]  = A[:,iens] - ffmean

        Pfxx = np.dot(V,V.T)/float(nens)
        Pfyy = np.dot(e,e.T)/float(nens)

        # get Kalman denominator (H*Pfxx*HT + Pfyy)^-1 by inversion
        PH, HPH = getHPH(dd,L,Pfxx)
        #print PH
        #print HPH
        Rinv = np.linalg.inv(HPH+Pfyy)
        Adj = np.dot(np.dot(PH,Rinv), Res)

        # update all ensemble members
        A = A + Adj

        # check for negative fluxes
        idx = np.where(A<minPSD)
        A[idx] = minPSD

        #
        # average A from analysis step and save in results dictionary
        #print y, np.mean(HA,axis=1), np.mean(getHA(init,L,A),axis=1)
        dd['results']['fa']['PSD'][:,tcnt] = np.mean(A, axis=1)
        dd['results']['fa']['Tgrid'][tcnt] = Tnow

        # print message
        print 'Tnow: ', rbtools.TAInum2ISOdate(Tnow)

        #print np.log10(np.mean(A,axis=1))[30]

    # copy some results into ['kalman'] key
    dd['kalman']['fa'] =  dd['results']['fa']
    dd['kalman']['ff'] =  dd['results']['ff']
    return dd

# -----------------------------------------------
def addmodelerror_old2(dd, A, y, L):
    """
    this routine will add a standard error to the ensemble states
    """
    from numpy import where, random, log10, mean, reshape, abs, max, arange
    from numpy import linspace, min

    Lgrid = dd['model']['Lgrid']
    dL = Lgrid[1]-Lgrid[0]
    nens = int(dd['kalman']['nens'])
    #print y, L

    radius = 1

    for Lcenter, yval in zip(L, y):
        L1 = max((Lcenter-radius, Lgrid[0]))
        L2 = min((Lcenter+radius, Lgrid[-1]))
        #print Lcenter, Lcenter+radius, Lgrid[-1], min((Lcenter+radius,Lgrid[-1])), (L2-L1)/dL
        NLs = int(round( (L2-L1)/dL ) + 1 )
        for Lpos in linspace(L1, L2, NLs):
            #print dL, L1, L2, NLs
            index = where( feq(Lpos,Lgrid) ) # use float point comparison
            stdev = 0.5*abs( yval - mean(A[index,:]) )
            #print Lpos
            center = reshape( A[index,:], (nens) )
            A[index,:] = random.normal( center, scale=stdev)
            # now check if any below 1e-99 and repeat
            for i in range(nens):
                icnt = 0
                while A[index,i] < 1e-99:
                    A[index,i] = random.normal( center[i], scale=stdev)
                    icnt += 1
                    if icnt > 1000: print 'too many iterations'
    return A

# -----------------------------------------------
def addmodelerror_old(dd, A, y, L):
    """
    this routine will add a standard error to the ensemble states
    """
    Lgrid = dd['model']['Lgrid']
    dL = Lgrid[1]-Lgrid[0]
    nens = int(dd['kalman']['nens'])
    #print y, L

    radius = 1

    for Lcenter, yval in zip(L, y):
        L1 = np.max((Lcenter-radius, Lgrid[0]))
        L2 = np.min((Lcenter+radius, Lgrid[-1]))
        #print Lcenter, Lcenter+radius, Lgrid[-1], min((Lcenter+radius,Lgrid[-1])), (L2-L1)/dL
        NLs = int(round( (L2-L1)/dL ) + 1 )
        for Lpos in np.linspace(L1, L2, NLs):
            #print dL, L1, L2, NLs
            index = np.where( feq(Lpos,Lgrid) ) # use float point comparison
            stdev = 0.5*np.abs( yval - np.mean(A[index,:]) )
            #print Lpos
            center = np.reshape( A[index,:], (nens) )
            A[index,:] = np.random.normal( center, scale=stdev)
            # now check if any below 1e-99 and repeat
            for i in range(nens):
                icnt = 0
                while A[index,i] < 1e-99:
                    #A[index,i] = np.random.normal( center[i], scale=stdev)
                    A[index,i] = 10**np.random.normal( np.log10(center[i]), scale=0.001)
                    icnt += 1
                    if icnt > 1000: print 'too many iterations'
    return A

# -----------------------------------------------
def forecast(dd, f, Tnow):
    """
    call the model to calculate update and return result
    fnew = f(Tnow+Twindow)
    """
    import diff1d
    fnew = diff1d.diff_LL(dd, f, Tnow)
    return fnew
