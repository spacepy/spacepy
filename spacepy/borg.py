#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
Data assimilation functions

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
	"""
	Version:
	========
	V1: 
	
	"""
	
	def __init__(self, ensembles = 50):
		self.Nens = ensembles # number of ensemble members
		
		return
		
	# -----------------------------------------------    
	def __str__(self):
		"""
		"""
 
		return '<enKF instance, Nens='+self.Nens+'>'
	__repr__ = __str__

	# -----------------------------------------------
	def analysis(self, A, Psi, Inn, HAp):
		"""
		analysis subroutine after code example in Evensen 2003
		this will take the prepared matrices and calculate
		the analysis most efficiently, A will be returned
		"""
	
		import sys
		from numpy import shape, dot, eye, set_printoptions, diag, sum
		#from numpy import shape, dot, eye, set_printoptions, diag, savetxt, sum
		import numpy as n
		from numpy.linalg import svd
	
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
	
		print n.shape(sig), n.shape(U)
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
		
		#C[n.where(C<0)] = 1e-99
		if 2*NL*nobs > nens*(nobs+NL):
			print '\nERROR: add routine when nobs is larger\n'
			#sys.exit(1)
	
		# check A for negative values
		#print n.log10(n.mean(A,axis=1))
	
		return C



	# -----------------------------------------------
        def analysis_Evensen(self, A, Psi, Inn, HAp):
		"""
                analysis subroutine after code example in
                Evensen 2003 this will take the prepared
                matrices and calculate the analysis most
                efficiently, A will be returned
		"""

                # A   = ensemble of model states
                # Psi = E = observation perturbations
                # Inn = D = innovations
                # HAp = S = HA'

                import numpy
		import sys
	
		nobs, nens = HAp.shape
		nrmin = min(nobs,nens)
		NL = A.shape[0]

                # compute HAp + Psi
                ES = HAp + Psi

                # Compute SVD of HAp + Psi  ->  U and sig,
                # using Eispack
                [U,sig,V] = numpy.linalg.svd(ES)
                U = U[:,0:nrmin]

                # convert eigenvalues
                sig = sig**2

                # Compute number of significant eigenvalues
                sigsum = numpy.sum(sig)
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
                X2 = numpy.dot(X1,Inn)

                # Compute X3=U*X2
                X3 = numpy.dot(U,X2)

                # Compute representers Reps=A'*S^T
                Reps =numpy.dot(A,HAp.transpose()) 

                # Compute A=A+Reps*X3
                C = A + numpy.dot(Reps,X3)

                return C

	# -----------------------------------------------
        def analysis_oneobs(self, A, Psi, Inn, HAp):

		"""
                analysis subroutine for a single observations
                with the EnKF. This is a special case.
		"""

                # A   = ensemble of model states
                # Psi = E = observation perturbations
                # Inn = D = innovations
                # HAp = S = HA'

                import numpy
		import sys
                import pdb
	
                # debugging command 
                #pdb.set_trace()

		nobs, nens = HAp.shape
		nrmin = min(nobs,nens)
		NL = A.shape[0]

                # compute HPH^t
                HPHt = numpy.sum(HAp**2)

                # compute R
                R = numpy.sum(Psi**2)

                # get A'
		N1 = numpy.ones( (nens,nens) )/nens
		I = numpy.eye(nens)
		Ap = numpy.dot(A,I-N1)

                # compute PHt
                PHt = numpy.dot(Ap,HAp.T)

                # compute K
                K = 1/(HPHt + R)*PHt

                # Compute A=A+K*Inn
                C = A + numpy.outer(K,Inn)

                return C

	# -----------------------------------------------
	def getHAprime(self, HA):
		"""
		calculate ensemble perturbation of HA
		HA' = HA-HA_mean
		"""
	
		import numpy as n
	
		m, nens = n.shape(HA)
	
		N1 = n.ones( (nens,nens) )/nens
		I = n.eye(nens)
	
		S = n.dot(HA,I-N1)
	
		return S
	
	
	# -----------------------------------------------
	def getInnovation(self, y, Psi, HA):
		"""
		compute innovation ensemble D'
		"""
	
		import numpy as n
	
		m = n.size(y)
		nens = self.Nens
	
		# broadcoast y over ensembles to get Y
		Yens = n.ones( (m,nens) ) * y[:,n.newaxis]
	
		D = Yens + Psi - HA
		
		return D
	
	# -----------------------------------------------
	def getperturb(self, model, y):
		"""
		compute perturbations of observational vector
		"""
	
		import numpy as n
	
		m = n.size(y)
		nens = self.Nens
		relstd = model.PSDerr
		
		Psi = n.random.standard_normal((m,nens)) * y[:,n.newaxis] * relstd
		#Psi = n.random.standard_normal((m,nens)) * n.mean(y)*relstd
		#Psi = n.zeros((m,nens))
		return Psi
	
	# -----------------------------------------------
	def getHA(self, model, Lobs, A):
		"""
		compute HA provided L vector of observations
		and ensemble matrix A
		"""
	
		import numpy as n
		from spacepy.toolbox import feq 
		import sys
	
		# number of observations
		m = n.size(Lobs)
		# number of gridpoint and ensembles
		NL, nens = n.shape(A)
		# extract the grid
		Lgrid = model.Lgrid
		
		# dimension of HA will be (m x NL)(NL x nens) = m x nens
		HA = n.zeros( (m,nens) )

		i = 0
		for Lpos in Lobs:        
                        # closest grid point
			index = n.argmin(n.abs(Lpos - Lgrid))

			if n.size(index) == 0:
				print '\nERROR in getHA\n'
				sys.exit(1)
			HA[i,:] = A[index,:]
			i += 1
	
		return HA
			
	# -----------------------------------------------
	def getHPH(self,Lobs,Pfxx):
	
		import numpy as n
	
		# dimension(L) = nobs(t_i)
		# dimension(Pfxx) = (NL,NL)
	
		Lgrid = dd['model']['Lgrid']
		NL = dd['model']['ngrid']
	
		HP = n.zeros((len(Lobs),NL))
		PH = n.zeros((NL,len(Lobs)))
		HPH = n.zeros((NL,NL))
	
		idx = n.zeros(len(Lobs), dtype=int)
		
		for Lval,iL in zip(Lobs, range(len(Lobs))):
			try:
				idx[iL] = n.argmin(n.abs(Lval - Lgrid))[0]
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
		"""
	
                # from numpy import where, random, log10, mean, reshape, abs,
                # max, arange
		# from numpy import linspace, min
		import numpy as n
		
		Lgrid = model.Lgrid
		dL = Lgrid[1]-Lgrid[0]
		nens = int(self.Nens)
		#print y, L
	
                # in units of L (should not be smaller than dL)
		radius = dL*0.9
		radius = 1
		
		# create L vector where random noise will be added
		erridx = n.array([False]*len(Lgrid))
		for Lval in PSDdata['Lstar']:
			L1 = n.max((Lval-radius, Lgrid[0]))
			L2 = n.min((Lval+radius, Lgrid[-1]))
			Lidx = n.where( (Lgrid > L1) & (L2 > Lgrid) )[0]
			erridx[Lidx] = True
	
                # calculate the rel. average uncertainty based on (obs -
                # modelfcst)/obs
		fcst = n.mean(A, axis=1)
		rst = n.zeros( len(PSDdata['Lstar']) )
		for Lval, yval, i in zip(PSDdata['Lstar'], PSDdata['PSD'], range(len(PSDdata['Lstar'])) ):
			idx = n.argmin(n.abs(Lval - Lgrid))
			rst[i] = n.abs( yval - fcst[idx] )
	
		stdev = n.mean(rst)
		for iens in range(nens):
			f = A[erridx, iens]
			A[erridx, iens] = n.random.normal( f, scale = stdev )

	
		return A


	# -----------------------------------------------
	def add_model_error_obs(self, model, A, Lobs, y):
		"""
		this routine will add a standard error to the ensemble states
		"""
	
		import numpy as n
		from spacepy.toolbox import feq
		
		Lgrid = model.Lgrid
		dL = Lgrid[1]-Lgrid[0]
		nens = int(self.Nens)
		#print y, L
	
                # in units of L (should not be smaller than dL)
		radius = dL*0.9
		radius = 1
		
		# create L vector where random noise will be added
		erridx = n.array([False]*len(Lgrid))
		for Lval in Lobs:
                    Lidx = n.where( feq(Lval,Lgrid) )[0]
                    erridx[Lidx] = True
	
                # calculate the rel. average uncertainty based on (obs -
                # modelfcst)/obs
		fcst = n.mean(A, axis=1)
		rst = n.zeros( len(Lobs) )
		for Lval, yval, i in zip(Lobs, y, range(len(Lobs)) ):
			idx = n.argmin(n.abs(Lval - Lgrid))
			rst[i] = n.abs( yval - fcst[idx] )
	
		stdev = n.mean(rst)
		for iens in range(nens):
                    A[erridx, iens] = n.random.normal( y, scale = stdev )
		return A

# -----------------------------------------------
# Data Assimilation class: enKF
# -----------------------------------------------    
class direct(object):
	"""
	Version:
	========
	V1: 
	
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
		"""
		import sys
		import numpy as n
		from numpy.linalg import svd
	
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
	
		print n.shape(sig), n.shape(U)
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
		
		#C[n.where(C<0)] = 1e-99
		if 2*NL*nobs > nens*(nobs+NL):
			print '\nERROR: add routine when nobs is larger\n'
			#sys.exit(1)
	
		# check A for negative values
		#print n.log10(n.mean(A,axis=1))
	
		return C

# -----------------------------------------------
def average_window(PSDdata, Lgrid):
    """
    combine observations on same L shell in 
    """
    import numpy as n
    from spacepy.toolbox import feq

    # sort observations first in L
    idx = PSDdata['Lstar'].argsort()
    Lobs = PSDdata['Lstar'][idx]
    y = PSDdata['PSD'][idx]

    # map Lobs onto closest grid point
    for i, obs in enumerate(y):
       Lobs[i] = Lgrid[ n.argmin(n.abs(Lobs[i]-Lgrid)) ]

    # identify unique grid-points of Lobs
    tmpLobs = n.unique(Lobs)

    # declare average observation array
    tmpy = n.zeros_like(tmpLobs)

    # run through all unique grid-points and compute average observation
    for i, iL in enumerate(tmpLobs):

        # identify idex of each unique grid-point
        idx = n.where(feq(iL,Lobs))

        # compute average observation for each unique grid-point
        tmpy[i] = n.average(y[idx])

    # walk through all grid points and find how many obs available
    #for i, iL in enumerate(Lgrid):
    #   idx = n.where( iL == Lobs[:])[0]

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
    """

    import numpy as n
    from matplotlib.mlab import find

    Tstart = Tnow - dd['model']['Twindow']
    Tend = Tnow

    # read all PSD observation and append into y,L
    y = []; L = []
    for satkey in dd['data']['satids']: # for every s/c
        PSD = dd['data'][satkey]['PSD']
        TAI = dd['data'][satkey]['TAI']
        Lvals = dd['data'][satkey]['L']
        indices = n.where( (Tstart < TAI) & (TAI <= Tend) )
        y = n.append(y, PSD[indices])
        L = n.append(L, Lvals[indices])
            
    # sort them and average if they fall into same L range
    ind = L.argsort() # return indices to a sorted y
    y = y[ind]
    L = L[ind]

    #check if observation available
    if len(y) == 0: return dd, n.array([]), n.array([])

    i = 0
    Lmean = []; ymean = []
    Lgrid = dd['model']['Lgrid']
    for iL in Lgrid[1:-1]:
        i += 1
        L_lower = n.mean([ Lgrid[i], Lgrid[i-1] ])
        L_upper = n.mean([ Lgrid[i+1], Lgrid[i] ])
        indices = n.where( (L_lower < L) & (L <= L_upper) )
        if n.size(indices) is 0: continue
        
        Lmean = n.append(Lmean, iL)
        ymean = n.append(ymean, n.mean(y[indices]))

    # add Lmean, ymean to dd
    try: # first time executing this
        dd['data']['dwindow'].append( (Tnow, Lmean, ymean) )
    except:
        dd['data']['dwindow'] = []
        dd['data']['dwindow'].append( (Tnow, Lmean, ymean) )

    # these line are only for averaging over L* as well!
    
    Lmean = n.mean(Lmean)
    # move to actual grid point
    print Lmean
    idx = find(Lgrid>Lmean)[0]
    
    Lmean = Lgrid[idx]
    ymean = n.mean(ymean)     
    return dd, n.array([Lmean]), n.array([ymean])

                     
# -----------------------------------------------
def output(init,result):
    """
    write results to file and be done
    """

    import sys


    return


# -----------------------------------------------
# ----- older stuff (not used) ------------------
# -----------------------------------------------
def assimilate_JK(dd):
    """
    this version is currently not working
    main function to assimilate all data provided in init
    """
    
    import sys
    import numpy as n

    n.random.seed(123)
    
    Lgrid = dd['model']['Lgrid']
    dd['model']['initPSD'] = n.array( 4e-7*n.exp( -(Lgrid-5)**2/2) )
    
    nens = dd['kalman']['nens']
    NL = dd['model']['ngrid']
    f0 = n.array(dd['model']['initPSD'])
    Tgrid = dd['model']['Tgrid']

    dd['results'] = {}
    dd['results']['fa'] = {}
    dd['results']['ff'] = {}
    
    # --- initialize some stuff
    # broadcast all f0 into A and randomize
    A = n.ones( (NL,nens) )*f0[:,n.newaxis]
    #A = 10**n.random.normal(n.log10(A),0.3)  # watch out for negative values

    
    dd['results']['fa']['Tgrid'] = n.zeros(len(Tgrid))
    dd['results']['ff']['Tgrid'] = n.zeros(len(Tgrid))
    dd['results']['fa']['PSD'] = n.zeros( (f0.shape[0],Tgrid.shape[0]) )
    dd['results']['fa']['PSD'][:,0] = f0
    dd['results']['fa']['Tgrid'][0] = Tgrid[0]
    dd['results']['ff']['PSD'] = n.zeros( (f0.shape[0],Tgrid.shape[0]) )
    dd['results']['ff']['PSD'][:,0] = f0    
    dd['results']['ff']['Tgrid'][0] = Tgrid[1]

    minPSD = n.min(A)/100

    # time loop to assimilate everything
    for Tnow, tcnt in zip(Tgrid, range(len(Tgrid)-1)):

        # make forecast using all ensembles in A
        icnt = 0
        for f in A.T:
            fnew = forecast(dd, f, Tnow)
            A[:,icnt] = fnew[:]
            icnt += 1

        #print n.log10(n.mean(A,axis=1))
        Tnow = Tnow + dd['model']['Twindow']
        tcnt = tcnt + 1
        dd['results']['ff']['PSD'][:,tcnt] = n.mean(A, axis=1)
        dd['results']['ff']['Tgrid'][tcnt] = Tnow

        # get observations for time window ]Tnow-Twindow,Tnow]        
        dd, L,y = getobs4window(dd, Tnow)

        # check if len(y) ==0
        if len(y) ==0:
            dd['results']['fa']['PSD'][:,tcnt] = n.mean(A, axis=1)
            dd['results']['fa']['Tgrid'][tcnt] = Tnow
            continue
            
        # perturb observations so that ensemble get sufficient spread
        HA = getHA(dd,L,A) # project ensemble states to obs. grid: HA(nobs,nens)
        e = n.zeros( (len(y),nens) )  # this is Psi in Evensen 2003
        D = n.zeros( (len(y),nens) )
        err_z = 0.3
        for yval, iobs in zip(y, range(len(y))):
            relstd = yval*err_z
            rnd = n.random.normal( yval*n.ones(nens), relstd)
            rnd[n.where(rnd<minPSD)] = minPSD
            e[iobs,:] = rnd - yval
            D[iobs,:] = yval + e[iobs,:]
            
        # add model error
        relstd = n.zeros(nens)
        for yval, iobs in zip(y, range(len(y))):
            idx = n.where( feq(L[iobs],Lgrid) )
            relstd[:] = 0.5*( yval - HA[iobs,:] )
            A[idx,:] = n.random.normal( HA[iobs,:], n.abs(relstd))
            idx2 = n.where(A[idx,:] < minPSD); A[idx,idx2] = minPSD # add min flux here
            
        # get residual or innovation in the classical sense
        Res = D - HA # dim = (nobs,nens)

        # error of covariant matrix
        ffmean = n.mean(A, axis=1)
        V = n.zeros((NL,nens))
        for iens in range(nens):
            V[:,iens]  = A[:,iens] - ffmean
            
        Pfxx = n.dot(V,V.T)/float(nens)
        Pfyy = n.dot(e,e.T)/float(nens)

        # get Kalman denominator (H*Pfxx*HT + Pfyy)^-1 by inversion
        PH, HPH = getHPH(dd,L,Pfxx)
        #print PH
        #print HPH
        Rinv = n.linalg.inv(HPH+Pfyy)
        Adj = n.dot(n.dot(PH,Rinv), Res)

        # update all ensemble members
        A = A + Adj

        # check for negative fluxes
        idx = n.where(A<minPSD)
        A[idx] = minPSD

        #
        # average A from analysis step and save in results dictionary
        #print y, n.mean(HA,axis=1), n.mean(getHA(init,L,A),axis=1)
        dd['results']['fa']['PSD'][:,tcnt] = n.mean(A, axis=1)
        dd['results']['fa']['Tgrid'][tcnt] = Tnow

        # print message
        print 'Tnow: ', rbtools.TAInum2ISOdate(Tnow)
        
        #print n.log10(n.mean(A,axis=1))[30]

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

    #from numpy import where, random, log10, mean, reshape, abs, max, arange
    #from numpy import linspace, min
    import numpy as n

    Lgrid = dd['model']['Lgrid']
    dL = Lgrid[1]-Lgrid[0]
    nens = int(dd['kalman']['nens'])
    #print y, L

    radius = 1
    
    for Lcenter, yval in zip(L, y):
        L1 = n.max((Lcenter-radius, Lgrid[0]))
        L2 = n.min((Lcenter+radius, Lgrid[-1]))
        #print Lcenter, Lcenter+radius, Lgrid[-1], min((Lcenter+radius,Lgrid[-1])), (L2-L1)/dL 
        NLs = int(round( (L2-L1)/dL ) + 1 )
        for Lpos in n.linspace(L1, L2, NLs):
            #print dL, L1, L2, NLs
            index = n.where( feq(Lpos,Lgrid) ) # use float point comparison
            stdev = 0.5*n.abs( yval - n.mean(A[index,:]) )
            #print Lpos
            center = n.reshape( A[index,:], (nens) )
            A[index,:] = n.random.normal( center, scale=stdev)
            # now check if any below 1e-99 and repeat
            for i in range(nens):
                icnt = 0
                while A[index,i] < 1e-99:
                    #A[index,i] = n.random.normal( center[i], scale=stdev)
                    A[index,i] = 10**n.random.normal( n.log10(center[i]), scale=0.001)
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
    



