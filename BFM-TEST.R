###########################################################
#
# FUNCTION
# KERN.k computes for x1 and x2 a value of the kernel k
# with a given value of the tuning parameter
#
# x1 and x2 are d-dimensional vectors
# two tuning parameters at most: theta1 and theta2 
#
# OUTPUT -> scalar
# theta should be equal to 2*sigma^2 (see Gaussian RBF Kernel, 
# https://en.wikipedia.org/wiki/Radial_basis_function_kernel)
# JD uses in his Dresden-presentation \gamma^2 instead of \theta
KERN.k = function(x1,x2, theta1, theta2=NULL)
{
  # vectorwise operations
  out= exp(-sum((x1-x2)^2)/theta1)
  return(out)
}

#########################################################
#
# FUNCTION
# STAND.h computes for z1 and z2 a value of the kernel h
#
# See Equation (3) of submission and the subsequent formula, Page 4 
# DEPENDS on KERN.k
# 
# d  is dimension  
# z1 is a vector of length 2*d
# z2 is a vector of length 2*d
# theta1 and theta2=NULL -> tuning parameters of KERN.k
# output -> scalar

STAND.h = function(d, z1, z2, theta1, theta2=NULL)
{
  x1 = z1[1:d]
  y1 = z1[(d+1):(2*d)]
  x2 = z2[1:d]
  y2 = z2[(d+1):(2*d)]
  
  t1 = KERN.k(x1,x2, theta1=theta1, theta2=NULL)
  t2 = KERN.k(y1,y2, theta1=theta1, theta2=NULL)
  t3 = KERN.k(x1,y2, theta1=theta1, theta2=NULL)
  t4 = KERN.k(x2,y1, theta1=theta1, theta2=NULL)
  
  h  = t1+t2-t3-t4 
  return(h)
}


#########################################################
#
# FUNCTION
# MMD.stand, see Equation (3) in the Submission
#
# d is dimension of X-data and Y-data
# n is sample size
# X is a matrix of dimension nxd. Rows are iid
# If d=1 then transform X with as.matrix
# Y is a matrix of dimension nxd
# If d=1 then transform Y with as.matrix
#
MMD.stand = function(d, n,  X, Y, theta1, theta2=NULL)
{
  ustat = 0
  
  # Create Z.i vector consisting of X.i and Y.i
  # Create Z with 'cbind'
  Z = cbind(X, Y)
  
  # Two  loops for U-stat
  for(i in 1:(n-1))
  {
    for(j in (i+1):n)
    {
      ustat = ustat + STAND.h(d, Z[i, ], Z[j,], theta1, theta2=NULL)
    }
  }
  
  ustat = 2*ustat/(n*(n-1))
  return(ustat)
}  


#####################################################################
# FUNCTION
# NEW.h 
# Equation (4) in the submission, Page 8
# Computes for z1 and z2 a value of the kernel h
# DEPENDS on KERN.k
# 
# d  is dimension  
# z1 is a vector of length 4*d, 
#       (X_{2i-1}, Y_{2i-1}, X_{2i}, Y_{2i})
# z2 is a vector of length 4*d, #
#       (X_{2j-1}, Y_{2j-1}, X_{2j}, Y_{2j}) 
# theta1 and theta2=NULL -> tuning parameters of KERN.k
# output -> scalar

NEW.h = function(d, z1, z2, theta1, theta2=NULL)
{
  # First argument of the kernel q
  # See Equation (5) in pdf
  x1 = z1[1:d]
  y1 = z1[(d+1):(2*d)]
  x2 = z1[(2*d+1):(3*d)]
  y2 = z1[(3*d+1):(4*d)]
  
  # Second argument of the kernel q
  # See Equation (5) in pdf
  x3 = z2[1:d]
  y3 = z2[(d+1):(2*d)]
  x4 = z2[(2*d+1):(3*d)]
  y4 = z2[(3*d+1):(4*d)]
  
  # Terms with "+"
  t1 = KERN.k(x1,x3, theta1=theta1, theta2=NULL)
  t2 = KERN.k(y1,y3, theta1=theta1, theta2=NULL)
  # Terms with "-"
  t3 = KERN.k(x2,y4, theta1=theta1, theta2=NULL)
  t4 = KERN.k(x4,y2, theta1=theta1, theta2=NULL)
  
  h  = t1+t2-t3-t4 
  return(h)
}

################################################################
# FUNCTION
# U-statistics
# This is MMD_q, Equation (4)
#
# DEPENDS on KERN.k
# DEPENDS on NEW.h
#
# d is dimension of X-data and Y-data
# n is sample size, it should be even 
# X is a matrix of dimension nxd. 
# If d=1 then transform X with as.matrix
# Y is a matrix of dimension nxd
# If d=1 then transform Y with as.matrix
#
U.stat = function(d, n,  X, Y, theta1, theta2=NULL)
{
  ustat = 0
  #
  nn = n/2
  
  # Create Z.i vector consisting of X[2*i-1, ], Y[2*i-1, ], X[2*i, ], Y[2*i, ]
  
  # Separate 'even'and 'odd' rows.
  IND <- 2*(1:nn) - 1
  
  # Create Z with 'cbind'
  Z = cbind(X[IND,], Y[IND,], X[IND+1,], Y[IND+1,])
  
  
  for(i in 1:(nn-1))
  {
    for(j in (i+1):nn)
    {
      # two observations is one argument
      #Z.i =  c(X[2*i-1, ], Y[2*i-1, ], X[2*i, ], Y[2*i, ])
      #Z.j =  c(X[2*j-1, ], Y[2*j-1, ], X[2*j, ], Y[2*j, ])
      
      ustat = ustat + NEW.h(d, Z[i, ], Z[j,], theta1, theta2=NULL)
    }
  }
  
  ustat = 2*ustat/(nn*(nn-1))
  return(ustat)
}  

########################################################

# Simulate data


SIM.mu.fig.1 = function(d, n, SD.temp=1.2)
{
  # matrices of dimension nxd
  X = matrix(rnorm(n*d, mean=0), ncol = d)
  # SD of Y-data is equal to 3
  Y = matrix(rnorm(n*d, mean=0, sd=SD.temp), ncol = d)
  
  # Estimate mean using X-data
  hat.mu.X = apply(X, 2, mean)
  
  # Matrix, each row is hat.mu
  matr.hat.mu =  matrix(hat.mu.X, ncol=d, nrow=n, byrow = T)
  # Y + \alpha
  hat.Y  = Y + matr.hat.mu
  
  out = list("X"=X, "Y"=Y, "hat.Y"=hat.Y, "hat.mu.X"=hat.mu.X )
  return(out)
}

##
##  Original data X from the multiavraite normal with expectation zero and 
##  identity matrix as the covariance
##  SD.theor1=c(1,1) is the standard deviation of the first and the second half
##  of  Sample 1
##  SD.theor1=c(1,1) is the standard deviation of the first and the second half
##  of Sample 2


SIM.mu.2sampl = function(d, n, SD.theor1=c(1,1), SD.theor2=c(1,1))
{
  
  # matrices of dimension nxd
  # Theoretical SD of X is always 1
  X = matrix(rnorm(n*d, mean=0), ncol = d)
  
  dd = d/2
  
  ######
  # First sample
  # SD of Y-data should be defined
  Y1 = matrix(rnorm(n*dd, mean=0, sd=SD.theor1[1]), ncol = dd)
  Y2 = matrix(rnorm(n*dd, mean=0, sd=SD.theor1[2]), ncol = dd)
  
  Y = cbind(Y1, Y2)
  
  # Estimate mean using X-data
  hat.mu.X = apply(X, 2, mean)
  
  # Matrix, each row is hat.mu
  matr.hat.mu =  matrix(hat.mu.X, ncol=d, nrow=n, byrow = T)
  # Y + \alpha
  hat.Y1  = Y + matr.hat.mu
  
  ########
  # Second sample
  # SD of Y-data should be defined
  # SD are switched
  Y1 = matrix(rnorm(n*dd, mean=0, sd=SD.theor2[1]), ncol = dd)
  Y2 = matrix(rnorm(n*dd, mean=0, sd=SD.theor2[2]), ncol = dd)
  
  Y = cbind(Y1, Y2)
  
  # Matrix, each row is hat.mu and  matr.hat.mu are computed above 
  hat.Y2  = Y + matr.hat.mu
  
  out = list("X"=X, "Y"=Y, "hat.Y1"=hat.Y1, 
             "hat.Y2"=hat.Y2, "hat.mu.X"=hat.mu.X )
  return(out)
}


###############################################################
### This function should be used only internally in
# BFM.test function since it requires "MMD.hat"

SIGMA.spec.a =function(d, n,  X, Y, theta1, theta2=NULL, MMD.hat)
{
  #############################
  # Create Z.i vector consisting of X.i and Y.i
  # Create Z with 'cbind'
  Z = cbind(X, Y)
  
  ## Compute \sigma_\alpha
  h1 <- array(0, c(n))
  for(i in 1:n)
  {
    for(j in 1:n)
    {
      h1[i] = h1[i] + STAND.h(d, Z[i, ], Z[j,], theta1, theta2=NULL) 
    }
    h1[i] = h1[i]/n - MMD.hat
  }
  
  s.a.squared   = 4*var(h1)
  s.a           = sqrt(s.a.squared)
  return(list("s.a"=s.a, "h1.a"=h1))
}

############################################################
### This function should be used only internally in
# BFM.test function since it requires "MMD.hat"

SIGMA.spec.q =function(d, n,  X, Y, theta1, theta2=NULL, MMD.hat)
{
#############################
## Compute \sigma_q
nn = n/2

# Create Z.i vector consisting of X[2*i-1, ], Y[2*i-1, ], X[2*i, ], Y[2*i, ]

# Separate 'even'and 'odd' rows.
IND <- 2*(1:nn) - 1

# Create Z with 'cbind'
Z = cbind(X[IND,], Y[IND,], X[IND+1,], Y[IND+1,])

h1 <- array(0, c(nn))
for(i in 1:nn)
{
  for(j in 1:nn)
  {
    h1[i] = h1[i] + NEW.h(d, Z[i, ], Z[j,], theta1, theta2=NULL) 
  }
  h1[i] = h1[i]/nn - MMD.hat
}

s.q.squared   = 2*4*var(h1)

s.q = sqrt(s.q.squared)

out = list("s.q"=s.q,  "h1.q"=h1)
return(out)
}

##########################################

# depends on KERN.k, STAND.h, MMD.stand, NEW.h, U.stat, 
# d = dimension
# n = sample size, should be n=2*k, where k is a natural number
#
# Y = Y1 is the first sample
# Y2 is the second sample
# epsilon can be also a vector
# spec.reject.1 = 1 means reject

##### 
BFM.test = function(d, n, epsilon,  X, Y, 
                    Y2=NULL,  theta1, theta2=NULL, 
                    model.selection=FALSE, level=0.05)
{
  
  ###########################
  # model.selection=FALSE
  ###########################
  if(!model.selection)
  {
    # Compute "classical" MMD
    # Equation (3) on Page 4 of the submission
    MMD.hat = MMD.stand(d, n,  X, Y, theta1, theta2=NULL)
    
    ###########################
    # Compute MMD.eps
    # Equation prior to Equation (5) on Page 8
    
    # MMD_q, Equation (4) on Page 8 of the submission
    MMD.q = U.stat(d=d, n=n,  X=X, Y=Y, theta1, theta2=NULL)
    
    # MMD_eps, Equation prior to Equation (5) on Page 8
    MMD.eps = MMD.hat + epsilon*MMD.q
    
    #############################
    s.a = SIGMA.spec.a(d, n,  X, Y, theta1, theta2=NULL, MMD.hat=MMD.hat)$s.a
    s.q = SIGMA.spec.q(d, n,  X, Y, theta1, theta2=NULL, MMD.hat=MMD.hat)$s.q
    
    # this is a vector, Standard error from Theorem 3.1
    sd.mmd.eps = s.a + epsilon*s.q 
    
    #############################
    # this is a vector, Test statistics from Theorem 3.1
    spec.test.stat.1 = sqrt(n)*MMD.eps/sd.mmd.eps
    spec.reject.1    = (abs(spec.test.stat.1) >  qnorm(1-level/2))
    
    spec.test.stat.q.1 =  sqrt(n)*(MMD.q)/s.q
    spec.reject.q.1    =  ( abs(spec.test.stat.q.1) > qnorm(1-level/2) )
    
    #############################
    # Output 
    out = list("selec.test.stat" = NULL,
               "selec.reject" = NULL,
               "selec.test.stat.q"=NULL,
               "selec.reject.q" = NULL,
               "sd.mmd.sel"=NULL,
               "sd.mmd.sel.q"= NULL,
               #
               "spec.test.stat.1" = spec.test.stat.1, 
               "spec.reject.1" = spec.reject.1,
               "MMD.eps.1" = MMD.eps,
               "MMD.hat.1" = MMD.hat, 
               "MMD.q.1" = MMD.q,  
               "spec.test.stat.q.1" = spec.test.stat.q.1, 
               "spec.reject.q.1" = spec.reject.q.1,
               "s.a.1" = s.a, "s.q.1" = s.q,
               #
               "spec.test.stat.2" = NULL,
               "spec.reject.2" = NULL,
               "MMD.eps.2" = NULL,
               "MMD.hat.2" = NULL, 
               "MMD.q.2" = NULL,  
               "spec.test.stat.q.2" = NULL,
               "spec.reject.q.2" = NULL,
               "s.a.2" = NULL, "s.q.2" = NULL)
    
  }
  
  ###########################
  # model.selection=TRUE
  ###########################
  
  if(model.selection)
  {
    ### SAMPLE  1, i.e. Y
    # Compute "classical" MMD
    # Equation (3) on Page 4 of the submission
    MMD.hat.1 = MMD.stand(d, n,  X, Y, theta1, theta2=NULL)
    
    ###########################
    # Compute MMD.eps
    # Equation prior to Equation (5) on Page 8
    
    # MMD_q, Equation (4) on Page 8 of the submission
    MMD.q.1 = U.stat(d=d, n=n,  X=X, Y=Y, theta1, theta2=NULL)
    
    # Equation prior to Equation (5) on Page 8
    MMD.eps.1 = MMD.hat.1 + epsilon*MMD.q.1
    
    out.spec.a.1 = SIGMA.spec.a(d, n,  X, Y, theta1, theta2=NULL, MMD.hat=MMD.hat.1) 
    s.a.1        =  out.spec.a.1$s.a
    
    out.spec.q.1 = SIGMA.spec.q(d, n,  X, Y, theta1, theta2=NULL, MMD.hat=MMD.hat.1) 
    s.q.1        =  out.spec.q.1$s.q
    
    #############################
    # this is a vector, Standard error from Theorem 3.1
    sd.mmd.eps.1 = s.a.1 + epsilon*s.q.1 
    
    #############################
    # this is a vector, Test statistics from Theorem 3.1
    spec.test.stat.1 = sqrt(n)*MMD.eps.1/sd.mmd.eps.1
    spec.reject.1    = ( abs(spec.test.stat.1) > qnorm(1-level/2) )
    
    spec.test.stat.q.1 = sqrt(n)*(MMD.q.1)/s.q.1
    spec.reject.q.1    =  ( abs(spec.test.stat.q.1) > qnorm(1-level/2) )
    
    #############################################################
    #############################################################
    ### SAMPLE  2, i.e. Y2
    # Compute "classical" MMD
    # Equation (3) on Page 4 of the submission
    MMD.hat.2 = MMD.stand(d, n,  X, Y2, theta1, theta2=NULL)
    
    ###########################
    # Compute MMD.eps
    # Equation prior to Equation (5) on Page 8
    
    # MMD_q, Equation (4) on Page 8 of the submission
    MMD.q.2 = U.stat(d=d, n=n,  X=X, Y=Y2, theta1, theta2=NULL)
    
    # Equation prior to Equation (5) on Page 8
    MMD.eps.2 = MMD.hat.2  + epsilon*MMD.q.2
  
    out.spec.a.2 = SIGMA.spec.a(d, n,  X, Y2, theta1, theta2=NULL, MMD.hat=MMD.hat.2)
    s.a.2        =  out.spec.a.2$s.a
    
    out.spec.q.2 = SIGMA.spec.q(d, n,  X, Y2, theta1, theta2=NULL, MMD.hat=MMD.hat.2)
    s.q.2        =  out.spec.q.2$s.q
    
    #############################
    # this is a vector, Standard error from Theorem 3.1
    sd.mmd.eps.2 = s.a.2 + epsilon*s.q.2 
    
    #############################
    # this is a vector, Test statistics from Theorem 3.1
    spec.test.stat.2 = sqrt(n)*MMD.eps.2/sd.mmd.eps.2
    spec.reject.2    = ( abs(spec.test.stat.2) > qnorm(1-level/2) )
    
    spec.test.stat.q.2 = sqrt(n)*(MMD.q.2)/s.q.2
    spec.reject.q.2 =  ( abs(spec.test.stat.q.2) > qnorm(1-level/2) )
    
    #################################################################
    ##########  MODEL SELECTION
    # Estimate varainces 
    
    h1.a.1 = out.spec.a.1$h1.a
    h1.a.2 = out.spec.a.2$h1.a
    
    s.a.squared   = 4*var(h1.a.1-h1.a.2)
    
    h1.q.1 = out.spec.q.1$h1.q
    h1.q.2 = out.spec.q.2$h1.q
    
    s.q.squared   = 2*4*var(h1.q.1-h1.q.2)
    # Here we ignore the covariance
  
    sd.mmd.sel = sqrt(s.a.squared) + epsilon*sqrt(s.q.squared)
    
    selec.test.stat = sqrt(n)*(MMD.eps.1-MMD.eps.2)/sd.mmd.sel
    selec.reject    = ( abs(selec.test.stat) > qnorm(1-level/2) )
    
    selec.test.stat.q = sqrt(n)*(MMD.q.1-MMD.q.2)/sqrt(s.q.squared)
    selec.reject.q =  ( abs(selec.test.stat.q) > qnorm(1-level/2) )
    
    #############################
    # Output 
    out = list("selec.test.stat"=selec.test.stat,
               "selec.reject"= selec.reject,
               "selec.test.stat.q"=selec.test.stat.q,
               "selec.reject.q" = selec.reject.q,
               "sd.mmd.sel"=sd.mmd.sel,
               "sd.mmd.sel.q"= sqrt(s.q.squared),
               #
               "spec.test.stat.1"=spec.test.stat.1,
               "spec.reject.1"= spec.reject.1,
               "MMD.eps.1"=MMD.eps.1,
               "MMD.hat.1"=MMD.hat.1, 
               "MMD.q.1"= MMD.q.1,  
               "spec.test.stat.q.1" = spec.test.stat.q.1, 
               "spec.reject.q.1" = spec.reject.q.1,
               "s.a.1"=s.a.1, "s.q.1"=s.q.1,
               #
               "spec.test.stat.2"=spec.test.stat.2,
               "spec.reject.2"= spec.reject.2,
               "MMD.eps.2"=MMD.eps.2,
               "MMD.hat.2"=MMD.hat.2, 
               "MMD.q.2"= MMD.q.2,  
               "spec.test.stat.q.2" = spec.test.stat.q.2, 
               "spec.reject.q.2" = spec.reject.q.2,
               "s.a.2"=s.a.2, "s.q.2"=s.q.2)
    
  }
  
  return(out)
}



n      = 250
d      = 16
theta1 = d

SD.temp = 1.2

pow     = -1/c(2.5, 4.5, 6.5)
epsilon = c( n^(pow))

# Model specification
dat = SIM.mu.fig.1(d, n, SD.temp)
    
X     = dat$X
hat.Y = dat$hat.Y
    
res = BFM.test(d=d, n=n, epsilon= epsilon,  
             X=X, Y=hat.Y, Y2=NULL,  
             theta1, theta2=NULL, 
             model.selection=FALSE, level=0.05)
    
res    


## Model comparison
dat.mod.comp=SIM.mu.2sampl(d, n, SD.theor1=c(1.2,1.3), SD.theor2=c(1.2,1.3))


X     = dat.mod.comp$X
hat.Y1 = dat.mod.comp$hat.Y1
hat.Y2 = dat.mod.comp$hat.Y2

res.mod.comp = BFM.test(d=d, n=n, epsilon= epsilon,  
               X=X, Y=hat.Y1, Y2=hat.Y2,  
               theta1, theta2=NULL, 
               model.selection=TRUE, level=0.05)

res.mod.comp  

 

