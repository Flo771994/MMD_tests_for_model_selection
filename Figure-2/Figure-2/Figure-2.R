

MU.mean <- seq(0, 0.6, by=0.1)


SAMPLE.SIZE = c(100, 250, 500,1000)

DIM = c(2^c(1,4))

pow = -1/c(2.5, 4.5, 6.5)

l.epsilon = length(pow) 

res    <- array(0, c(length(DIM), length(SAMPLE.SIZE),  l.epsilon, length(MU.mean)))
res.SW <- array(0, c(length(DIM), length(SAMPLE.SIZE), length(MU.mean)) )


for(i.sd.new in 1:length(MU.mean) )
{
  
  load(paste("Worspace-MU-of-Y-is-", MU.mean[i.sd.new],".RData", sep=""))
  
  PVALUE    <- array(-1, c(length(DIM), length(SAMPLE.SIZE), length(epsilon)) ) 
  PVALUE.SW <- array(-1, c(length(DIM), length(SAMPLE.SIZE) ) )

  for(iii in 1:length(DIM))
  {
    for(i in 1: length(SAMPLE.SIZE))
    {
        PVALUE.SW[iii, i] <- mean(REJECT.q[iii,i, ]) 
        
          for(j in 1:length(epsilon) )
          {
            PVALUE[iii, i, j] <- mean( REJECT[iii,i, , j]) 
          }
    }
  }

  res[ , , , i.sd.new]  = PVALUE

  res.SW[ , , i.sd.new] = PVALUE.SW
}


dim(res.SW)
dim(res)


DIM = c(2^c(1, 4))


for (iii in 1:length(DIM))
{
  pdf(paste("x-mean-pvalue-dim-", DIM[iii], ".pdf", sep =""))
  plot(MU.mean[1:length(MU.mean)], res.SW[iii, 1,  ],  pch=20, type="b",  
       lwd=2, col=1,  main= paste("p=", DIM[iii], sep=""), 
       ylim=c(0,1), xlim=c(0,1.1),       xlab="Mean", ylab="Emp. Level/Power")
  for(i in 1: length(SAMPLE.SIZE))
  {
    if( i!=1 ) 
    {
      points(MU.mean[1:length(MU.mean)], res.SW[iii, i,  ], pch=20, col=i)
      lines(MU.mean[1:length(MU.mean)], res.SW[iii, i, ], col=i, lwd=2, lty=i)
      
    }
    
    
    # only one epsilon
    # for(j in c(1,3,5))
    for(j in c(1))
    {
      points(MU.mean[1:length(MU.mean)], res[iii, i, j, ], pch=20+j, col=i)
      lines(MU.mean[1:length(MU.mean)], res[iii, i, j, ], col=i, lwd=2, lty=i)
    }

    leg.txt <- c("MMD_q", "eps=n^(-1/2.5)") 
    legend(0.8, 1, leg.txt, pch=c(20, 21), col=rep(5,2))
    
    leg.txt <- c("n=100", "n=250", "n=500", "n=1000") 
    legend(0.8, 0.5, leg.txt, lty=c(1, 2, 3, 4), col=c(1,2,3,4),  lwd=c(2,2,2,2))
    
    
  }
  abline(h=0.05, lty=2)
  abline(h=0, lty=1, col=1)
  
  dev.off()
}

