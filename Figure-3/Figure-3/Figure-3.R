

# First change the working directory in R

ST.ERROR <- round(seq(1, 3, by=0.1), digits=1)


SAMPLE.SIZE = c(100, 250, 500,1000)

DIM = c(2^c(1,4))

pow = -1/c(2.5, 4.5, 6.5)

l.epsilon = length(pow) 

res    <- array(0, c(length(DIM), length(SAMPLE.SIZE),  l.epsilon, 5))
res.SW <- array(0, c(length(DIM), length(SAMPLE.SIZE), 5) )


for(i.sd.new in 1:5 )
{
  

  load(paste("Worspace-Fig-3-SD-of-Y-is-", ST.ERROR[i.sd.new],".RData", sep=""))
  
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
  pdf(paste("2-models-DEGENER-Fig3-dim-", DIM[iii], ".pdf", sep =""))
  plot(ST.ERROR[1:5], res.SW[iii, 1,  ],  pch=20, type="b",  
       lwd=2, col=1,  main= paste("p=", DIM[iii], sep=""), 
       ylim=c(0,1), xlim=c(1,1.7),       xlab="Stand. dev.", ylab="Emp. Level/Power")
  for(i in 1: length(SAMPLE.SIZE))
  {
    if( i!=1 ) 
    {
      points(ST.ERROR[1:5], res.SW[iii, i,  ], pch=20, col=i)
      lines(ST.ERROR[1:5], res.SW[iii, i, ], col=i, lwd=2, lty=i)
      
    }
    
    
    # only one epsilon
    # for(j in c(1,3,5))
    for(j in c(1))
    {
      points(ST.ERROR[1:5], res[iii, i, j, ], pch=20+j, col=i)
      lines(ST.ERROR[1:5], res[iii, i, j, ], col=i, lwd=2,lty=i)
    }

    leg.txt <- c("MMD_q", "eps=n^(-1/2.5)") 
    legend(1.5, 1, leg.txt, pch=c(20, 21), col=rep(5,2))
    
    leg.txt <- c("n=100","n=250", "n=500", "n=1000") 
    legend(1.5, 0.5, leg.txt, lty=c(1,2,3,4), col=c(1,2,3,4),  lwd=c(2,2,2,2))
    
    
  }
  abline(h=0.05, lty=2)
  abline(h=0, lty=1, col=1)
  
  dev.off()
}

