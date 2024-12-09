
# Set working directory correspondingly

ST.ERROR <- round(seq(1, 3, by=0.1), digits=1)

SAMPLE.SIZE =c(100, 250, 500, 1000)

DIM = c(2^c(1,4))

pow = -1/c(2.5, 4.5, 6.5)

l.epsilon = length(pow) 

res    <- array(0, c(length(DIM), length(SAMPLE.SIZE),  l.epsilon, 5))
res.SW <- array(0, c(length(DIM), length(SAMPLE.SIZE), 5) )


for(i.sd.new in 1:5)
{
  

  load(paste("Worspace-MP-0-SD-of-Y-is-", ST.ERROR[i.sd.new],".RData", sep=""))

    PVALUE    <- array(-1, c(length(DIM), length(SAMPLE.SIZE), length(epsilon)) ) 
  PVALUE.SW <- array(-1, c(length(DIM), length(SAMPLE.SIZE) ) )

  for(iii in 1:length(DIM))
  {
    for(i in 1: length(SAMPLE.SIZE))
    {
        PVALUE.SW[iii, i] <- mean(REJECT.q[iii,i, ]) # abs(SW.STAT[iii, i, ]) > qnorm(0.975) )
        
          for(j in 1:length(epsilon) )
          {
            PVALUE[iii, i, j] <- mean( REJECT[iii,i, , j]) #  abs(TEST.STAT[iii, i, , j]) > qnorm(0.975) 
          }
    }
  }

  res[ , , , i.sd.new]  = PVALUE

  res.SW[ , , i.sd.new] = PVALUE.SW
}


DIM = c(2^c(1, 4))


for (iii in 1:length(DIM))
{
  pdf(paste("pvalue-dim-", DIM[iii], "-black-lines.pdf", sep =""))
  
  # only sample size n=500, it is the third entry
  plot(ST.ERROR[1:5], res.SW[iii, 3,  1:5],  pch=20, type="b",  
       lwd=2, col=1, lty=1,  main= paste("p=", DIM[iii], sep=""), 
       ylim=c(0,1), xlim=c(1,1.7),       xlab="Stand. dev.", ylab="Emp. Level/Power")
  for(i in 3:3) # only n=500
  {
    for(j in 1:l.epsilon)
    {
      points(ST.ERROR[1:5], res[iii, i, j,1:5 ], pch=20 + 2*(j-1)+1, col=1)
      lines(ST.ERROR[1:5], res[iii, i, j, 1:5], col=1, lty=1,  lwd=2)
    }
    leg.txt <- c("MMD_q", "eps=n^(-1/2.5)", "eps=n^(-1/4.5)", "eps=n^(-1/6.5)") 
    legend(1.5, 1, leg.txt, pch=c(20, 21, 23, 25), col=rep(1,4))
    
  }
  abline(h=0.05, lty=2)
  abline(h=0, lty=1, col=1)
  
  dev.off()
}

