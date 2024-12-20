library(psych)
cat("\014")
rm(list=ls())

set.seed(100)


# Number of Kernels
NK = 20

# NR = 100 
# NR = 2
NR = 50

ST.ERROR <- seq(1.0, 1.4, by=0.1)

for (i.sd in 1:length(ST.ERROR))
{
  set.seed(100)
  
  #Global variable SD.temp
  SD.temp <<- ST.ERROR[i.sd]

  library(foreach)
  library(doParallel)
  library(parallel)
  library(iterators)
  library(snow)


  source("BFM-TEST.R")

  
  SAMPLE.SIZE = c( 100, 250,  500, 1000)

  #DIM = c(2^(1:8))
  DIM = c(2^c(1,4))

  THETA1 = DIM

  SEED <- sample(1:10^6, NK)

  pow = -1/c(2.5, 4.5, 6.5)

  TEST.STAT   <- array(0, c(length(DIM), length(SAMPLE.SIZE),  NR*NK, length(pow) ) )
  TEST.STAT.q <- array(0, c(length(DIM), length(SAMPLE.SIZE),  NR*NK  ) )
  
  REJECT     <- array(0, c(length(DIM), length(SAMPLE.SIZE),  NR*NK, length(pow) ) )
  REJECT.q   <- array(0, c(length(DIM), length(SAMPLE.SIZE),  NR*NK  ) )
  

  for(iii in 1:length(DIM))
  {  
  print("DIM")
  print(iii)
  
    d      = DIM[iii]
    theta1 = THETA1[iii]
  
    for(ijk in 1:length(SAMPLE.SIZE))
    {
      n = SAMPLE.SIZE[ijk] 
      print("n")
      print(ijk)
      
    
      epsilon = c( n^(pow))
    
    
      cl <- makeCluster(NK, "SOCK") #Use 20 Cores
      registerDoParallel(cl)
    
      #Load all needed libraries into cluster
      clusterEvalQ(cl, library(psych))
      clusterEvalQ(cl, library(foreach))
      clusterEvalQ(cl, library(doParallel))
      clusterEvalQ(cl, library(parallel))
      clusterEvalQ(cl, library(iterators))
      clusterEvalQ(cl, library(snow))
    
      # Export R-functions
      clusterExport(cl, "Ex.Funct.Fig.3")
      clusterExport(cl, "fun.parallel.Fig.3")
    
      clusterExport(cl, "KERN.k")
      clusterExport(cl, "STAND.h")
      clusterExport(cl, "MMD.stand")
      clusterExport(cl, "NEW.h")
      clusterExport(cl, "U.stat")
      clusterExport(cl, "BFM.test")
      
      clusterExport(cl, "SIM.mu.2sampl")
    
      clusterExport(cl, "SIGMA.spec.a")
      clusterExport(cl, "SIGMA.spec.q")
    
      # Import data
      clusterExport(cl, "NR")
      clusterExport(cl, "d")
      clusterExport(cl, "n")
      clusterExport(cl, "SD.temp")
      clusterExport(cl, "epsilon")
      clusterExport(cl, "SEED")
      clusterExport(cl, "theta1")
    
      # Parallelization!
      outT <- foreach(x = 1:NK) %dopar% {fun.parallel.Fig.3(x)}
    
      stat.eps   <- c()
      reject     <- c() 
      stat.q     <- c()
      reject.q   <- c()
      
      
      #combine output from each cluster
      for(i in 1:NK)
      {
      
        stat.eps   <- rbind(stat.eps, outT[[i]]$selec.test.stat)
        reject     <- rbind(reject, outT[[i]]$selec.reject)
        
        
        stat.q   <- c(stat.q, outT[[i]]$selec.test.stat.q)
        reject.q <- c(reject.q, outT[[i]]$selec.reject.q)
         
      }
    
      stopCluster(cl)
    
      TEST.STAT[iii,ijk, , ]  = stat.eps
      TEST.STAT.q[iii,ijk, ]  = stat.q
      
      REJECT[iii,ijk, , ]  = reject
      REJECT.q[iii,ijk, ]  = reject.q 
      
      print("SD")
      print(i.sd)

      print( c("SD", "DIM", "n") )
      print( c(i.sd, iii, ijk) )
    }
    
  }

save.image(paste("Worspace-Fig-3-SD-of-Y-is-", SD.temp, ".RData",sep="" ))

}



