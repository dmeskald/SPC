############################################################
############################################################
############## Accoradance and Discordance #################
############################################################
############################################################

compute_con=function(X, u, l){   #create function, X is the design matrix, u and l are the upper and lower thresholds
  N=nrow(X)  #N: regions (in X, rows represent the regions)
  T=ncol(X)  #T:time (in X,columns represent the time)
  
  normalize=function(X){  #create the normalization function
    X = X - mean(X)
    X = X/sd(X)
    return(X)
  }
  
  for(i in 1:N){  #call normalization function for each line of X matrix
    X[i,]=normalize(X[i,])
  }
  ######## X matrix is now normalized ###########
  
  Xu=array(0, dim=dim(X))  # create new Xu,Xl matrices with zeros and dimensions same as X
  Xl=array(0, dim=dim(X))  # same
  
  Xu[X>=u]=1  # apply the thresholds(u) on X and put ones(=1) in Xu for the corresponding positions
  Xl[X<=l]=-1 # apply the thresholds(l) on X and put ones(=1) in Xu for the corresponding positions
  
  Scav=1/T*(Xu%*%t(Xu))
  Scdv=1/T*(Xl%*%t(Xl))
  Sa=Scav+Scdv
  E=diag(Sa)^(-0.5)           # Calculate the energy
  Acc=diag(E)%*%Sa%*%diag(E)  # Calculate the Accordance matrix
  
  Sd=1/T*(Xu%*%t(Xl)+Xl%*%t(Xu))
  Dis=diag(E)%*%Sd%*%diag(E)  # Calculate the Accordance matrix
  
  con=NULL
  
  con$acc=Acc
  con$dis=Dis
  return(con) # return the con (type : con$acc for Acc matrix and con$Dis for Dis matrix)
}