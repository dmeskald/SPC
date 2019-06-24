################     SERAFEIM LOUKAS     ##################
################ Accordance and Discordance ###############
###########################################################
# To call it type:"python conn.py" and you should get the
# Acc and Dis matrices
###########################################################

from numpy import * #import numpy module

def my_func(X,u,l): #define main function
  X=array(X) #type of X array
  N=X.shape[0] #N:regions (in X,the rows represent the regions)
  T=X.shape[1] #T:time (in X,the columns represent the time)

  def normalization(X): #define normalization function
    X= X - mean(X) 
    X= X / std(X)
    return X
  for i in range(N): #call normalization func for each line of X
    X[i:] = normalization(X[i:])
  # X is now normalized
  Xu = array([[0]*T]*N) # create zero arrays
  Xl = array([[0]*T]*N)
  
  Xu[X >= u] = 1 # apply thresholds on X and put 1 in Xu 
  Xl[X <= l] = -1 #for thecorresponding positions 
  
  Txu=transpose(Xu) #transpose Xu and Xl
  Txl=transpose(Xl) #we use them in the next lines
  
  Scav =(float(1) / float(T))*(dot(Xu,Txu))
  Scdv =(float(1) / float(T))*(dot(Xl,Txl))
  Sa = Scav + Scdv
  E =diag(Sa) ** (- 0.5) #Calculate energy
  d_E=diag(E)
  Acc =d_E.dot(Sa).dot(d_E) #Calculate Accordance array
  
  Sd =(float(1) /float(T))*(dot(Xu,Txl) + dot(Xl,Txu))
  Dis=d_E.dot(Sd).dot(d_E) #Calculate Discordance array

  return Acc, Dis

################################################################
# CALL MAIN FUNCTION USING X AS INPUT

random.seed(1) # control the random function-just get same randoms
# create array from normal distr with mean=0 and std=1
X = random.normal(0, 1, (10,20)) 

# call the main function for the defined X

acc, dis = my_func(X,1,1)

print "Accordance is \n",acc
print "and the Discordance is \n", dis
print "finished"
