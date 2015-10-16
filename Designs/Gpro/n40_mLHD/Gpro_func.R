##########################################


yeastODE=function(paramValuesMatrix,u)
{
  library(deSolve)
  place_holder=dim(paramValuesMatrix)
  response=mat.or.vec(1,place_holder[2])
  for(i in 1:place_holder[2])
  {
    k3=paramValuesMatrix[3,i]
    k5=paramValuesMatrix[5,i]
    Gt=paramValuesMatrix[place_holder[1],i]
    x=lsoda(y = c(k5/k3,0,Gt,0), times = c(0,60), func = LOCAL_ode, parms =c(paramValuesMatrix[,i],u[i]),atol = 1e-8, rtol = 1e-5)
    
    Gt=paramValuesMatrix[place_holder[1],i]
    place_holder1=dim(x)
    y=(Gt-x[place_holder1[1],4])/Gt
    
    response[1,i]=y
  }
  return(response)  
}


LOCAL_ode=function(t,x,k)
{
  u=k[10]
  x1d = -k[1]*x[1]*u + k[2]*x[2]-k[3]*x[1]+k[5]
  x2d = k[1]*x[1]*u-k[2]*x[2]-k[4]*x[2]
  x3d = -k[6]*x[2]*x[3]+k[8]*(k[9]-x[3]-x[4])*(k[9]-x[3])
  x4d = k[6]*x[2]*x[3]-k[7]*x[4]
  
  x_dot=list(c(x1d,x2d,x3d,x4d))
  
  return(x_dot)	
}

##########################
##########################

GProtein =function(x)
{xd=x
 place_holder=dim(xd)
 nr=place_holder[1]
 d=place_holder[2]
 
 u1=xd[,1]
 u6=xd[,2]
 u7=xd[,3]
 
 x=xd[,4]
 #Parmeters scaled to orginal values
 k1=2e6*(2e7/2e6)^u1
 k6=3e-5*(3e-4/3e-5)^u6
 k7=0.3*(3/0.3)^u7
 #Scaled x vairbale
 u=1e-9*(1e-6/1e-9)^x
 
 #parameters held constant.
 k2=5e-3*(mat.or.vec(nr,1)+1)
 k3=1e-3*(mat.or.vec(nr,1)+1)
 k4=2e-3*(mat.or.vec(nr,1)+1)
 k5=8*(mat.or.vec(nr,1)+1)
 k8=(mat.or.vec(nr,1)+1)
 Gt=10000*(mat.or.vec(nr,1)+1)
 #Construct the Parameter matrix which is #9XN (N=33)
 parameter=rbind(t(k1), t(k2), t(k3), t(k4), t(k5), t(k6), t(k7), t(k8), t(Gt))
 
 
 y=yeastODE(parameter,u);
 y=t(y)
 
 return(y)
}


