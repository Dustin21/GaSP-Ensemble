source("borehole_func.R")

#Set Dimension for problem
d = 8
#Set the function Name here 
funcname = 'bore.function'
  
#Read in the base design.
xf = read.table('XD8N40.txt',header=F)   
colnames(xf) = paste('x',1:d,sep='')
  
  
#Compute test set data
xp = data.matrix(read.table('test8.csv',header=F,sep=',')[,1:8])
colnames(xp) = colnames(xf)
yp = data.frame(y=do.call(funcname,list(xp)))

####### Permutation on the base design
ocn = colnames(xf)
cp.m = data.frame(read.csv("order_borehole.csv",  h=F)[, 1:8])

#####for the i~th permutation, i=1, 2,...., 25. 
cp = as.numeric(cp.m[i, ]) #i=1, 2, 3,....., 25. 
#####Can do a for loop, such as for(i in 1:25)

####actual training set
xf_train = data.matrix(xf[,cp])
colnames(xf_train) = ocn  
yf_train = data.frame(y=do.call(funcname,list(xf_train)))

########xf_train, yf_train: training X and training Y. 
########xp and yp: test X and true output Y. 
######## can call mlegp, etc. 
library(mlegp)
fit <- mlegp(X=xf_train, Z=yf_train)
yp_pred <- predict(fit, xp) 
( RMSE<- sqrt(sum((yp_pred-yp)^2)/length(yp_pred)) )



