ptw=function(x){

#rg <- as.matrix(read.table("cont_ranges",row.names=1))
rg=matrix(c(113.0,575.0,0.0001,3900.0, 0.0,0.4,0.027,0.054,0.53, 3.1, 0.11, 0.19,12.0,16.0,0.0011,0.002,0.0000014,0.00072,0.0065,0.011,0.000097,     0.0017),ncol=2,byrow=T)

rownames(rg)=c('temp','strain_rate','strain', 'ptw_be_theta0' ,  'ptw_be_p' ,'ptw_be_kappa','ptw_be_mloggamma','ptw_be_y0','ptw_be_yinf','ptw_be_s0','ptw_be_sinf')


des=as.matrix(x)

#des <- as.matrix(read.table("s-lhs.120.11"))
for (ii in 1:nrow(rg)) {
  des[,ii] <- (rg[ii,2]-rg[ii,1])*des[,ii]+rg[ii,1]
}
colnames(des) <- rownames(rg)

iis <- which(des[,"ptw_be_yinf"] > des[,"ptw_be_sinf"])
for (ii in iis) {
  yinf <- runif(1,min=rg["ptw_be_yinf",1],max=des[ii,"ptw_be_sinf"])
  des[ii,"ptw_be_yinf"] <- yinf
} 

des=signif(des,6)

# Constants

barperdyne <- 1e6
avogadro <- 6.023e23

# Defined quantities (Be)

G0 <- 1524
Tm <- 1560
rho <- 1.85
W <- 9.013
alphaP <- 0.32

# Error function

erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

# Compute tau.s and tau.y

ctau <- function(a,b,c,d,psidot,xidot) {
   if (log(b) >= log(psidot/xidot)) {
      c-(c-d)*erf(a*That*log(b*xidot/psidot))
   } else {
      NaN
   }
}

# Compute ptw output

ptw.o <- function(psi, ptw) {
   theta0 <- ptw[1]; p <- ptw[2]; kappa <- ptw[3]; gamma <- ptw[4];
   y0 <- ptw[5]; yinf <- ptw[6]; s0 <- ptw[7]; sinf <- ptw[8];
   y1 <- ptw[9]; y2 <- ptw[10]; beta <- ptw[11];

   G0 <- G0*1000; G0 <- G0*barperdyne;
   G <- G0 * (1 - alphaP*That)
   xidot <- 0.5*((4*pi*rho*avogadro)/(3*W))^(1/3)*sqrt(G/rho)

# low strain rate

   tauhatyl <- ctau(kappa,gamma,y0,yinf,psidot,xidot)
   tauhatsl <- ctau(kappa,gamma,s0,sinf,psidot,xidot)

# medium strain rate

   tauhatym <- y1*(psidot/(gamma*xidot))^y2

# high strain rate

   tauhatsh <- s0*(psidot/(gamma*xidot))^beta
   tauhatyh <- tauhatsh

   tauhats <- max(tauhatsl,tauhatsh)
   tauhaty <- max(tauhatyl, min(tauhatym,tauhatyh))   

# stress calculation

   term1 <- s0 - tauhaty
   term2 <- tauhats - tauhaty
   term3 <- 1 - exp((-p*term2)/term1)
   term4 <- (s0 - tauhaty)*(exp(p*term2/term1)-1)
   if (term4 != 0) {
      term5 <- term3*exp((-p*theta0*psi)/term4)
      tauhat <- ifelse(term5 != 1, tauhats + term1*log(1-term5)/p,
                       -Inf)
   } else if (p != 0) {
      tauhat <- tauhats
   } else {
      tauhat <- tauhats - term2*exp((-theta0*psi)/term2)
   }
   tau <- tauhat*G
   list(psi=psi,stress=2*tau*1e-7,tauhat=tauhat,tauhaty=tauhaty,
        tauhats=tauhats)
}

# Read in design
#des <- as.matrix(read.table("design.txt"))
des[,7] <- exp(-des[,7])
ptw.fix <- c(0.012, 0.4, 0.23)

sim <- numeric(nrow(des))
for (ii in 1:nrow(des)) {
  # Temperature/strain rate/strain specification

  T <- des[ii,1]
  psidot <- des[ii,2]
  psi <- des[ii,3]

  # Compute That

  That <- T/Tm

  ptw <- c(des[ii,-c(1:3)],ptw.fix)

  sim[ii] <- ptw.o(psi, ptw)$stress
}

#write.table(signif(sim,6),"../sim_outputs",quote=F,row.names=F,
#            col.names=F)
return(sim)  
}