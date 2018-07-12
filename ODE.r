# Differentials - change in time

# load packages
require(deSolve)

# define the equation - simple population growth
LogisticGrowth <- function(t, y, parms) {
  N <- y;
  dN <- parms[1]*N*(1-N/parms[2]);
  return(list(dN));  # return as list for when returning multiple outcomes
}

# set up the model
  x0 <- 10 # initial condition
  r <- 1 # growth rate
  k <- 100 # carrying capacity
  times <- seq(0,10, by=.1) # vector of time values
  
# solve the equation for how population changes over time  
# starting with an initial condition of x0
  out <- lsoda(x0, times, LogisticGrowth, parms=c(r,k)) 
  # out here is a matrix, with the first col a time, and each subsequent column the other outcomes
  plot(out[,1], out[,2], type='l', xlab="time", ylab="population")
  
# Lotka Volterra Predator Prey model, with carrying capacity
  LotkaVolteraPredPrey = function (t, y, parms){
    h = y[1]; p = y[2] # herbivores(prey), predators(p)
    dh = parms[1]*h*(1-h/parms[5]) - parms[2]*h*p
    dp = -parms[3]*p + parms[4]*h*p
    return( list( c(dh, dp) ) )
  }
  
  h0 <- 2
  p0 <- 2
  r <- 1
  c <- 1
  m <- 5
  a <- 0.5
  k <- 50
  times <- seq(0,20, by=.2)
  
  out.lv <- lsoda(c(h0,p0), times, LotkaVolteraPredPrey, parms=c(r,c,m,a,k)) # out here is a matrix
  plot(out.lv[,1], out.lv[,2], type="l", xlab="time (t)", ylab= "Population Size", ylim=c(0,max(out.lv[,2:3])))
  lines (out.lv[,1], out.lv[,3], col="red")
  legend("topright", c("prey (h)", "predator (p)"), lty=1, col=c("black","red"))
  
  plot(out.lv[,2], out.lv[,3], xlab="prey", ylab="predator", type="l")
  
  
# Gene Regulatory Network - "Toggle-switch model"  
  Toggle = function(t,y,parms) {
    u=y[1]; v=y[2];
    du = -u + parms[1]/(1+v^parms[2]);
    dv = -v + parms[3]/(1+u^parms[4]);
    dY = c(du, dv);
    return(list(dY));
  }
  
  u0 <- .6
  v0 <- .8
  a1 <- 2
  a2 <- 2
  b <- 3
  g <- 3
  times <- seq(0,20, by=.2)
  out <- lsoda(c(u0,v0), times, Toggle, parms=c(a1,a2,b,g)) # out here is a matrix
  
  plot(out[,1], out[,2], type="l", xlab="time (t)", ylab= "Protein conc.", ylim=c(0,max(out[,2:3])))
  lines(out[,1], out[,3], col="red")
  
# SIR models
  SIR <- function(t, y, parms) {
    s=y[1]; i=y[2]; r=y[3]; # susceptible, infected, recovered
    ds= -parms[1]*s*i;
    dr= parms[2]*i;
    di= parms[1]*s*i - parms[2]*i;
    dY = c(ds, di, dr);
    return(list(dY));
  }
  
  S0 <- 7900000
  I0 <- 10
  R0 <- 0
  n <- S0+I0+R0 # total population size
  s0 <- S0/n
  i0 <-I0/n
  r0 <- R0/n
  b <- 1/2 # contacts per day
  k <- 1/3 # fraction of the infected group that will recover every day.  
  times <- seq(0,150,by=1)
    
  out <- lsoda(c(s0,i0,r0), times, SIR, parms=c(b,k))
  
  plot(out[,1], out[,2], type="l", xlab="time (t)", ylab= "prop.", ylim=c(0,max(out[,2:4])), col="blue")
  lines(out[,1], out[,3], col="green")
  lines(out[,1], out[,4], col="red")
