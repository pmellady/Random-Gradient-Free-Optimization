#Booth function
f<-function(x){
  x1<-x[1]
  x2<-x[2]
  f_x<-(x1+2*x2-7)^2+(2*x1+x2-5)^2
  return(f_x)
}

gradf<-function(x){
  x1<-x[1]
  x2<-x[2]
  return(c(10*x1+8*x2-34, 8*x1+10*x2-38))
}

hess_f<-function(x){
  hf<-matrix(c(10,8,8,10), byrow=TRUE, ncol=2)
  return(hf)
}

# Gradient free random method ##################################################

g_mu<-function(x, u, mu){
  gmu_k<-(1/mu)*(f(x+mu*u)-f(x))*u
  return(gmu_k)
}

RG_mu<-function(guess, K, sol, epsilon=1e-3){
  n<-length(guess)
  
  H<-hess_f(x0)
  L<-max(eigen(H)$value)
  mu<-(3/(3*(n+4)))*sqrt((epsilon)/(2*L))
  
  x_k<-guess
  h<-1/(4*(n+4)*L)
  k<-0
  x_trace<-c(f(x_k), x_k)
  while(k<K){
    u_k<-matrix(rnorm(n, 0, 1), ncol=1)
    g_k<-g_mu(x_k, u_k, mu)
    x0<-x_k
    x_k<-x_k-h*g_k
    k<-k+1
    
    x_trace<-rbind(x_trace, c(f(x_k), x_k))
    
    if(sqrt(sum(sol-x_k)^2)<epsilon){
      break
    }
  }
  x_trace<-as.data.frame(x_trace)
  names(x_trace)<-c("f",paste('X', 1:n, sep=" "))
  rownames(x_trace)<-NULL
  return(list(x_k=x_k, x_trace=x_trace, k=k))

}

# Gradient descent for comparison ##############################################

back_track<-function(rho, c, alpha0, x_k, dir){
  alpha_k<-alpha0
  while(f(x_k+alpha_k*dir)>f(x_k)+c*alpha_k*gradf(x_k)%*%dir){
    alpha_k<-alpha_k*rho
  }
  
  return(alpha_k)
}

steep_descent<-function(guess, epsilon=1e-3, rho=.3, c=1e-4, alpha0=1){
  x_k<-guess
  n<-length(x_k)
  k<-0
  x_trace<-c(f(x_k), x_k)
  while(sqrt(sum(gradf(x_k)^2))>epsilon){
    dir<--gradf(x_k)
    alpha<-back_track(rho, c, alpha0, x_k, dir)
    x_k<-x_k+alpha*dir
    k<-k+1
    x_trace<-rbind(x_trace, c(f(x_k), x_k))
    if(alpha<1e-15){
      break
    }
  }
  x_trace<-as.data.frame(x_trace)
  names(x_trace)<-c("f",paste('X', 1:n, sep=" "))
  rownames(x_trace)<-NULL
  return(list(x_k=x_k, x_trace=x_trace, k=k))
}

# Finding the solutions ########################################################
set.seed(405)
x0<-matrix(c(-1.2, 1), ncol=1)
sol<-matrix(c(1,3), ncol=1)
RG<-RG_mu(x0, 1e6, sol)

x0<-c(-1.2,1)
GD<-steep_descent(x0)

# Plot of f vs iteration number for a single value of epsilon
plot(1:RG$k, RG$x_trace[1:RG$k, 1], type="l", xlab="iteration number", ylab="criterion f")
lines(1:GD$k, GD$x_trace[1:GD$k, 1], col="blue")


# Plot of iterations required to reach certain accuracy values
m<-20
epsilon<-seq(1e-4,1, length=m)
RG_ep<-c()
GD_ep<-c()
for(i in 1:m){
  RG_ep<-rbind(RG_ep, c(RG_mu(x0, 1e6, sol, epsilon[i])$k)) 
  GD_ep<-rbind(GD_ep, c(steep_descent(x0, epsilon[i])$k)) 
}
plot(epsilon[m:1], RG_ep[m:1], type="l", xlab="epsilon", ylab="iteration number", ylim=c(0,250))
lines(epsilon[m:1], GD_ep[m:1], col="blue")









