
#library(Rfast)




Compare.tr.vec <- function() {}
#### Comparing tr and vec operations. 

tr <- function(A=NULL) { sum(diag(A)) }


pred = function(gpars=NULL, newdata=NULL){
  n = dim(newdata)
  G = length(gpars$pi)
  
  dens = sapply(1:G, FUN=function(g){logden(gpars[[g]], newdata)})
  
  vals = exp(log(gpars$pi) + t(dens))
  
  vals = rowMaxs(t(vals))
  return(vals)
}


logden <- function(par = NULL, datax = NULL) {
  ## datax is an array with dimension (p1 x p2 x n) 
  
  logc = -prod(par$p)*log(2*pi*par$s2)/2
  
  logDet1 = -par$p[2]/2*( sum(log(par$del1)) + (par$p[1]-par$q[1])*log(par$eta1) )  
  logDet2 = -par$p[1]/2*( sum(log(par$del2)) + (par$p[2]-par$q[2])*log(par$eta2) )  
  
  r = sweep(datax, c(1,2), par$mu, "-")
  
  trterm = -1/(2*par$s2)*tr.term(r=r, par=par)
  
  lval <- logc + logDet1 + logDet2 + trterm 
  lval
}

logden1 = function(par = NULL, datax = NULL) {
  ## datax is an array with dimension (p1 x p2 x n) 
  
  logc = -prod(par$p)*log(2*pi*par$s2)/2
  
  logDet1 = -par$p[2]/2*( sum(log(par$del1)) + (par$p[1]-par$q[1])*log(par$eta1) )  
  logDet2 = -par$p[1]/2*( sum(log(par$del2)) + (par$p[2]-par$q[2])*log(par$eta2) )  
  
  r = datax-par$mu
  
  phi1 = diag(1/par$del1 - 1/par$eta1, nrow=length(par$del1))
  phi2 = diag(1/par$del2 - 1/par$eta2, nrow=length(par$del2))
  
  S1 = par$lam1 %*% phi1 %*% t(par$lam1) +  1/par$eta1*diag(par$p[1])
  S2 = par$lam2 %*% phi2 %*% t(par$lam2) +  1/par$eta2*diag(par$p[2])
  
  trterm = -1/(2*par$s2)*sum(diag(S1 %*% r %*% S2 %*% t(r)))
  
  lval <- logc + logDet1 + logDet2 + trterm 
  lval
}



loglik <- function(par = NULL, datax = NULL) {
  sum(logden(par = par, datax = x))
}  

norm_eign <- function(del=NULL, eta=NULL, p=NULL, q=NULL) { 
  sum.logval = sum( c(log(del), (p-q)*log(eta)) )
  logval = c(log(del) - sum.logval/(2*q), ( log(eta)- sum.logval/(2*(p-q)) ) )
  exp(logval)
}

#Works in all cases.
norm_eig2 = function(del=NULL, eta=NULL, p=NULL)
{
  q = length(del)
  logdet = sum(  c(log(del),(p-q)*log(eta))  )
  logvals = c(  log(del)-logdet/p, log(eta)-logdet/p  )
  exp(logvals)
}


ipar.diag <- function(p=c(3,5), q=c(1,2)) {
  val = list()
  sigma1 = diag(p[1])
  sigma2 = diag(p[2])
  
  val$mu = matrix(rnorm(prod(p), sd=sqrt(0.5)), p[1], p[2])
  
  temp1 = eigen(sigma1)
  val$lam1  = temp1$vectors[,1:q[1],drop=FALSE]
  
  temp2 = eigen(sigma2)
  val$lam2  = temp2$vectors[,1:q[2],drop=FALSE]
  
  tv1 = norm_eign(del=temp1$values[1:q[1]], eta=temp1$values[q[1]+1], p=p[1], q=q[1])
  val$del1 = tv1[1:q[1]]
  val$eta1 = tv1[q[1]+1]
  
  tv2 = norm_eign(del=temp2$values[1:q[2]], eta=temp2$values[q[2]+1], p=p[2], q=q[2])
  val$del2 = tv2[1:q[2]]
  val$eta2 = tv2[q[2]+1]
  
  val$s2 = 1
  val$p  = p
  val$q  = q
  
  val
}

#Initialize q's automatically
rparq <- function(p=c(3,5), thres=0.75) {
  val = list()
  sigma1 = cov(matrix( rnorm(2*p[1]^2), 2*p[1], p[1] )   )
  sigma2 = cov(matrix( rnorm(2*p[2]^2), 2*p[2], p[2] )   )
  
  val$mu = matrix(rnorm(prod(p), sd=0.1), p[1], p[2])
  
  temp1 = eigen(sigma1)
  eigv = temp1$values/sum(temp1$values)
  props = cumsum(eigv)
  q1 = ifelse( which(props >= thres)[1] < p[1], which(props >= thres)[1], p[1]-1)
  
  
  val$lam1  = temp1$vectors[,1:q1,drop=FALSE]
  
  temp2 = eigen(sigma2)
  eigv = temp2$values/sum(temp2$values)
  props = cumsum(eigv)
  
  q2 = ifelse( which(props >= thres)[1] < p[2], which(props >= thres)[1], p[2]-1)
  
  val$lam2  = temp2$vectors[,1:q2,drop=FALSE]
  
  tv1 = norm_eig2(del=temp1$values[1:q1], eta=temp1$values[q1+1], p=p[1])
  val$del1 = tv1[1:q1 ]
  val$eta1 = tv1[q1+1]
  
  tv2 = norm_eig2(del=temp2$values[1:q2], eta=temp2$values[q2+1], p=p[2])
  val$del2 = tv2[1:q2]
  val$eta2 = tv2[q2+1]
  
  
  val$s2 = 1
  val$p  = p
  val$q = c(q1,q2)
  
  val
}
rgparq = function(p=c(3,5), g=2, n=5, data=NULL, v=1, ...){
  val = lapply(1:g, FUN=function(gr){ rparq(p=p, ...) } )
  val$pi = rep(1/g, times=g)
  
  if ( n > 0) {
    if (is.null(data)) stop("n > 0 and data is null; we need to run iterations ")
    else val = EMqn(data=data, gpar0=val, G=g, n=n, v=v)$gpar
  }
  
  return(val)
}

#Specify qs manually
rgpargen <- function(g=1, p=c(3,5), q=c(2,2), v=1, only.mu=FALSE, sigma="random", init=TRUE) {
  #q converted to a 2xG matrix
  if(is.null(q)) q = matrix(rep(c(1,2), times=g), nrow=g, byrow=TRUE)
  if(length(q) == length(p)) q = matrix(rep(q, times=g), nrow=g, byrow=TRUE)
  if(!(dim(q)[1] == g) ) stop("q's must be in the form of a 2xg matrix")
  
  if(init){
    val = lapply(1:g, FUN=function(gp){rpargen(p=p, q=as.numeric(q[gp,]))} )
    val$pi = rep(1/g,g)
  }else{
    val = lapply(1:g, FUN=function(gp){rpar(p=p, q=as.numeric(q[gp,]))} )
    val$pi = rep(1/g,g) 
  }
  return(val)
}	
rpargen <- function(p=c(3,5), q=c(2,2)) {
  val = list()
  sigma1 = cov(matrix(  rnorm(prod(2*p[1]^2), sd=0.1), 2*p[1], p[1]  ))
  sigma2 = cov(matrix(  rnorm(prod(2*p[2]^2), sd=0.1), 2*p[2], p[2]  ))
  
  val$mu = matrix( rnorm(prod(p), mean=10, sd=2), p[1], p[2])
  
  temp1 = eigen(sigma1)
  val$lam1  = temp1$vectors[,1:q[1],drop=FALSE]
  
  temp2 = eigen(sigma2)
  val$lam2  = temp2$vectors[,1:q[2],drop=FALSE]
  
  #Makes determinant equal to 1
  tv1 = norm_eig2(del=temp1$values[1:q[1]], eta=temp1$values[q[1]+1], p=p[1])
  val$del1 = tv1[1:q[1]]
  val$eta1 = tv1[q[1]+1]
  
  tv2 = norm_eig2(del=temp2$values[1:q[2]], eta=temp2$values[q[2]+1], p=p[2])
  val$del2 = tv2[1:q[2]]
  val$eta2 = tv2[q[2]+1]
  
  
  val$s2 = 1
  val$p  = p
  val$q  = q
  val$sigma1 = diag( c(val$del1, rep(val$eta1, times=p[1]-q[1])) )
  val$sigma2 = diag( c(val$del2, rep(val$eta2, times=p[2]-q[2])) )
  
  val
}


#handles qi=1
tr.term <- function(r=NULL, par=NULL) {
  e12 =  array(apply(r, c(2,3), crossprod, y=par$lam1 )^2, c(par$q[1], dim(r)[2:3]))
  f1 = array( apply(r, c(1,3), "%*%", y=par$lam2 ), c(par$q[2], dim(r)[c(1,3)]))
  f12 = f1^2
  f22 = array( apply(f1, c(1,3), "%*%", y=par$lam1 )^2, c(par$q, dim(r)[3]))
  t1 = apply(f22, c(3), function(z, v1, v2) { v1 %*% (z) %*% v2 }, v1=1/par$del1, v2=1/par$del2 )
  t21 = array(apply(e12, 3, rowSums), c(par$q[1], dim(r)[3]))
  t22 = apply(f22, 3, rowSums)
  t2 = as.numeric( (1/par$del1) %*% (t21-t22) )/par$eta2


  t31 = array(apply(f12, 3, rowSums), c(par$q[2], dim(r)[3]))
  t32 = colSums(f22)
  t3 = as.numeric( (1/par$del2) %*% (t31-t32) )/par$eta1

  t43 = colSums(t21, dims=1) 
  t42 = colSums(t31, dims=1) 
  t44 = colSums(colSums(f22, dims=1), dims=1)
  t41 = colSums(colSums(r^2, dims=1), dims=1)
  t4 = as.numeric( t44 - t43 - t42 + t41 )/(par$eta1*par$eta2)

  val = t1+t2+t3+t4
  val
}

combinewk <- function(weights=NULL, label=NULL)	{
  # known is a numeric with 
  # 0 if unknown group membership 
  # 1,2,3,.. for label of known group
  if (!is.null(label)) {
    kw     = label !=0
    for (j in 1:ncol(weights)) weights[kw,j] = (label == j)[kw]
  }  
  return(weights)	
}


MAP <- function(data, gpar, label=NULL) {
  w = weights(data=data, gpar=gpar)
  if (!is.null(label)) w = combinewk(weights=w, label= label)
  z = apply(w, 1, function(z) { z=(1:length(z))[z==max(z)]; return(z[1]) })
  z = as.numeric(z)
  return( z)	
}


weights <- function(data=NULL, gpar=NULL) {
  G = length(gpar$pi)	
  if (G > 1) {
    zlog = matrix(0, nrow=dim(data)[3], ncol=length(gpar$pi))
    for (k in 1:G ) zlog[,k] =  logden(par=gpar[[k]], datax=data)
    w = t(apply(zlog, 1, function(z,wt) { 
      x = z + log(wt)
      x = x - max(x)
      x = exp(x)
      x =  x/sum(x)
      return( x ) 
    }, wt=gpar$pi))
  } else w = matrix(1,nrow=nrow(data), ncol=G)
  return(w)
}


update.mu <- function(par=NULL, datax=NULL, wt=NULL) {
  ## datax is an array with dimension (p1 x p2 x n) 
  if (is.null(wt)) wt = rep(1, n)
  
  n.wt   = sum(wt)
  par$mu = rowSums( sweep(datax, 3, wt, "*",  check.margin = FALSE), dims=2)/n.wt
  
  par
}

EMn.mu <- function(data=NULL, gpar0=NULL, G=2, n=10, label=NULL, start="kmeans", v=1) 
{
  if (is.null(gpar0)) gpar = igpar(data=data, g=G, start=start)
  else gpar  = gpar0
  
  loglik = numeric(n)
  
  for (i in 1:n) {
    w = weightsv(data=data, gpar=gpar, v=v) # weights update
    w = combinewk(w,label=label)            # make known weights equal to 1
    
    for (k in 1:G ) {
      n.wt   = sum(w[,k])
      
      gpar[[k]]$mu = rowSums( 
                      sweep(data, 3, w[,k],"*", check.margin = FALSE), 
                      dims=2
      )/n.wt
      
      gpar$pi = apply(w,2,mean)
    }  
    
    loglik[i] = mix.loglik(data, gpar)
  }
  
  val = list(loglik=loglik, 
             gpar=gpar, 
             z=weights(data=data, gpar=gpar), 
             map=MAP(data=data, gpar=gpar, label=label)
  )
  
  return(val)
}


#Add emEM initialization to this function
EMn <- function(data=NULL, gpar0=NULL, G=2, q=c(2,3), n=10, label=NULL, v=1, ...) {
  p = dim(data)[1:2]
  
  if (is.null(gpar0)) gpar = rgpar(p=p,  q=q, data=data, g=G, n=10, 
                                   only.mu=FALSE, sigma="identity", v = 1/prod(p))
  else gpar  = gpar0
  
  loglik = numeric(n)
  for (i in 1:n) {
    w = weightsv(data=data, gpar=gpar, v=v)
    w = combinewk(w,label=label)
    for (k in 1:G ) gpar[[k]] = update.par0(par=gpar[[k]], datax=data, wt=w[,k])
    
    gpar$pi = apply(w,2,mean)
    loglik[i] = mix.loglik(data, gpar)
  }
  
  val = list(loglik=loglik, 
             gpar=gpar, 
             z=weights(data=data, gpar=gpar), 
             map=MAP(data=data, gpar=gpar, label=label) 
  )
  
  return(val)
}


#Add emEM initialization to this function
EMqn <- function(data=NULL, gpar0=NULL, G=2, n=10, label=NULL, v=1, ...) {
  p = dim(data)[1:2]
  
  if (is.null(gpar0)) gpar = rgparq(p=p, data=data, g=G, n=5)
  else gpar = gpar0
  
  loglik = numeric(n)
  
  for ( i in 1:n ) {
    w = weightsv(data=data, gpar=gpar, v=v)
    w = combinewk(w,label=label)
    for ( k in 1:G ) gpar[[k]] = update.parq(par=gpar[[k]], datax=data, wt=w[,k], ...)
    
    gpar$pi = apply(w,2,mean)
    loglik[i] = mix.loglik(data, gpar)
  }
  
  val = list(loglik=loglik, 
             gpar=gpar, 
             z=weights(data=data, gpar=gpar), 
             map=MAP(data=data, gpar=gpar, label=label)
  )
  
  return(val)
}


#Allows for different qs across groups
rgpar <- function(g=2, p=c(3,5), q=NULL, n=10, v=1, only.mu=TRUE, sigma="identity", data=NULL ) {
  #q must be 2xg matrix
  if(is.null(q)) q = matrix(rep(c(2,2), times=g), nrow=g, byrow=TRUE)
  if(length(q) == length(p)) q = matrix(rep(q, times=g), nrow=g, byrow=TRUE)
  if(!(dim(q)[1] == g) ) stop("q's must be in the form of a 2xg matrix")

  #stop()
  
  val = list()
  if (sigma == "identity") {
    val = lapply(1:g, FUN=function(gp) {ipar.diag( p=p, q=as.numeric(q[gp,]) )} )
  } else {
    val = lapply(1:g, FUN=function(gp){rpar( p=p, q=as.numeric(q[gp,]) )} )
  }
  val$pi = rep(1/g,g)

  if ( n > 0) {
    if (is.null(data)) stop("n > 0 and data is null; we need to run iterations ")
    if (only.mu) val = EMn.mu(data=data, gpar0=val, G=g, n=n, v=v)$gpar
    else val = EMn(data=data, gpar0=val, G=g, n=n, v=v)$gpar
  }
  return(val)
}

rpar <- function(p=c(3,5), q=c(2,2)) {
  val = list()
  sigma1 = cov( matrix( rnorm(prod(2*p[1]^2)), 2*p[1], p[1] ) )
  sigma2 = cov( matrix( rnorm(prod(2*p[2]^2)), 2*p[2], p[2] ) )
  
  val$mu = matrix(rnorm(prod(p), sd=0.1), p[1], p[2])
  
  temp1 = eigen(sigma1)
  val$lam1  = temp1$vectors[,1:q[1],drop=FALSE]
  
  temp2 = eigen(sigma2)
  val$lam2  = temp2$vectors[,1:q[2],drop=FALSE]
  
  tv1 = norm_eig2(del=temp1$values[1:q[1]], eta=temp1$values[q[1]+1], p=p[1])
  val$del1 = tv1[1:q[1]]
  val$eta1 = tv1[q[1]+1]
  
  tv2 = norm_eig2(del=temp2$values[1:q[2]], eta=temp2$values[q[2]+1], p=p[2])
  val$del2 = tv2[1:q[2]]
  val$eta2 = tv2[q[2]+1]
  
  
  val$s2 = 1
  val$p  = p
  val$q  = q
  
  val
}


wgpar <- function(data, g=2, p=c(3,5), q=c(1,2), w=NULL) {
  if (is.null(w) ) w = matrix(1/g, nrow=dim(data)[3], ncol=g)
  val = list()
  for (k in 1:g) val[[k]] = update.par0(rpar(p,q), data, w[,k])
  val$pi = rep(1/g,g)
  return(val)
}	

mix.loglik <- function(data, gpar) {
  
  logz = matrix(0, nrow=dim(data)[3], ncol=length(gpar$pi))
  
  for (k in 1:length(gpar$pi) ) logz[,k] = logden(par=gpar[[k]], datax=data)
  
  val = sum(apply(logz,1,function(z,wt=NULL) {
    z = z + log(wt)
    logmax = max(z)
    z = z - logmax
    return( log(sum(exp(z))) + logmax )  
  },wt=gpar$pi))
  return(val)
} 





#Allows for qi=1
update.par0 <- function(par=NULL, datax=NULL, wt=NULL) {
  ## datax is an array with dimension (p1 x p2 x n) 
  if (is.null(wt)) wt = rep(1, n)
  
  n.wt   = sum(wt)
  par$mu = rowSums( sweep(datax, 3, wt, "*",  check.margin = FALSE), dims=2)/n.wt
  
  r = sweep(datax, c(1,2), par$mu, "-",  check.margin = FALSE)
  r = sweep(r, 3, sqrt(wt), "*",  check.margin = FALSE)
  
  par$s2 = sum(tr.term(r=r, par=par))/(prod(par$p)*n.wt )
  
  B1 = par$lam2 %*% diag(1/par$del2 - 1/par$eta2, nrow=length(par$del2)) %*% t(par$lam2) +  1/par$eta2*diag(par$p[2])
  A1 = diag(par$p[1]) - par$lam1 %*% t(par$lam1)
  
  ceta1 = mean(apply(r, 3, function(z, A=NULL, B=NULL) { sum( (t(z) %*% A %*% (z) )*B) }, A=A1, B=B1 ))/( par$p[1] - par$q[1])
  
  
  cdel1 = rowMeans(
    array(
      apply(r, 3, function(z, A=NULL, B=NULL) { diag( t(A) %*% z %*% B %*% t(z) %*% A) }, A=par$lam1, B=B1 ),
      c(par$q[1], dim(r)[3]) #dims of array
    )
  )
  
  par$del1 = cdel1 * ( prod(cdel1) * (ceta1)^(par$p[1] - par$q[1]) )^(-1/ par$p[1] )  #del1 update
  par$eta1 = ceta1 * ( prod(cdel1) * (ceta1)^(par$p[1] - par$q[1]) )^(-1/ par$p[1] )  #eta1 update
  
  B0 = diag(par$p[2]) - par$lam2 %*% t(par$lam2)
  A0 = par$lam1 %*% diag(1/par$del1 - 1/par$eta1, nrow=length(par$del1)) %*% t(par$lam1) + 1/par$eta1*diag(par$p[1])
  ceta2 = mean(apply(r, 3, function(z, A=NULL, B=NULL) { sum( (z %*% B %*% t(z) )*A) }, B=B0, A=A0 ))/( par$p[2] - par$q[2])
  cdel2 = rowMeans(
    array(
      apply(r, 3, function(z, A=NULL, B=NULL) { diag( t(B) %*% t(z) %*% A %*% (z) %*% B) }, B=par$lam2, A=A0 ),
      c(par$q[2], dim(r)[3])
    )
  )
  
  par$del2 = cdel2 * ( prod(cdel2) * (ceta2)^(par$p[2] - par$q[2]) )^(-1/ par$p[2] )  
  par$eta2 = ceta2 * ( prod(cdel2) * (ceta2)^(par$p[2] - par$q[2]) )^(-1/ par$p[2] ) 
  
  B1 = par$lam2 %*% diag(1/par$del2 - 1/par$eta2, nrow=length(par$del2)) %*% t(par$lam2) +  1/par$eta2*diag(par$p[2])
  W1= rowMeans(apply(r, 3, function(z, B=NULL) { z %*% B %*% t(z)  }, B=B1 ))
  dim(W1) = rep(par$p[1],2)
  par$lam1 = newLam1.MM(lam=par$lam1, W=W1, A = 1/par$del1 -1/par$eta1 ) 
  
  A1 = par$lam1 %*% diag(1/par$del1 - 1/par$eta1, nrow=length(par$del1)) %*% t(par$lam1) +  1/par$eta1*diag(par$p[1])
  W2= rowMeans(apply(r, 3, function(z, A=NULL) { t(z) %*% A %*% (z)  }, A=A1 ))
  dim(W2) = rep(par$p[2],2)
  par$lam2 = newLam1.MM(lam=par$lam2, W=W2, A = 1/par$del2 -1/par$eta2 ) 
  
  par
}


#Allows for qi=1
update.parq <- function(par=NULL, datax=NULL, wt=NULL, thres=0.75) {
  ## datax is an array with dimension (p1 x p2 x n) 
  if (is.null(wt)) wt = rep(1, n)
  
  n.wt   = sum(wt)
  par$mu = rowSums( sweep(datax, 3, wt, "*",  check.margin = FALSE), dims=2)/n.wt
  
  r = sweep(datax, c(1,2), par$mu, "-",  check.margin = FALSE)
  r = sweep(r, 3, sqrt(wt), "*",  check.margin = FALSE)
  
  par$s2 = sum(tr.term(r=r, par=par))/(prod(par$p)*n.wt )
  
  #invpsi2
  B1 = par$lam2 %*% diag(1/par$del2 - 1/par$eta2, nrow=length(par$del2)) %*% t(par$lam2) +  1/par$eta2*diag(par$p[2])
  A1 = diag(par$p[1]) - par$lam1 %*% t(par$lam1)
  
  W2  = matrix(
    rowmeans( apply( r, 3, function(z, A=NULL){ z%*%tcrossprod(A,z) }, A=B1 )), 
    nrow=par$p[1] 
  )
  
  #Get eigenvalues of W2
  eigv = eigen(W2)$values/sum(eigen(W2)$values)
  props = cumsum(eigv)
  
  #Choose q1 using thresholding
  par$q[1] = ifelse( which(props >= thres)[1] < par$p[1], which(props >= thres)[1], par$q[1])
  
  #estimate eigenvalues delta1
  ceta1 = mean(apply(r, 3, function(z, A=NULL, B=NULL) { sum( (t(z) %*% A %*% (z) )*B ) }, A=A1, B=B1 ))/(par$p[1] - par$q[1])
  cdel1 = rowMeans(
    array(
      apply(r, 3, function(z, A=NULL, B=NULL) { diag( t(A) %*% z %*% B %*% t(z) %*% A) }, A=par$lam1, B=B1 ),
      c(par$q[1], dim(r)[3]) #dims of array
    )
  )
  
  par$del1 = cdel1 * ( prod(cdel1) * (ceta1)^(par$p[1] - par$q[1]) )^(-1/ par$p[1] )  #del1 update
  par$eta1 = ceta1 * ( prod(cdel1) * (ceta1)^(par$p[1] - par$q[1]) )^(-1/ par$p[1] )  #eta1 update
  
  par$lam1 = newLam1.MM(lam=par$lam1, W=W2, A = 1/par$del1 -1/par$eta1, q=par$q[1] ) 

  
  #invpsi1
  A0 = par$lam1 %*% diag(1/par$del1 - 1/par$eta1, nrow=length(par$del1)) %*% t(par$lam1) + 1/par$eta1*diag(par$p[1])
  B0 = diag(par$p[2]) - par$lam2 %*% t(par$lam2)
  W1 = matrix(
    rowmeans( apply( r, 3, function(z, A=NULL){ t(z)%*%A%*%z }, A=A0 )), 
    nrow=par$p[2] 
  )
  
  #Get eigenvalues of W1
  eigv = eigen(W1)$values/sum(eigen(W1)$values)
  props = cumsum(eigv)
  
  #estimate q2 with thresholding
  par$q[2] = ifelse( which(props >= thres)[1] < par$p[2], which(props >= thres)[1], par$q[2])
  

  #estimate eigenvalues delta2
  ceta2 = mean(apply(r, 3, function(z, A=NULL, B=NULL) { sum( (z %*% B %*% t(z) )*A) }, B=B0, A=A0 ))/( par$p[2] - par$q[2])
  cdel2 = rowMeans(
    array(
      apply(r, 3, function(z, A=NULL, B=NULL) { diag( t(B) %*% t(z) %*% A %*% (z) %*% B) }, B=par$lam2, A=A0 ),
      c(par$q[2], dim(r)[3])
    )
  )
  
  par$del2 = cdel2 * ( prod(cdel2) * (ceta2)^(par$p[2] - par$q[2]) )^(-1/ par$p[2] )  
  par$eta2 = ceta2 * ( prod(cdel2) * (ceta2)^(par$p[2] - par$q[2]) )^(-1/ par$p[2] ) 
  
  
  par$lam2 = newLam1.MM(lam=par$lam2, W=W1, A = 1/par$del2 -1/par$eta2, q=par$q[2] ) 
  
  par
}



weightsv <- function(data=NULL, gpar=NULL, v=1) {
  G = length(gpar$pi)	
  if (G > 1) {
    zlog = matrix(0, nrow=dim(data)[3], ncol=length(gpar$pi))
    for (k in 1:G ) zlog[,k] =  logden(par=gpar[[k]], datax=data)
    w = t(apply(zlog, 1, function(z,wt) { 
      x = z + log(wt)
      x = x - max(x)
      x = exp(v*x)
      x =  x/sum(x)
      return( x ) 
    }, wt=gpar$pi))
  } else w = matrix(1,nrow=nrow(data), ncol=G)
  return(w)
}

newLam1.MM <- function(lam=NULL, W=NULL, A=NULL, q=NULL) {
  p = nrow(lam)
  lam.new = eigen(W, symmetric=TRUE)$vectors[,1:q, drop=FALSE]
  #lam.new = eigen(-1*W, symmetric=TRUE)$vectors[, p:(p-q+1), drop=FALSE]
  return( lam.new )
}

newLam1.MM.old <- function(lam=NULL, W=NULL, A=NULL ) {
  alpha = max(eigen(W, only.values = TRUE, symmetric=TRUE)$values)
  #z = diag(A) %*% t(lam) %*% W  - alpha *(diag(A)%*% t(lam) )
  z =  (A * t(lam)) %*% W  - alpha *(A * t(lam) )
  
  z1 = svd(z)
  Xk1 = (z1$v) %*% t(z1$u) 
  return( Xk1 )
}



## Combining function for foreach to split results into matrices
custcomb = function(r1=NULL, r2=NULL)
{
  val = list("diff" = rbind(r1[["diff"]], r2[["diff"]]),
             "map"  = rbind(r1[["map"]],  r2[["map"]]))
  return(val)
}  

