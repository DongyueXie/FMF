# Bayesian Multiscale models for Poisson Process
BMSM = function(x,reflect=T,shape=NULL,point_mass=T,nullweight=100,weight,
                g_init = NULL, fix_g = FALSE, control =  NULL,
                return_pos_only=F,s=1){
  if(min(x) < 0){stop ("negative values in x not permitted")}
  if(is.null(shape)){
    shape = c(100, 50, 20, 10, 5, 2, 1)
  }

  # whether reflect data so that it has a length of powers of 2

  if(reflect){
    extended_len <- 2^{ceiling(log(length(x), base=2))}
    if(extended_len > length(x)){
      pos_to_fill <- extended_len - length(x)
      pos_to_fill_1 <- floor(pos_to_fill/2)
      pos_to_fill_2 <- pos_to_fill - pos_to_fill_1
      if(pos_to_fill_1 >= 1){
        x_ext <- c(rev(head(x, pos_to_fill_1)), x, rev(tail(x, pos_to_fill_2)))
      }else{
        x_ext <- c(x, rev(tail(x, pos_to_fill_2)))
      }
    }else if(extended_len == length(x)){
      pos_to_fill_1 <- 0
      x_ext <- x
    }else{
      stop("error in extending the vector to make its size a power of 2")
    }
  }else{
    extended_len <- 2^{ceiling(log(length(x), base=2))}
    if(extended_len > length(x)){
      pos_to_fill <- extended_len - length(x)
      pos_to_fill_1 <- floor(pos_to_fill/2)
      pos_to_fill_2 <- pos_to_fill - pos_to_fill_1
      if(pos_to_fill_1 >= 1){
        x_ext <- c(rep(head(x,1), pos_to_fill_1), x, rep(tail(x, 1), pos_to_fill_2))
      }else{
        x_ext <- c(x, rep(tail(x, 1), pos_to_fill_2))
      }
    }else if(extended_len == length(x)){
      pos_to_fill_1 <- 0
      x_ext <- x
    }else{
      stop("error in extending the vector to make its size a power of 2")
    }
  }

  # construct Translation-Invariant table

  x_odd <- x_ext[c(TRUE,FALSE)]
  x_even <- x_ext[c(FALSE, TRUE)]

  titable=ParentTItable(x_ext)

  tit=titable$TItable
  ptit=titable$parent
  n=dim(tit)[2]
  J=dim(tit)[1]-1
  # ns: left child
  # nf: right child
  nt=tit[-1,]
  ns=ptit[,((1:(2*n))%%2==1)]
  nf=nt-ns
  loglik = 0

  # for each scale, perform EB shrinkage.

  post_mean = matrix(0,dim(nt)[1],dim(nt)[2])
  #post_mean_log = matrix(0,dim(nt)[1],dim(nt)[2])
  pi_weights = matrix(0,dim(nt)[1],length(shape)+point_mass)

  for(s in 1:dim(ns)[1]){
    fit = ebbp_beta_mixture(ns[s,],nt[s,],shape,point_mass,nullweight,weight=rep(1,dim(nt)[2]),
                            g_init, fix_g, control,
                            return_pos_only)
    loglik = loglik + fit$log_likelihood
    post_mean[s,] = fit$posterior$mean
    #post_mean_log[s,] = fit$posterior$mean_log
    if(point_mass){
      pi_weights[s,] = c(fit$fitted_g$pi0,fit$fitted_g$pi)
    }else{
      pi_weights[s,] = fit$fitted_g$pi
    }
  }

  # get estimate of lambda
  est=reverse.pp(tit,post_mean)

  est2 = est - mean(est) + mean(x_ext)

  results = list("estimate" = est2[(pos_to_fill_1+1):(pos_to_fill_1 + length(x))],
             "full_est" = est2,
             "pi_weights" = pi_weights,
             "loglik" = loglik,
             'x_ext' = x_ext)

  return(results)


}




#############  Shift operator functions  ##########################

rshift = function(x){L=length(x); return(c(x[L],x[-L]))}

lshift = function(x){return(c(x[-1],x[1]))}


##############  Parent TI table maker (R version)  ######################

ParentTItable=function(sig){
  n = length(sig)
  J = log2(n)

  # Create decomposition table of signal, using pairwise sums,
  # keeping just the values that are *not* redundant under the
  # shift-invariant scheme.  This is very similar to TI-tables
  # in Donoho and Coifman's TI-denoising framework.
  dmat = matrix(0, nrow=J+1, ncol=n)
  dmat[1,] = sig
  #dmat[1,] = as.matrix(sig)
  dmat2 = matrix(0, nrow=J, ncol=2*n) #the parent table

  for(D in 0:(J-1)){
    nD = 2^(J-D);
    nDo2 = nD/2;
    twonD = 2*nD;
    for(l in 0:(2^D-1)){
      ind = (l*nD+1):((l+1)*nD)
      ind2 = (l*twonD+1):((l+1)*twonD)
      x = dmat[D+1,ind]
      lsumx = x[seq(from=1,to=nD-1, by=2)] + x[seq(from=2,to=nD,by=2)]
      rx = rshift(x);
      rsumx = rx[seq(from=1,to=nD-1, by=2)] + rx[seq(from=2,to=nD,by=2)]
      dmat[D+2,ind] = c(lsumx,rsumx)
      dmat2[D+1,ind2] = c(x,rx)
    }
  }
  return(list(TItable=dmat,parent=dmat2))
}


# dmat: TI table
# pp: prior weight
# qq: shape

#now
# pm: post_mean
reverse.pp=function(dmat,pm){
  n=dim(dmat)[2]
  J=log2(n)
  # Beginning with the total number of counts, working from coarse
  # scales (i.e., deep depths) upwards, gradually build up the MSPB
  # shrinkage estimate by multiplying by appropriate shrinkage factors
  # at each level.
  est = dmat[J+1,]
  for (D in J:1){
    nD = 2^(J-D+1)
    nDo2 = nD/2
    for (l in 0:(2^(D-1)-1)){
      # Set indexing so as to pick off blocks of size 2^(J-D+1)
      # when shrinking estimates at depth D+1 down to finer
      # scale at depth D.
      ind = (l*nD+1):((l+1)*nD)
      estvec = est[ind]

      # In the first half of the vector of D+1-depth estimates,
      # we can shrink using the D-depth counts in the order
      # in which they appear.
      estl = estvec[1:nDo2]
      ss = pm[D,ind[1:nDo2]]
      nestl = interleave(estl*ss,estl*(1 - ss))

      # In the second half of the vector of D+1-depth counts,
      # we right-shift the D-depth counts, compute the shrunken
      # values, and left-shift these values back to the order
      # of the above.
      estr = estvec[(nDo2+1):nD]
      ss=pm[D,ind[(nDo2+1):nD]]
      nestr = interleave(estr*ss,estr*(1 - ss))
      nestr = lshift(nestr)

      # Combine the estimates from both halves of the D+1-depth
      # counts, and store.
      est[ind] = 0.5*( nestl + nestr )
    }
  }
  return(est)
}

interleave=function(x,y){
  return(as.vector(rbind(x,y)))
}
