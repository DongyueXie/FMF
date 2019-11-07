
update_smooth = function(x,sf,reflect=T,shape=NULL,point_mass=T,nullweight=100,weight,
                    g_init = NULL, fix_g = FALSE, control =  NULL,
                    return_pos_only=F){
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
  post_mean_log = matrix(0,dim(nt)[1],dim(nt)[2])
  post_mean_log1_p = matrix(0,dim(nt)[1],dim(nt)[2])
  pi_weights = matrix(0,dim(nt)[1],length(shape)+point_mass)

  for(s in 1:dim(ns)[1]){
    fit = ebbp_beta_mixture(ns[s,],nt[s,],shape,point_mass,nullweight,weight=rep(1,dim(nt)[2]),
                            g_init, fix_g, control,
                            return_pos_only)
    loglik = loglik + fit$log_likelihood
    post_mean[s,] = fit$posterior$mean
    post_mean_log[s,] = fit$posterior$mean_log
    post_mean_log1_p[s,] = fit$posterior$mean_log1_p
    if(point_mass){
      pi_weights[s,] = c(fit$fitted_g$pi0,fit$fitted_g$pi)
    }else{
      pi_weights[s,] = fit$fitted_g$pi
    }
  }

  # get estimate of lambda
  est=reverse_pp_log(tit,post_mean,post_mean_log,post_mean_log1_p,sf)

  est2 = est$est

  results = list("E" = est2[(pos_to_fill_1+1):(pos_to_fill_1 + length(x))],
                 'Elog' = log(est2[(pos_to_fill_1+1):(pos_to_fill_1 + length(x))]),
                 #'Elogf' = est$est_log[(pos_to_fill_1+1):(pos_to_fill_1 + length(x))],
                 "pi_weights" = pi_weights,
                 "loglik" = loglik)

  return(results)

}

reverse_pp_log=function(dmat,pm,pmlogp,pmlog1_p,s){
  n=dim(dmat)[2]
  J=log2(n)
  # Beginning with the total number of counts, working from coarse
  # scales (i.e., deep depths) upwards, gradually build up the MSPB
  # shrinkage estimate by multiplying by appropriate shrinkage factors
  # at each level.
  est = dmat[J+1,]/s
  est_log = log(dmat[J+1,]/s)
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

      ss_logp = pmlogp[D,ind[1:nDo2]]
      ss_lop1_p = pmlog1_p[D,ind[1:nDo2]]
      nestl_log = interleave(ss_logp,ss_lop1_p)



      # In the second half of the vector of D+1-depth counts,
      # we right-shift the D-depth counts, compute the shrunken
      # values, and left-shift these values back to the order
      # of the above.
      estr = estvec[(nDo2+1):nD]
      ss=pm[D,ind[(nDo2+1):nD]]
      nestr = interleave(estr*ss,estr*(1 - ss))
      nestr = lshift(nestr)

      ss_logp = pmlogp[D,ind[(nDo2+1):nD]]
      ss_lop1_p = pmlog1_p[D,ind[(nDo2+1):nD]]
      nestr_log = interleave(ss_logp,ss_lop1_p)
      nestr_log = lshift(nestr_log)

      # Combine the estimates from both halves of the D+1-depth
      # counts, and store.
      est[ind] = 0.5*( nestl + nestr )
      est_log[ind] = 0.5*(nestl_log+nestr_log)
    }
  }
  return(est=list(est=est,est_log=est_log))
}

interleave=function(x,y){
  return(as.vector(rbind(x,y)))
}
