
devtools::load_all('~/Documents/Rpackages/flashr')
devtools::load_all('~/Documents/Rpackages/smashr')
library(wavethresh)
#'@param Y
smooth_flash = function(Y,S,Kmax=1000,tol=0.01,
                        init_fn = "udv_si",
                        ebnm_fn = "ebnm_pn",
                        var_type = 'by_row',
                        ebnm_param = NULL,
                        f_init = NULL,
                        filter.number = 1,
                        family = 'DaubExPhase',
                        nullcheck=FALSE,
                        maxiter=100){

  f = handle_f(f_init, init_null_f = TRUE)
  data = flash_set_data(Y,S)
  data = handle_data(data,f)
  if(!is.null(S)){
    f$tau = 1/S^2
  }
  var_type = handle_var_type(var_type, data)
  init_fn = handle_init_fn(init_fn)
  ebnm_fn = handle_ebnm_fn(ebnm_fn)
  ebnm_param = handle_ebnm_param(ebnm_param, ebnm_fn, Kmax)

  history = list()
  prev_K = flash_get_k(f)

  for(k in 1:Kmax){


    old_f = f
    res = smooth_flash_r1(data,
                   f,
                   var_type,
                   init_fn,
                   tol,
                   ebnm_fn$l,
                   ebnm_param$l[[k]],
                   filter.number,
                   family,
                   nullcheck,
                   maxiter)

    f = res$f
    next_history = res$history

    # Test whether the factor/loading combination is effectively zero.
    if (is_tiny_fl(f, flash_get_k(f))) {
      # Remove zero factor as long as it doesn't create an empty object.
      if (flash_get_k(f) > 1) {
        f = old_f
        next_history$zeroed_out = prev_K + k
      }
      history = c(history, list(next_history))
      break
    }

    history = c(history, list(next_history))

  }

  flash_object = construct_flash_object(data = data,
                                        fit = f,
                                        history = history,
                                        f_init = f_init)

  return(flash_object)

}


smooth_flash_r1 = function(data,
                    f_init,
                    var_type,
                    init_fn,
                    tol,
                    ebnm_fn_l,
                    ebnm_param_l,
                    filter.number,
                    family,
                    nullcheck,
                    maxiter
                    ) {

  f = add_factors_from_data(data, K = 1, f_init, init_fn)

  opt_res = smooth_flash_optimize_single_fl(data,
                                     f,
                                     flash_get_k(f),
                                     var_type,
                                     tol,
                                     ebnm_fn_l,
                                     ebnm_param_l,
                                     filter.number,
                                     family,
                                     maxiter
                                     )

  f = opt_res$f

  if (nullcheck) {
    null_res = perform_nullcheck(data,
                                 f,
                                 flash_get_k(f),
                                 var_type,
                                 verbose = F)
    f = null_res$f
    # zeroed_out field is handled in flash_greedy_workhorse
  }

  return(list(f = f, history = opt_res$history))
}



smooth_flash_optimize_single_fl = function(data,f,k,
                                           var_type,
                                           tol,
                                           ebnm_fn_l,
                                           ebnm_param_l,
                                           filter.number,
                                           family,
                                           maxiter){

  R2 = flash_get_R2(data, f)
  # Expected residuals and squared residuals with factor k excluded:
  Rk = flash_get_Rk(data, f, k)
  R2k = (R2 + 2 * outer(f$EL[, k], f$EF[, k]) * Rk
         - outer(f$EL2[, k], f$EF2[, k]))

  iter = 0
  diff = Inf
  diff_track = rep(NA, maxiter)
  obj_track = rep(NA, maxiter)

  while ((iter < maxiter) && (diff > tol)) {
    iter = iter + 1

    f = smooth_flash_update_single_lf(data,
                               f,
                               k,
                               var_type,
                               ebnm_fn_l,
                               ebnm_param_l,
                               filter.number,
                               family,
                               Rk,
                               R2)

    R2 = (R2k - 2 * outer(f$EL[, k], f$EF[, k]) * Rk
          + outer(f$EL2[, k], f$EF2[, k]))

    obj_track[iter] = (sum(unlist(f$KL_l)) + sum(unlist(f$KL_f)) +
                         e_loglik_from_R2_and_tau(R2, f$tau, data))
    obj_diff = calc_obj_diff(obj_track, iter)

    diff = calc_diff("objective", obj_diff)
    #diff = abs(calc_diff("objective", obj_diff))
    diff_track[iter] = diff

  }

  history = list(type = "greedy",
                 kset = k,
                 niter = iter,
                 obj_track = obj_track[1:iter],
                 diff_track = diff_track[1:iter])

  return(list(f = f, history = history))

}


smooth_flash_update_single_lf = function(data,f,k,
                                  var_type,
                                  ebnm_fn_l,
                                  ebnm_param_l,
                                  filter.number,
                                  family,
                                  Rk = NULL,
                                  R2){

  #browser()

  if (is.null(Rk)) {
    Rk = flash_get_Rk(data, f, k)
  }

  f = smooth_flash_update_single_factor(data,
                                        f,
                                        k,
                                        filter.number,
                                        family,
                                        Rk,
                                        calc_obj = TRUE)

  f = flash_update_single_loading(data,
                                  f,
                                  k,
                                  ebnm_fn_l,
                                  ebnm_param_l,
                                  Rk,
                                  calc_obj = TRUE)


  return(f)

}



smooth_flash_update_single_factor = function(data,f,k,filter.number,family,Rk,calc_obj = TRUE){

  #browser()

  subset = which(!f$fixf[, k])
  any_fixed = any(f$fixf[, k])
  ebnm_args = calc_ebnm_f_args(data, f, k, subset, any_fixed, Rk)
  if (is.null(ebnm_args)) {
    return(NULL)
  }

  a = smash_wave(ebnm_args$x, ebnm_args$s,filter.number,family)
  #a = do.call(ebnm_fn, list(ebnm_args$x, ebnm_args$s, ebnm_param))

  res = list(EX = a$mu.est,
             EX2 = a$mu.est.var+(a$mu.est)^2,
             g = a$g_hat)
  if (calc_obj) {
    KL = a$loglik - NM_posterior_e_loglik(ebnm_args$x, ebnm_args$s,
                                             res$EX, res$EX2)
    res = c(res, KL = KL)
  }

  if (!is.null(res)) {

    f$EF[subset, k] = res$EX
    f$EF2[subset, k] = res$EX2
    f$gf[[k]] = res$g
    if (calc_obj) {
      f$KL_f[[k]] = res$KL
    }

  }

  return(f)

}


