
## first understand the code of flash


#'
#' The input of flash function is a data matrix, or a data object returned from flash_set_data(). If
#' the standard errors S is known, can use flash_set_data().
#'
#' The greedy algorithm fits a rank-1 factor model, subtract it, then repeat.
#' The function performs this greedy approach is called flash_add_greedy(), which is a wrapper function
#' of flash_greedy_workhorse()
#'
#' In flash_greedy_workhorse(), The first step is to create an object, using handle_f() function,
#' it creates an object of flash_fit, called f. Then,
#' it iteratively fits rank-1 flash via function flash_r1(), then check
#' whether the factor/loading combination is effectively zero and stops the fitting procedure.
#'
#' In the function flash_r1(), the first step is to add a new pair of factor and loading, by using function add_factors_from_data().
#' The function initialize the new factor and loadings by using the preset init method. Then the factors and loadings are updated by
#' the function flash_optimize_single_fl(), which runs flash_update_single_fl() iteratively. In flash_update_single_fl(), it first compute R2, the expected
#' squared residual, for updating the precision tau, then update factor and loading, by running flash_update_single_loading() or ..._factor().
#'
#' I'll use flash_update_single_loading as an example. All updates are done by calc_update_vals(). The function calc_ebnm_l_args() prepares input (x,s)
#' for ebnm functuons. After fitting ebnm models, the function outputs El, El2, g, and objective function value.
#'

wave_flash = function(data,
                      Kmax = 100,
                      f_init = NULL,
                      var_type = c("by_row", "constant", "zero"),
                      init_fn = "udv_si",
                      tol = 0.01,
                      ebnm_fn = "ebnm_pn",
                      ebnm_param = NULL,
                      verbose = TRUE,
                      nullcheck = TRUE,
                      seed = 123,
                      greedy = TRUE,
                      backfit = FALSE){

  if (greedy) {
    f = flash_greedy_workhorse(data,
                         Kmax,
                         f_init,
                         var_type,
                         init_fn,
                         tol,
                         ebnm_fn,
                         ebnm_param,
                         verbose_output = "odn",
                         nullcheck,
                         seed)
  }

  return(f)
}


wave_flash_greedy_workhorse = function(data,
                                       Kmax = 100,
                                       f_init = NULL,
                                       var_type = c("by_column",
                                                    "by_row",
                                                    "constant",
                                                    "zero",
                                                    "kroneker"),
                                       init_fn = "udv_si",
                                       tol = 1e-2,
                                       ebnm_fn = "ebnm_pn",
                                       ebnm_param = NULL,
                                       verbose_output = "odn",
                                       nullcheck = TRUE,
                                       seed = 123,
                                       maxiter = 5000,
                                       stopping_rule = c("objective",
                                                         "loadings",
                                                         "factors",
                                                         "all_params")) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  f = handle_f(f_init, init_null_f = TRUE)
  data = handle_data(data, f)
  var_type = handle_var_type(var_type, data)
  init_fn = handle_init_fn(init_fn)
  ebnm_fn = handle_ebnm_fn(ebnm_fn)
  ebnm_param = handle_ebnm_param(ebnm_param, ebnm_fn, Kmax)
  verbose_output = unlist(strsplit(verbose_output, split=NULL))
  stopping_rule = match.arg(stopping_rule)

  history = list()

  prev_K = flash_get_k(f)
  for (k in 1:Kmax) {
    if (length(verbose_output) > 0) {
      verbose_greedy_next_fl(prev_K + k, stopping_rule, tol)
    }

    old_f = f
    res = flash_r1(data,
                   f,
                   var_type,
                   init_fn,
                   tol,
                   ebnm_fn$l,
                   ebnm_param$l[[k]],
                   ebnm_fn$f,
                   ebnm_param$f[[k]],
                   verbose_output,
                   nullcheck,
                   maxiter,
                   stopping_rule)

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



