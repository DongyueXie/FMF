library(funflash)
datax = readRDS('/project2/mstephens/dongyue/gtex/V8/data/gtex_v8.rds')
# Y_tilde = biwhitening(datax$counts)
# rm(datax)
# gc()
# fit_sf = scaledflash(Y_tilde$Y,Y_tilde$u,Y_tilde$v,
#                      S2 = NULL,
#                      var.type = 'by_column',
#                      Kmax=20,
#                      tol=0.01,
#                      maxiter = 1000,
#                      ebnm_fn = 'ebnm_pe',
#                      init_fn = 'nnmf_r1',
#                      ebnm_param=NULL,
#                      verbose=TRUE,
#                      nullcheck=TRUE,
#                      sigma2 = NULL,
#                      seed=12345)
# saveRDS(fit_sf,'biwhite_ebnmf_fit_K20.rds')


bw = biwhitening(datax$counts)
S = (1/tcrossprod(bw$u,bw$v))
library(flashier)
fit_sf <- flash.init(datax$counts,var.type = 2,S=S) %>%
  flash.add.greedy(Kmax=20,
                   maxiter=500,
                   ebnm.fn = c(ebnm::ebnm_point_exponential, ebnm::ebnm_point_exponential),
                   init.fn = function(f) init.fn.default(f, dim.signs = c(1, 1))
  )%>%
  flash.backfit() %>%
  flash.nullcheck()
saveRDS(fit_sf,'biwhite_ebnmf_flashier_K20.rds')



library(fastTopics)
datax = readRDS('/project2/mstephens/dongyue/gtex/V8/data/gtex_v8.rds')
fit = fit_topic_model(datax$counts,k=20,verbose="detailed")
saveRDS(fit,'topic_model_fit.rds')
