
update_nsmooth = function(x,s = 1, scale = c(1e-10,1e-5,0.01,0.1,1,10,100,1e5), point_mass=F,
                    nullweight=1000, shapel=1,weight = rep(1,length(x)),
                    g_init = NULL, fix_g = FALSE,
                    m = 2, control =  NULL, low = NULL,d=NULL){

  ebpm_exp_mixture(x,s,scale,point_mass,nullweight,weight,g_init,fix_g,m,control,low,d,shapel)
}
