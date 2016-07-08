chou_synergy_helper = function( D, fu, fct ){
    # m is the slope of the median effect plot
    # measures the sigmoidicity of the dose effect curve; m=1 is hyperbolic
    fa = 1-fu
    fu[fu<0] = 0.001
    fa[fa<0] = 0.001
    y_m = rep(0, length(fu) )
    not_zero = !fa==0 & !fu==0
    y_m[ not_zero ] = log( fa[not_zero] / fu[not_zero] )
    
    x_m = rep(0, length(D))
    not_zero = D != 0
    x_m[not_zero] = log( D[not_zero] )

    m = as.numeric(coefficients( lm( y_m ~ x_m ) )[2])
    
    # Dm is the IC50 for this curve
    data = data.frame(fu, D)
    ed = drc::ED( drc::drm(fu~D, data=data, fct=fct), c(50), display=FALSE )
    Dm = ed[1]
    Dm_stderr = ed[2]
    list( D=D, fa=fa, fu=fu, y_m=y_m, x_m=x_m, m=m, Dm=Dm, Dm_stderr=Dm_stderr )
}

#' Calculate statistics needed for Chou Synergy
#' Curve fitting is performed by the \code{drm()} function in the \code{drc} 
#' library. To fit the curve, you need to select a non-linear function. To 
#' estimate the slope, upper asymptote, lower asymptote, and EC50, pass 
#' drc::LL.4(). To fix the lower asymptote at 1 and estimate the other 
#' parameters, pass drc::LL.3(). To fix the upper asympotote at 1 and the lower 
#' asymptote at 0, pass dcr::LL.2. For a list of available functions, see 
#' \code{drc::getMeanFunctions()}. 
#'
#' Ting-Chao Chou Pharmacological Reviews 2006
#'
#' @param ds dataset
#' @param sample_type sample type in ds
#' @param treatment_1 treatment in ds
#' @param treatment_2 treatment in ds
#' @param treatment_12 treatment in ds
#' @param hour hour in ds
#' @param proportion_1 value between 0 and 1, indicating the fixed proportional
#' relationship between concentrations of treatment_1 and treatment_2 in the 
#' combined observations. If 50/50, pass 0.5.
#' @param fct Non-linear function to fit, e.g. drc::LL.3(). See summary.
#' @param summary_method mean or median
#' @export
chou_synergy = function( ds, sample_type, treatment_1, treatment_2, hour, 
                         proportion_1, fct, summary_method ){
    # m is the slope of the median effect plot
    # measures the sigmoidicity of the dose effect curve; m=1 is hyperbolic
    
    if( summary_method != "mean" & summary_method != "median"){
        stop("summary_method parameter must be either mean or median")
    }
    if( sum(ds$sample_type==sample_type)==0 ){
        stop("sample_type not found in ds")   
    }
    if( sum(ds$treatment==treatment_1)==0 ){
        stop("treatment_1 not found in ds")   
    }
    if( sum(ds$treatment_2==treatment_2)==0 ){
        stop("treatment_2 not found in ds")   
    }
    if( sum( names(ds)=="treatment_2" )==0 ){
        stop(paste("dataset ds must be a synergy dataset with columns called",
            "treatment_2 and concentration_2"))   
    }
    ds_1 = ds[ds$sample_type==sample_type & ds$hours==hour & 
                  ds$treatment==treatment_1 &
                  ds$concentration_2==0 & !ds$is_negative_control,]
    ds_2 = ds[ds$sample_type==sample_type & ds$hours==hour & 
                  ds$treatment_2==treatment_2 &
                  ds$concentration==0 & !ds$is_negative_control,]
    ds_c = ds[ds$sample_type==sample_type & ds$hours==hour & 
                  ds$treatment==treatment_1 & ds$treatment_2==treatment_2 &
                  ds$concentration != 0 & ds$concentration_2 != 0 & 
                  !ds$is_negative_control,]
    ds_1$conc_final = ds_1$concentration
    ds_2$conc_final = ds_2$concentration_2
    ds_c$conc_final = ds_c$concentration + ds_c$concentration_2

    get_mu = function(x){ Fa=mean(x$value_normalized, na.rm=TRUE) }
    get_median = function(x){ Fa=median(x$value_normalized, na.rm=TRUE) }
    
    if( summary_method=="mean" ){
        Dsum1 = plyr::ddply(ds_1, c("conc_final"), get_mu )
        Dsum2 = plyr::ddply(ds_2, c("conc_final"), get_mu )
        Dsumc = plyr::ddply(ds_c, c("conc_final"), get_mu )
    }else{
        Dsum1 = plyr::ddply(ds_1, c("conc_final"), get_median )
        Dsum2 = plyr::ddply(ds_2, c("conc_final"), get_median )
        Dsumc = plyr::ddply(ds_c, c("conc_final"), get_median )
    }
    names(Dsum1) = c("D", "Fu")
    names(Dsum2) = c("D", "Fu")
    names(Dsumc) = c("D", "Fu")
    cs_1 = chou_synergy_helper( Dsum1$D, Dsum1$Fu, fct )
    cs_2 = chou_synergy_helper( Dsum2$D, Dsum2$Fu, fct )
    cs_c = chou_synergy_helper( Dsumc$D, Dsumc$Fu, fct )
    
    Fa = cs_c$fa
    Fu = 1-cs_c$fa
    proportion_2 = 1-proportion_1
    CI_1 = (cs_c$D * proportion_1) / (cs_1$Dm * (Fa/Fu)^(1/cs_1$m) )
    CI_2 = (cs_c$D * proportion_2) / (cs_2$Dm * (Fa/Fu)^(1/cs_2$m) )
    CI = data.frame( Fa=Fa, Fu=Fu, CI=CI_1+CI_2 )
    L=list( treatment_1 = cs_1, 
            treatment_2 = cs_2, 
            treatment_12 = cs_c, 
            CI=CI)
}


#' Plot single-value Chou Combination Index at IC50
#' 
#' @param CS chou statistics calculated by \code{\link{chou_synergy_helper}}
#' @param proportion proportion of dose 1 vs. dose 2, numeric between 0 and 1
#' @return Combination Index
#' @export
chou_synergy_CI_median = function( CS, proportion_1 ){
    D1 = CS$treatment_12$Dm * proportion_1
    D2 = CS$treatment_12$Dm * (1-proportion_1)
    Dm1 = CS$treatment_1$Dm
    Dm2 = CS$treatment_2$Dm
    (D1/Dm1) + (D2/Dm2)
}

#' Construct confidence intervals for observed effects at combination doses 
#' having observed effects.
#' 
#' Copied directly from code published in 
#' Lee and Kong Statistics in Biopharmaceutical Research 2012
#' 
#' @param ds dataset
#' @param sample_type sample type in ds
#' @param treatment_1 treatment in ds
#' @param treatment_2 treatment in ds
#' @param proportion_1 value between 0 and 1, indicating the fixed proportional
#' relationship between concentrations of treatment_1 and treatment_2 in the 
#' combined observations. If 50/50, pass 0.5.
#' @param E fixed effects, their corresponding interaction indices and 
#' confidence intervals are estimated.
#' @param hour hour in ds. Default 0. 
#' @param alpha 1-alpha is the size of the confidence intervals, default 0.05
#' @return list ii: estimated interaction indices corresponding to the 
#' observations (c.d1, c.d2, E); ii.low, ii.up: estimated lower and upper CI
#' @references Lee & Kong Statistics in Biopharmaceutical Research 2012
#' @export
synergy_interaction_CI = function(ds, sample_type, 
                                   treatment_1, treatment_2, proportion_1,
                                  E, hour=0, alpha=0.05,
                                  summary_method="mean"){
    
    if(! (summary_method=="mean" | summary_method=="median" ) ){
        stop("parameter summary_method must be either mean or median")   
    }
    if( sum(ds$sample_type==sample_type)==0 ){
        stop("sample_type not found in ds")   
    }
    if( sum(ds$treatment==treatment_1)==0 ){
        stop("parameter treatment_1 with passed value not found in ds")   
    }
    if( sum(ds$treatment_2==treatment_2)==0 ){
        stop("parameter treatment_2 with passed value not found in ds")   
    }
    if( sum(ds$hours==hour)==0 ){
        stop("parameter hour with passed value not found in ds")   
    }
    idx_1=ds$sample_type==sample_type & ds$treatment==treatment_1 & 
          ds$concentration>0 & ds$concentration_2==0 & ds$hours==hour
    idx_2=ds$sample_type==sample_type & ds$treatment_2==treatment_2 & 
          ds$concentration==0 & ds$concentration_2 > 0 & ds$hours==hour
    idx_12=ds$sample_type==sample_type & 
        ds$treatment == treatment_1 & ds$treatment_2==treatment_2 & 
        ds$concentration>0 & ds$concentration_2 > 0 & ds$hours==hour
    if( sum(idx_1)==0 )
        stop("No samples meet the criteria for treatment_1")
    if( sum(idx_2)==0 )
        stop("No samples meet the criteria for treatment_2")
    if( sum(idx_12)==0 )
        stop("No samples meet the criteria for combined treatments")
    if( summary_method=="mean" ){
        sumfunc = function(po){ data.frame( 
            value=mean(po$value_normalized, na.rm=TRUE)) }
    }
    else{
        sumfunc = function(po){ data.frame( 
            value=median(po$value_normalized, na.rm=TRUE)) }
    }
    E[ E > 1 ] = 0.999
    
    ds$conc_final = ds$concentration + ds$concentration_2
    E1 = plyr::ddply( ds[idx_1,], c("conc_final"), sumfunc )
    e1 = E1$value
    d1 = E1$conc_final
    
    E2 = plyr::ddply( ds[idx_2,], c("conc_final"), sumfunc)
    e2 = E2$value
    d2 = E2$conc_final
    
    E12 = plyr::ddply( ds[idx_12,], c("conc_final"), sumfunc)
    e12 = E12$value
    d12 = E12$conc_final
    
    e1[e1>1] = 0.999
    e2[e2>1] = 0.999
    e12[e12>1] = 0.999
    
    d2.d1 = proportion_1

    lm1 <- lm(log(e1/(1-e1))~log(d1))
    dm1 <- exp(-summary(lm1)$coef[1,1]/summary(lm1)$coef[2,1])
    lm2 <- lm(log(e2/(1-e2))~log(d2))
    dm2 <- exp(-summary(lm2)$coef[1,1]/summary(lm2)$coef[2,1])
    lmcomb <- lm(log(e12/(1-e12))~log (d12))
    dm12 <- exp(-summary(lmcomb)$coef[1,1]/summary(lmcomb)$coef[2,1]) 
    Dx1 <- dm1*(E/(1-E))^(1/summary(lm1)$coef[2,1])
    Dx2 <- dm2*(E/(1-E))^(1/summary(lm2)$coef[2,1])
    dx12 <- dm12*(E/(1-E))^(1/summary(lmcomb)$coef[2,1])
    iix <- (dx12/(1+d2.d1))/Dx1+(dx12*d2.d1/(1+d2.d1))/Dx2
    lm1.s <-summary(lm1)
    lm2.s <-summary(lm2)
    lm12.s <-summary(lmcomb)
    c1 <- 1.0/lm1.s$coef[2,1]^2*lm1.s$coef[1,2]^2
    temp <- - mean(log(d1))*lm1.s$coef[2,2]^2
    ### temp <- lm1.s$coef[1,2]*lm1.s$coef[2,2]*lm1.s$cor[1,2]   ### covariance of b0 and b1
    c1 <- c1+2.0*(log(E/(1-E))-lm1.s$coef[1,1])/lm1.s$coef[2,1]^3*temp
    c1 <- c1+(log(E/(1-E))-lm1.s$coef[1,1])^2/lm1.s$coef[2,1]^4*lm1.s$coef[2,2]^2
    c2 <- 1.0/lm2.s$coef[2,1]^2*lm2.s$coef[1,2]^2
    temp <- - mean(log(d2))*lm2.s$coef[2,2]^2
    ### temp <- lm2.s$coef[1,2]*lm2.s$coef[2,2]*lm2.s$cor[1,2]   ### covariance of b0 and b1
    c2 <- c2+2.0*(log(E/(1-E))-lm2.s$coef[1,1])/lm2.s$coef[2,1]^3*temp
    c2 <- c2+(log(E/(1-E))-lm2.s$coef[1,1])^2/lm2.s$coef[2,1]^4*lm2.s$coef[2,2]^2
    c12 <- 1.0/lm12.s$coef[2,1]^2*lm12.s$coef[1,2]^2
    temp <- - mean(log(d12))*lm12.s$coef[2,2]^2
    ### temp <- lm12.s$coef[1,2]*lm12.s$coef[2,2]*lm12.s$cor[1,2]   ### covariance of b0 and b1
    c12 <- c12+2.0*(log(E/(1-E))-lm12.s$coef[1,1])/lm12.s$coef[2,1]^3*temp
    c12 <- c12+(log(E/(1-E))-lm12.s$coef[1,1])^2/lm12.s$coef[2,1]^4*lm12.s$coef[2,2]^2
    var.ii <-((dx12/Dx1)^2*c1+(dx12*d2.d1/Dx2)^2*c2+(1.0/Dx1+d2.d1/Dx2)^2*dx12^2*c12)/(1+d2.d1)^2 
    t975 <- qt(1-alpha/2,length(d1)+length(d2)+length(d12)-6)
    iix.low1 <- iix*exp(-t975*var.ii^0.5/iix)
    iix.up1 <- iix*exp(t975*var.ii^0.5/iix)
    
    iix.low1[ is.nan(iix.low1) ] = NA
    iix.up1[ is.nan(iix.up1) ] = NA
    iix.up1[ is.nan(iix.up1) ] = NA
    return(list(interaction_index=iix, cl_lower=iix.low1, cl_upper=iix.up1))

}