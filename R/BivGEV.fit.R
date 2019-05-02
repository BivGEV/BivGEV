

################################################################################################
##############################     BivGEV.fit       ############################################
################################################################################################

BivGEV.fit <- function (func.opt, start.v, rinit, rmax, iterlim, respvec, VC)
{

sc <- TRUE

parsc <- rep(VC$parscale, length(start.v))


fit <- try(trust(func.opt, start.v, rinit = rinit, rmax = rmax,
                 parscale = parsc, respvec = respvec, VC = VC,
                 blather = TRUE, iterlim = iterlim),
                 silent = sc)

if (class(fit) == "try-error" || is.null(fit$l)) {
  fit <- try(trust(func.opt, start.v, rinit = rinit, rmax = rmax,
                   parscale = parsc, respvec = respvec, VC = VC,
                   blather = TRUE, iterlim = iterlim/4),
                   silent = sc)
  
  if (class(fit) == "try-error" || is.null(fit$l)) {
    fit <- try(trust(func.opt, start.v, rinit = rinit, rmax = rmax,
                     parscale = parsc, respvec = respvec, VC = VC,
                     blather = TRUE, iterlim = iterlim/10),
                     silent = sc)
  }
}

iter.if    <- fit$iterations
iter.inner <- bs.mgfit <- post.f <- magpp <- NULL

post.f <- work.fit(fit) 


bs.mgfit <- magic(post.f$Z, post.f$X, numeric(0), list(), numeric(0))
magpp    <- magic.post.proc(post.f$X, bs.mgfit)

if(VC$gc.m == TRUE) 
  gc()


list(fit = fit, iter.if = iter.if, 
     bs.mgfit = bs.mgfit, post.f = post.f, magpp = magpp)
}



