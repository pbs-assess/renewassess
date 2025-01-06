
## Use Harck Data:
library(tmbstan)
library(RTMB)

load("../data/harck")

lamw <- ADjoint(
  function(x){gsl::lambert_W0(x)},
  function(x, y, dy) {dy / (x + exp(y))}
)

## Toy function to do some different types of reporting (matrix, vector, scalar)
fn <- function(par){
  getAll(par, harck)

  negll <- 0
  ## All priors are uniform...  
  Expected_logRS <- log(alpha) - beta*S
  negll <- negll - sum(dnorm(logRS, 
                        mean = Expected_logRS, sd = sigma, log = TRUE))

  Smsy <- (1 - lamw(exp(1 - log(alpha))))/beta
  umsy <- (1 - lamw(exp(1 - log(alpha))))
  Sgen <-  -1/beta*lamw(-beta*Smsy/alpha)  
  ## Make some fake reporting tables:
  benchmarks = c(Smsy, umsy, Sgen)
  
  ## Now a matrix for predicting other years:
  pred <-  matrix(0, 2, 2)
  pred[1,1] <- Expected_logRS[1]
  pred[1,2] <- Expected_logRS[2]
  pred[2,1] <- Expected_logRS[3]
  pred[1,2] <- Expected_logRS[4]

  dimcheck <- cbind(c(1,2,3), c(4, 5, 6))
  dimcheck2 <- array(1:8, c(2,2,2))
  
  REPORT(benchmarks)
  REPORT(pred)
  logalpha = log(alpha)
  REPORT(logalpha)
  REPORT(dimcheck)
  REPORT(dimcheck2)
  
  return(negll)
}

obj <- MakeADFun(fn, list(beta=0.001,sigma=1, alpha=12.5), silent = TRUE)
fit <- nlminb(obj$par, obj$fn, obj$gr, lower = c(0,0,0), upper = c(Inf, Inf, Inf))

fit.stan <- tmbstan(obj, chains=4, iter = 1000, 
  init = c(0,0,0), lower = c(0, 0, 0), upper = c(Inf, Inf, Inf))

## Report posterior function for tmbstan:
make_posterior_report <- function(fit, obj){
  ## Let's get posterior samples for all our REPORTS.
  post <- as.matrix(fit)
  post.report <- NULL
  for( i in 1:nrow(post) ) {
    report_i <- obj$report(post[i,])
      post.report <- rbind(post.report, unlist(report_i))
  }

  make_mcmc_names <- function(x){
    outernames <- names(x)
    names <- c()
    for( i in 1:length(x) ){
      if(is.null(dim(x[[i]]))){
        if(length(x[[i]]) == 1) 
          names <- c(names, outernames[i])
        else 
          names <- c(names, paste0(outernames[i], "[", 1:length(x[[i]]), "]"))
      }else{
        dims <- dim(x[[i]])
        ndims <- length(dims)
        j <- 3
        indices <- apply(expand.grid(1:dims[1], 1:dims[2]), 1, FUN = function(x) paste(x, collapse=","))
        while(j <= ndims){
          indices <- apply(expand.grid(indices, 1:dims[j]), 1, FUN = function(x) paste(x, collapse=","))
          j <- j+1
        }
        names <- c(names, paste0(outernames[i], "[", indices, "]"))
      }
    }
    names
  }

  ## Now name check
  colnames(post.report) <- make_mcmc_names(report_i)
  post.report
}

preds <- make_posterior_report(fit = fit.stan, obj = obj)

