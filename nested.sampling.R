### R script implementing a vanilla nested sampling with random walk MCMC kernel

library(MASS)

logsumexp <- function(a,b) {
    C <- max(c(a,b))
    x <- c(a,b) - C
    return(log(sum(exp(x))) + C)
}

## Requires:
# (1) [function] log.likelihood.fn : input (parameter:vector) : output (log likelihood:numeric)
# example:
 log.likelihood.fn <- function(theta) {return(dnorm(theta,1,0.05,log=T))}

# (2) [function] log.prior.density: input (parameter:vector) : output (log prior density:numeric)
# example:
 log.prior.density <- function(theta) {return(dnorm(theta,0,1,log=T))}

# (3) [function] sample.prior: input (void) : output (parameter:vector)
# example:
 sample.prior <- function() {return(rnorm(1,0,1))}

# (4) [integer] N: number of live particles to use in nested sampling algorithm
# example:
 N <- 50

# (5) [integer] Nrw: number of random walk MCMC steps to take in attempting particle refreshment
# example:
 Nrw <- 50

# (6) [function] set.rw.covariance: input (live_particles:list_of_vectors) : output (covariance_matrix:matrix)
# example:
set.rw.covariance <- function (live_particles) {return(cov(do.call(rbind,live_particles))/length(live_particles))}

# (7) [integer] maxNfact: maximum number of nested sampling draws in multiples of N 
# example:
 maxNfact <- 7

## Returns:


## Algorithm:
nested.sampling <- function(log.likelihood.fn,log.prior.density,sample.prior,N,Nrw,set.rw.covariance,maxNfact) {
    
# Initialize live particle population
    live.particles <- list()
    live.particle.log.likelihoods <- list()
    for (i in 1:N) {
        live.particles[[i]] <- sample.prior()
        live.particle.log.likelihoods[[i]] <- log.likelihood.fn(live.particles[[i]])
    }
    logX <- logL <- discarded.particles <- list()
    
# Perform discards of the current lowest likelihood point followed by RW MCMC-based refreshment until stopping condition reached
    stopping.condition <- F
    iter <- 1
    while (!stopping.condition) {

        logL[[iter]] <- min(as.numeric(live.particle.log.likelihoods))
        logX[[iter]] <- -iter/N
        
        discard.index <- which(as.numeric(live.particle.log.likelihoods)==logL[iter])[1]
        discarded.particles[[iter]] <- live.particles[[discard.index]]

        # RW MCMC for refreshment
        starting.index <- sample(which(!(1:N %in% discard.index)),1)
        current.position <- live.particles[[starting.index]]
        current.log.prior.density <- log.prior.density(current.position)
        current.log.likelihood <- live.particle.log.likelihoods[[starting.index]]
            
        rw.cov <- set.rw.covariance(live.particles)

        for (i in 1:Nrw) {
            proposal.position <- mvrnorm(1,rep(0,length(current.position)),rw.cov)
            proposal.log.prior.density <- log.prior.density(proposal.position)
            proposal.log.likelihood <- log.likelihood.fn(proposal.position)
            if (proposal.log.likelihood > logL[iter]) {
                if (log(runif(1)) < proposal.log.prior.density-current.log.prior.density) {
                    current.position <- proposal.position
                    current.log.prior.density <- proposal.log.prior.density
                    current.log.likelihood <- proposal.log.likelihood
                }
                
            }
        }
        live.particles[[discard.index]] <- current.position
        live.particle.log.likelihoods[[discard.index]] <- current.log.likelihood

        # Check stopping condition
        current.logZest <- -Inf
        auglogX <- c(0,as.numeric(logX))
        for (i in 1:length(logX)) {
            current.logZest <- logsumexp(current.logZest,logL[[i]]+log(exp(auglogX[i])-exp(auglogX[i+1])))
        }
        remaining.logZest <- max(as.numeric(live.particle.log.likelihoods))+logX[[iter]]

        cat(current.logZest,remaining.logZest,max(as.numeric(live.particle.log.likelihoods)),logX[[iter]],"\n")
        
        if (remaining.logZest-current.logZest < log(0.1) | iter/N < maxNfact) {
            stopping.condition <- T
        }  else {
            iter <- iter + 1
        }
        
    }

# Compile outputs
    output <- list()
    log.mean.remainder <- -Inf
    for (i in 1:N) {log.mean.remainder <- logsumexp(log.mean.remainder,live.particle.log.likelihoods[[i]])}
    log.mean.remainder <- log.mean.remainder - log(N)
    output$logZ <- logsumexp(current.logZest,log.mean.remainder-(iter+1)/N)

    output$logX <- logX
    output$logL <- logL
    
    output$particles <- c(discarded.particles,live.particles)
    output$log.particle.posterior.weights <- c(log(1-exp(logX[[1]])),log(-diff(exp(as.numeric(logX)))),-(iter+1)/N-log(N))

    output$iter <- iter
    
    return(output)
    
}
