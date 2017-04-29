# unmarked density analysis in Stan

# data sim
# M: number of populations
# N: population sizes
# p: detection rates per survey
# m: number of surveys

M <- 4
N <- c(10,20,50,100)
p <- c(0.5,0.5,0.5,0.5)
m <- rep(50,M)

y <- NULL
for ( i in 1:length(N) ) y <- c( y , rbinom(m[i],size=N[i],prob=p[i]) )
pop_id <- rep(1:length(N),times=m)

library(rethinking)

dat_list <- list(
    n = length(y),
    M = length(N),
    Y = y,
    pop_id = pop_id,
    Nmax = 200
)

m <- stan( file="unmarked.stan" , data=dat_list , chains=1 )

precis(m,2)

post <- extract.samples(m)

simplehist( post$N_est[,4] )
dens( post$lambda[,4] )

pairs(m)

