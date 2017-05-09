# unmarked density analysis in Stan

# data sim
# M: number of populations
# N: population sizes
# p: detection rates per survey
# m: number of surveys

M <- 10
N <- rep_len(c(10,20,50,75,100),length.out=M)
N <- sort(N)
p <- rep( 0.5 , M )
p <- runif( M , 0.2 , 0.8 )
m <- rep( 20 , M )

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

mp <- stan( file="unmarked_pooled.stan" , data=dat_list , chains=1 )

mp2 <- stan( file="unmarked_pooled2.stan" , data=dat_list , chains=1 )

precis(m,2,spark=10,pars="N_est")
precis(mp,2,spark=10,pars="N_est")

# show unpooled vs pooled N estimates

p1 <- extract.samples(m)
p2 <- extract.samples(mp)

med1 <- apply(p1$N_est,2,mean)
med2 <- apply(p2$N_est,2,mean)
ci1 <- apply(p1$N_est,2,PI,prob=0.5)
ci2 <- apply(p2$N_est,2,PI,prob=0.5)

blank()

# show true vs posterior mean
plot(N,med2,pch=16,xlim=c(0,100),ylim=c(0,max(med1,med2)) , xlab="True population size (N)" , ylab="posterior mean N")
points(N+2,med1)
for ( i in 1:M ) lines( c(N[i],N[i])+2 , ci1[,i] , lwd=0.75 )
for ( i in 1:M ) lines( c(N[i],N[i]) , ci2[,i] , lwd=1 )

abline(a=0,b=1,lwd=0.5)
points( 3 , 95 )
text( 3 , 95 , "unpooled" , cex=0.8 , pos=4 )
points( 3 , 89 , pch=16 )
text( 3 , 89 , "partially pooled" , cex=0.8 , pos=4 )  

l1 <- apply(p1$lambda,2,mean)
l2 <- apply(p2$lambda,2,mean)
pd1 <- apply(p1$p,2,mean)
pd2 <- apply(p2$p,2,mean)

plot( l1 , pd1 , ylim=c(0,1) , xlim=c(min(l1,l2),max(l1,l2)) , xlab="lambda (mean N hyperprior)" , ylab="Probability detect (p)" )
points( l2 , pd2 , pch=16 )
for ( i in 1:M ) lines( c(l1[i],l2[i]) , c(pd1[i],pd2[i]) , lwd=0.75 )
points( 19 , 0.95 )
text( 19 , .95 , "unpooled" , cex=0.8 , pos=4 )
points( 19 , .89 , pch=16 )
text( 19 , .89 , "partially pooled" , cex=0.8 , pos=4 ) 

