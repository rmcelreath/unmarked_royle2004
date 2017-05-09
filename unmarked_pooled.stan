// partial pooling version
// pools among population for estimating detection rates and population sizes
data{
    int n; // num surveys
    int M; // num populations
    int Y[n];
    int pop_id[n];
    int Nmax; // limit on likelihood approximation
}
transformed data{
    int Ymax[M];
    // get maximum observed Y in each population
    for ( i in 1:M ) Ymax[i] = 0;
    for ( i in 1:n ) if ( Y[i] > Ymax[pop_id[i]] ) Ymax[pop_id[i]] = Y[i];
    for ( i in 1:M ) print(Ymax[i]);
}
parameters{
    vector<lower=0>[M] lambda;
    vector<lower=0,upper=1>[M] p;
    real<lower=0> lambda_mean;
    real<lower=0,upper=1> p_mean;
    real<lower=0> p_theta;
}
model{
    lambda ~ exponential( 1.0/lambda_mean );
    lambda_mean ~ normal(50,20) T[0,];
    p ~ beta( p_mean*p_theta , (1-p_mean)*p_theta );
    p_mean ~ beta(1,1);
    p_theta ~ exponential(1);

    for ( m in 1:M ) { // loop over populations
        vector[Nmax-Ymax[m]+1] terms;
        int l;
        l = 0;
        // marginalize over pop size Ni
        for ( Ni in Ymax[m]:Nmax ) {
            // factor for probability size Ni
            l = l + 1;
            terms[l] = poisson_lpmf( Ni | lambda[m] );
            // now for each observation, binomial prob observed
            for ( i in 1:n ) {
                if ( pop_id[i]==m )
                    terms[l] = terms[l] + binomial_lpmf( Y[i] | Ni , p[m] );
            }//i
        }//Ni
        target += log_sum_exp( terms );
    }//m
}
generated quantities{
    vector[M] N_est;
    for ( m in 1:M ) { // loop over populations
        vector[Nmax-Ymax[m]+1] terms;
        real Z;
        int l;
        l = 0;
        // marginalize over pop size Ni
        for ( Ni in Ymax[m]:Nmax ) {
            // factor for probability size Ni
            l = l + 1;
            terms[l] = poisson_lpmf( Ni | lambda[m] );
            // now for each observation, binomial prob observed
            for ( i in 1:n ) {
                if ( pop_id[i]==m )
                    terms[l] = terms[l] + binomial_lpmf( Y[i] | Ni , p[m] );
            }//i
        }//Ni
        Z = log_sum_exp(terms);
        terms = exp(terms - Z);
        N_est[m] = categorical_rng( terms ) + Ymax[m] - 1;
    }//m
}
