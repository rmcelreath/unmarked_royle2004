// partial pooling version
// this one adds correlation between detection p and lambda rate
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
    vector[2] bmeans;
    matrix[2,M] zsite;
    vector<lower=0>[2] sigma_site;
    cholesky_factor_corr[2] L_Rho_site;
}
transformed parameters{
    vector<lower=0>[M] lambda;
    vector<lower=0,upper=1>[M] p;
    {
        matrix[M,2] bsite;
        bsite = (diag_pre_multiply(sigma_site,L_Rho_site) * zsite)';
        for ( m in 1:M ) {
            lambda[m] = exp(bmeans[1] + bsite[m,1]);
            p[m] = inv_logit(bmeans[2] + bsite[m,2]);
        }
    }
}
model{
    bmeans[1] ~ normal(4,10);
    bmeans[2] ~ normal(0,4);
    to_vector(zsite) ~ normal(0,1);
    L_Rho_site ~ lkj_corr_cholesky(4);
    sigma_site ~ exponential(1);

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
