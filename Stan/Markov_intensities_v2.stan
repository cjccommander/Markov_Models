//////////////////////////////////////////////////////////
//
// continuous-time Markov Transition Model estimator
// authors: DK Okamoto, CJC Commander, T Tinker
// last update: 26 May, 2021
// 
//////////////////////////////////////////////////////////
//
functions{
  // dirichlet multinomial 
  real dirichlet_multinomial_lpmf(array[] y, vector alpha) {
    real alpha_plus = sum(alpha);
  
    return lgamma(alpha_plus) + sum(lgamma(alpha + to_vector(y))) 
      - lgamma(alpha_plus + sum(y)) - sum(lgamma(alpha));
  }
}
//
data{
  int<lower=1> n_obs;
  int<lower=1> n_states;
  int<lower=1> n_pars;
  int<lower=0> dm;
  array[n_obs] real<lower=0> time_diff;
  array[n_obs, n_states] int y_obs;
  array[n_states,n_states-1] int indxs ; // 1:nstates - diag
  matrix[n_obs, n_pars] X;
}
//
parameters{
  // log-scale transition intensities for each parameter
  array[n_states, n_states - 1] vector[n_pars] beta; 
  // scale parameter if dirichlet-multinomial is applied
  vector<lower= 0>[dm] alpha; 
}
//
transformed parameters{
  array[n_obs] vector[n_states] y_probs;
  array[n_states, n_states] vector[n_pars] Beta; 
  profile("trans pars"){
    for (i in 1 : n_states){
      for (j in indxs[i]){
        Beta[i,j] = j < i ? beta[i, j] * 3 - 2 - log(n_states-1) : 
                            beta[i, j - 1] * 3 - 2 - log(n_states-1) ;
      }
    }
    for (t in 1 : n_obs){
      if (time_diff[t] > 0){
        matrix[n_states, n_states] q_mat = rep_matrix(0,n_states, n_states) ;
        for (i in 1 : n_states){
          for (j in indxs[i]){
            q_mat[i, j] =  exp(X[t,] * Beta[i,j] ) ;
          } 
          q_mat[i, i] = - sum(q_mat[i]);
        }
        y_probs[t] = (to_row_vector(y_obs[t - 1]) * matrix_exp(time_diff[t] * q_mat))';
        // Ensure sum-to-1 for y_probs (sometimes this can be off during warmup)
        y_probs[t] = (y_probs[t] + 0.000000001) ./ (sum(y_probs[t]) + n_states * 0.000000001) ;
      }
    }
  }
}
//
model{
  profile("prior"){
    // folded unit normal, prior for intensities
    for(i in 1 : n_states){
      for(j in 1 : n_states - 1){
        target += std_normal_lpdf(beta[i, j]);
      }
    }
  }
  // folded student-t prior for dirichlet-scale
  target += student_t_lpdf(alpha| 3, 0, 2.5);
  // likelihood
  profile("likelihood"){
    for (t in 1 : n_obs){
      if (time_diff[t] > 0){
        // apply dirichlet multinomial if overdispersed
        if (dm == 1){
          target += dirichlet_multinomial_lpmf(y_obs[t] | y_probs[t] * alpha[1]);
        } else {
          target += multinomial_lpmf(y_obs[t] | y_probs[t]);
        }
      } 
    }
  }
}
