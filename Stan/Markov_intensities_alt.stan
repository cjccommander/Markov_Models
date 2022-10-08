//////////////////////////////////////////////////////////
  //
  // continuous-time Markov Transition Model estimator
// authors: DK Okamoto, CJC Commander, T Tinker
// last update: 21 May, 2021
// 
  //////////////////////////////////////////////////////////
  
  functions{
    // dirichlet multinomial 
    real dirichlet_multinomial_lpmf(int[] y, vector alpha) {
      real alpha_plus = sum(alpha);
      
      return lgamma(alpha_plus) + sum(lgamma(alpha + to_vector(y))) 
      - lgamma(alpha_plus + sum(y)) - sum(lgamma(alpha));
    }
  }

data{
  int<lower=1> n_obs;
  int<lower=1> n_states;
  int<lower=1> n_pars;
  int<lower=0> dm;
  real<lower=0> time_diff[n_obs];
  int y_obs[n_obs, n_states];
  matrix[n_obs, n_pars] X;
}

parameters{
  // log-scale transition intensities for each parameter
  vector[n_pars] beta[n_states, n_states - 1]; 
  // scale parameter if dirichlet-multinomial is applied
  vector<lower= 0>[dm] alpha; 
}

transformed parameters{
  vector[n_states] y_probs[n_obs];
  
  //profile("trans pars"){
    for (t in 1 : n_obs){
      if (time_diff[t] > 0){
        
        matrix[n_states, n_states] q_mat;
        
        for (i in 1 : n_states){
          for (j in 1 : n_states){
            // if parameter is before the diagonal, place in order
            if (j < i){
              q_mat[i, j] =  exp(X[t,] * (beta[i, j]  * 3 - log(n_states-1)));
            } 
            // if diagonal, temporarily hold with zero
            if (j == i){
              q_mat[i, j] = 0;
            }
            // if parameter is after the diagonal, place one ahead
            if (j > i){
              q_mat[i, j] =  exp(X[t,] * (beta[i, j - 1] * 3 - log(n_states-1)));
            }
          } 
          q_mat[i, i] = - sum(q_mat[i, 1 : n_states]);
        }
        y_probs[t] = (to_row_vector(y_obs[t - 1]) * matrix_exp(time_diff[t] * q_mat))';
      }
    }
  //}
}

model{
  //profile("prior"){
    // folded unit normal, prior for intensities
    for(i in 1 : n_states){
      for(j in 1 : (n_states - 1)){
        target += std_normal_lpdf(beta[i, j]);
      }
    }
  //}
  
  // folded unit normal prior for dirichlet-scale
  target += student_t_lpdf(alpha| 3, 0, 2.5);
  
  // likelihood
  //profile("likelihood"){
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
  //}
}
