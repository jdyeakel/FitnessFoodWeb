#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List SDP_single_foreq(int tmax, NumericVector res_bs, int cons_bs, int xc, NumericVector rep_gain, 
NumericMatrix f_m, NumericVector mort, List dec_ls, NumericVector rho_vec, NumericMatrix c_learn, 
NumericVector g_forage, NumericVector c_forage, List W_nr, List istar_nr, int N_init, int tsim, 
double alpha, double beta, double comp) {
   
   //Note:
   //Here, state variables are in standard format
   double xc_state = (double) xc; //state
   
   //To use a state variable X as an INDEX, it must be an INTEGER, and X=X-1
   xc = xc - 1; //Index
   
   //Book-keeping
   //Establish iniital variables required for the SDP
   int num_res = res_bs.size();
   //How many decisions? Count columns on first decision matrix
   NumericMatrix dec_ex = dec_ls(0);
   //int num_dec = dec_ex.ncol();
   int max_enc = f_m.nrow(); //Rows = encounters; Columns = resources
   
   //This is a state, not an index, so should equal cons_bs
   double xmax = (double) cons_bs;
   
   //Vectors for output
   IntegerVector pop_size(tsim);
   //OTHERS
   
   //Number of starting individuals
   int N = N_init;
   
   //Initialize the record-keeping of states
   //Record the state vector for each individual
   
   //Which of these individuals are alive (they all start out alive)
   IntegerVector alive(N);
   for (int n=0; n<N; n++) {
     alive(n) = 1;
   }
   int num_alive = sum(alive); //so num_alive = N at first
   
   //This records the energetic state of each individual... values will change over time
   IntegerVector state_v(N);
   for (int i=0; i<N; i++) {
     state_v(i) = xmax; //The initial state is full health
   }
   
   //Record the time vector for each individual
   IntegerVector time_v(N);
   for (int i=0; i<N; i++) {
     time_v(i) = 1; //tmax-1;
   }
   
   //Record the current resource for each individual
   IntegerVector res_v(N);
   for (int i=0; i<N; i++) {
    NumericVector rdraw_v = runif(1);
    double rdraw = as<double>(rdraw_v);
    //The initial resource is randomly drawn for each individual
    //Random draw is between 0 and num_res-1... so already 'indexed'
    res_v(i) = (int) floor(rdraw*(num_res)); 
   }
   
   
   //What are the optimal decisions for the current individual-level states?
   IntegerVector dec_v(N);
   for (int i=0; i<N; i++) {
     //res_v is already indexed
     //state_v and time_v are states, not indices, so alter by -1
     IntegerMatrix istar_t = istar_nr(res_v(i)); //Grab the istar for the focal resource
     dec_v(i) = istar_t(state_v(i)-1,time_v(i)-1); //Grab the decision for a given state and time t=0;
   }
   
   //TO EXPORT
   //Record initial pop size
   pop_size(0) = num_alive;
   
   List trophic_int(tsim);
   List energetic_state(tsim);
   List temporal_state(tsim);
   List decision_state(tsim);
   List fitness(tsim);
   
   trophic_int(0) = res_v;
   energetic_state(0) = state_v;
   temporal_state(0) = time_v;
   decision_state(0) = dec_v;
   
   NumericVector ind_fit(N);
   for (int i=0; i<N; i++) {
     NumericMatrix W_matrix = W_nr(res_v(i));
     ind_fit(i) = W_matrix(state_v(i)-1,time_v(i)-1);
   }
   fitness(0) = ind_fit;
   
   //Begin time iteration (this is the simulation time)
   //Sim time iterations stop at tsim - 1
   for (int t=0; t<(tsim-1); t++) {
     
     //Rcpp::Rcout << "tsim = " << t << std::endl;
     
     //Book-keeping
     //Reset metrics (first iteration will be repetitive, but oh well)
     //How many individuals are currently alive?
     int N = num_alive;
     
     //Create vector to record reproductive events (1,0)
     IntegerVector rep_success(N);
     
     //To save the spawning stock biomass
     IntegerVector SSB_v(N);
     
     //To save the 'next' values for the energetic vector
     NumericVector x_next_rd(N);
     
     //To save the 'next' resouce consumed
     IntegerVector res_next(N);
     
     //To save the 'next' values for the time vector
     //IntegerVector t_next(N);
     
     //Individual fitness values
     NumericVector ind_fit(N);
     
     //Begin Individual iterations... Note: N will change with each t
     for (int n=0; n<N; n++) {
       
       
       
       //Rcpp::Rcout << "N = " << n << std::endl;
       
       //Define states for the individual
       //Current energetic state (state)
       int ind_x = state_v(n);
       //Current temporal state (state)
       //int ind_t = time_v(n);
       //Current (focal) resource (index)
       int ind_r = res_v(n);
       //Decision matrix associated with current resource (state)
       int ind_d = dec_v(n);
       
       //Reporting
       //Rcpp::Rcout << "n = " << n << std::endl;
       //Rcpp::Rcout << "ind_r = " << ind_r << std::endl;
       //Rcpp::Rcout << "ind_d = " << ind_d << std::endl;
       //Rcpp::Rcout << "ind_x = " << ind_x << std::endl;
       
       //REPRODUCTION
       //What is the probability that the individual will reproduce?
       double rep_prob = rep_gain(ind_x-1);
       
       
       //Does the individual reproduce in this timestep?
       //Take care of reproduction after n iterations... to include density dependence
       NumericVector rdraw_v = runif(1);
       double rdraw = as<double>(rdraw_v);
       if (rdraw < rep_prob) {
         rep_success(n) = 1;
         SSB_v(n) = ind_x;
       } else {
         rep_success(n) = 0;
         SSB_v(n) = 0;
       }
       
       //Longevity mortality + Stochastic mortality
       //If the individual has not reached it's max longevity
       int stochmort;
       if (time_v(n) < tmax) {
         double pr_stochmort = mort(ind_r); //ind_r is index already
         NumericVector rdraw_stochmort_v = runif(1);
         double rdraw_stochmort = as<double>(rdraw_stochmort_v);
         if (rdraw_stochmort < pr_stochmort) {
           //Individual dies
           stochmort = 1;
         } else {
           //Individual survives
           stochmort = 0;
         }
       } else {
         stochmort = 1; //All individuals that reach tmax die
       }
       
       //The following steps are only relevant for surviving individuals
       
       if (stochmort == 0) {
         
         //Choose next resource as a function of the preference probability distribution
         NumericMatrix pref_prob_m = dec_ls(ind_r); //ind_r is index already
         NumericVector pref_prob(num_res);
         IntegerVector pref_num(num_res);
         IntegerVector num_prob(num_res);
         for (int i=0;i<num_res; i++) {
           pref_prob(i) = pref_prob_m(i,ind_d-1);
           //The number of each 'resource index' to choose across
           //Zeros will screw this up
           pref_num(i) = (int) floor(pref_prob(i)*1000) + 1; //+1 to ensure no zeros
         } 
         int tot_num = sum(pref_num);
         //The p_bin starts at 0 and ends at num_res-1... so it is already an index
         IntegerVector p_bin(tot_num);
         int tic = 0;
         for (int i=0; i<num_res; i++) {
           int num = pref_num(i);
           for (int j=0; j<num; j++) {
             p_bin(tic) = i;
             tic = tic + 1;
           }
         }
         NumericVector rdraw_v = runif(1);
         double rdraw = as<double>(rdraw_v);
         //Random draw between 0 and tot_num...
         int draw = (int) floor(rdraw*(tot_num));
         //This defines the next resource (the index for the next resource)
         //NOTE: this is an index, as it has been for res
         res_next(n) = p_bin(draw);
         
         
         
         //How many of this resource is found?
         NumericVector k_prob(max_enc+1);
         IntegerVector k_num(max_enc+1);
         for (int i=0; i<(max_enc+1); i++) { //there are max_enc+1 indices
           k_prob(i) = f_m(i,res_next(n)); //res_next is already an index
           //Zeros will screw this up
           k_num(i) = (int) floor(k_prob(i)*1000) + 1; //+1 to ensure no zeros
         }
         tot_num = sum(k_num);
         IntegerVector k_bin(tot_num);
         tic = 0;
         //These k values are not 'k', but the k index
         //This index will grab the appropriate value in the rho vector
         for (int i=0; i<(max_enc+1); i++) {
           int num = k_num(i);
           for (int j=0; j<num; j++) {
             k_bin(tic) = i;
             tic = tic + 1;
           }
         }
         rdraw_v = runif(1);
         rdraw = as<double>(rdraw_v);
         draw = (int) floor(rdraw*tot_num);
         //This defines the index for the resource amt / rho
         int k = k_bin(draw);
         
        
         
         //What is the probability of catching the resource: rho_vec(k)
         //Rcpp::Rcout << "k= " << k << std::endl;
         //Must be k-1 because k is a state, and needs to be used as an index
         double rho = rho_vec(k-1);
         
          
         //Is the prey capture successfull?
         //Energetic dynamics
         rdraw_v = runif(1);
         double x_next; //state
         double d_ind_x = (double) ind_x; //double current state
         rdraw = as<double>(rdraw_v);
         
         if (rdraw < rho) {
           //If food is captured (res_next is an index)
           x_next = d_ind_x + g_forage(res_next(n)) - c_forage(res_next(n));
         } else {
           //If food is not captured
           x_next = d_ind_x - c_forage(res_next(n));
         }
         
         //Turn state x_next into nearest integer
         x_next_rd(n) = (double) round(x_next);
         
         //Boundary conditions and determine whether individual dies
         if (x_next_rd(n) < xc_state) {
           x_next_rd(n) = 0;
           alive(n) = 0;
         }
         //Rcpp::Rcout << "Made it here; t = " << t << std::endl;
         if (x_next_rd(n) > xmax) {
           x_next_rd(n) = xmax;
           alive(n) = 1;
         }
         
         //If individual suffers stochastic mortality
       } else { 
         
         x_next_rd(n) = 0;
         alive(n) = 0;
         
       } //end stochmort if statement
       
       //Record individual fitness value
       NumericMatrix W_matrix = W_nr(res_next(n));
       ind_fit(n) = W_matrix(x_next_rd(n)-1,time_v(n)-1);
       
       
       //Next individual time integer
       //Shouldn't matter if we advance time for dead individuals. Culling occurs later
       //t_next(n) = time_v(n) + 1;
      
     } //end individual iterations
     
     
     
     //Import metrics
     //Distribution of individuals consuming each resource at time t
     trophic_int(t+1) = res_v;
     //Distribution of individual energetic states over time (distribution)
     energetic_state(t+1) = state_v;
     //Distribution of individuals at each ind_timestep at time t (distribution)
     temporal_state(t+1) = time_v;
     //Distribution of decisions at time t
     decision_state(t+1) = dec_v;
     //Distribution of fitness values at time t
     fitness(t+1) = ind_fit;
     //Imports
     
     num_alive = sum(alive);
     
     
     //Density-dependent Reproduction :: TODO
     //Must add on new individuals to the current state vectors
     
     //This is the sum of biomass that is reproducing in this timestep
     //Regardless of mortality
     double SSB = sum(SSB_v); 
     double recruitB = (alpha*SSB) / (1 + beta*pow(SSB,1/comp));
     
     
     //The number of new individuals: Recruit biomass / mass of individual
     int num_recruits = (int) round(recruitB/xmax);
     //int num_recruits = 0;

     
     //Total individual count = surviving individuals + recruits
     int new_N = num_alive + num_recruits;
     
     //Rcpp::Rcout << "num_alive = " << num_alive << std::endl;
     //Rcpp::Rcout << "recruit = " << num_recruits << std::endl;
     //Rcpp::Rcout << "Total = " << new_N << std::endl;
     
     //Add on new individuals
     IntegerVector new_alive(new_N);
     IntegerVector new_state_v(new_N);
     IntegerVector new_res_v(new_N);
     IntegerVector new_dec_v(new_N);
     IntegerVector new_time_v(new_N);
     

     //Transcribe over values for states at NEXT time interval
     if (num_alive > 0) {
       int tic = 0;
       for (int i=0; i<N; i++) {
         //only do this for surviving individuals
         if (alive(i) == 1) {
           new_alive(tic) = alive(i);
           new_state_v(tic) = x_next_rd(i); //state
           new_res_v(tic) = res_next(i); //index
           new_time_v(tic) = time_v(i) + 1; //advance time interval
           IntegerMatrix istar_t = istar_nr(new_res_v(tic)); //Grab the istar for the focal resource
           new_dec_v(tic) = istar_t(new_state_v(tic)-1,new_time_v(tic)-1);
           tic = tic + 1;
         }
       }
     }
     
     
     //Initiate new values for recruits
     if (num_recruits > 0) {
       int tic = 0; //just for the rdraw_v index
       for (int i=num_alive; i<new_N; i++) {
         new_alive(i) = 1; //Recruits are all alive
         new_state_v(i) = xmax; //Initial energetic state is xmax
         new_time_v(i) = 1;
         //Randomly draw initial resource
         NumericVector rdraw_v = runif(1);
         double rdraw = as<double>(rdraw_v); 
         new_res_v(i) = (int) floor(rdraw*(num_res));
         IntegerMatrix istar_t = istar_nr(new_res_v(i)); //Grab the istar for the focal resource
         new_dec_v(i) = istar_t(new_state_v(i)-1,new_time_v(i)-1); //Grab the decision for a given state and time t=0;
         tic = tic + 1;
       }
     }

     //RESET
     //How many individuals are there?
     num_alive = sum(new_alive);
     
     //Redefine states only for alive individuals
     state_v = new_state_v;
     time_v = new_time_v;
     res_v = new_res_v;
     dec_v = new_dec_v;
     alive = new_alive;
     
     //EXPORT
     //Save variables
     pop_size(t+1) = num_alive;
     
     
     
     
   } //end simulation time iterations
   
   List output(6);
   output(0) = pop_size;
   output(1) = energetic_state;
   output(2) = temporal_state;
   output(3) = decision_state;
   output(4) = trophic_int;
   output(5) = fitness;
   //output(1) = ;
   
   return output;
   
}
