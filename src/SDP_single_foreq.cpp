#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List SDP_single_foreq(int tmax, NumericVector res_bs, int cons_bs, int xc, NumericVector rep_gain, 
NumericMatrix f_m, NumericVector mort, List dec_ls, NumericVector rho_vec, NumericMatrix c_learn, 
NumericVector g_forage, NumericVector c_forage, List W_nr, List istar_nr, int N_init, int tsim, double eta,
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
   int num_dec = dec_ex.ncol();
   int max_enc = f_m.nrow(); //Rows = encounters; Columns = resources
   
   //This is a state, not an index, so should equal cons_bs
   double xmax = (double) cons_bs;
   
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
     time_v(i) = 1;
   }
   
   //Record the current resource for each individual
   IntegerVector res_v(N);
   NumericVector rdraw_v = runif(N);
   for (int i=0; i<N; i++) {
    double rdraw = as<double>(rdraw_v(i));
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
   
   
   //Begin time iteration (this is the simulation time)
   for (int t=0; t<tsim; t++) {
     
     //Book-keeping
     //Reset metrics (first iteration will be repetitive, but oh well)
     //How many individuals are currently alive?
     int N = num_alive;
     
     //Create vector to record reproductive events (1,0)
     IntegerVector rep_success(N);
     
     //To save the spawning stock biomass
     IntegerVector SSB_v(N);
     
     //To save the 'next' values for the energetic vector
     IntegerVector x_next_rd(N);
     
     //To save the 'next' resouce consumed
     IntegerVector res_next(N);
     
     //To save the 'next' values for the time vector
     IntegerVector t_next(N);
     
     //Vector of stochastic mortality draws for the n iterations
     NumericVector rdraw_stochmort_v = runif(N);
     
     //Begin Individual iterations... Note: N will change with each t
     for (int n=0; n<N; n++) {
       
       //Define states for the individual
       //Current energetic state
       int ind_x = state_v(n);
       //Current temporal state
       int ind_t = time_v(n);
       //Current (focal) resource
       int ind_r = res_v(n);
       //Decision matrix associated with current resource
       int ind_d = dec_v(n);
       
       //REPRODUCTION
       //What is the probability that the individual will reproduce?
       rep_prob = rep_gain(ind_x);
       //Does the individual reproduce in this timestep?
       //Take care of reproduction after n iterations... to include density dependence
       NumericVector rdraw_v = runif(1);
       double rdraw = as<double>(rdraw_v);
       if (rdraw < rep_prob) {
         int rep_success(n) = 1;
         SSB_v(n) = ind_x;
       } else {
         int rep_success(n) = 0;
         SSB_v(n) = 0;
       }
       
       //Stochastic mortality
       double pr_stochmort = mort(ind_r);
       double rdraw_stochmort = as<double>(rdraw_stochmort_v(n));
       int stochmort;
       if (rdraw_stochmort < pr_stochmort) {
         //Individual dies
         stochmort = 1;
       } else {
         //Individual survives
         stochmort = 0;
       }
       
       //The following steps are only relevant for surviving individuals
       
       if (stochmort == 0) {
         
         //Choose next resource as a function of the preference probability distribution
         NumericMatrix pref_prob_m = dec_ls(ind_r);
         NumericVector pref_prob(num_res);
         IntegerVector num_prob(num_res);
         for (int i=0;i<num_res; i++) {
           pref_prob(i) = pref_prob_m(i,ind_d);
           pref_num(i) = floor(pref_prob(i)*1000) + 1; //+1 to ensure no zeros
         } 
         int tot_num = sum(pref_num);
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
         //Random draw between 0 and tot_num
         int draw = (int) floor(rdraw*(tot_num));
         //This defines the next resource (the index for the next resource)
         res_next(n) = p_bin(draw);
         
         //How many of this resource is found?
         NumericVector k_prob(max_enc);
         for (int i=0; i<max_enc; i++) {
           k_prob(i) = f_m(i,res_next);
           k_num(i) = floor(k_prob(i)*1000) + 1; //+1 to ensure no zeros
         }
         int tot_num = sum(k_num);
         IntegerVector k_bin(tot_num);
         int tic = 0;
         for (int i=0; i<max_enc; i++) {
           int num = k_num(i);
           for (int j=0; j<num; j++) {
             k_bin(tic) = i;
             tic = tic + 1;
           }
         }
         NumericVector rdraw_v = runif(1);
         double rdraw = as<double>(rdraw_v);
         int draw = (int) floor(rdraw*tot_num);
         //This defines the amount of the next resource
         int k = k_bin(draw);
         
         //What is the probability of catching the resource: rho_vec(k)
         double rho = rho_vec(k);
         
         //Is the prey capture successfull?
         //Energetic dynamics
         NumericVector rdraw_v = runif(1);
         int x_next;
         double rdraw = as<double>(rdraw_v);
         if (rdraw < rho) {
           //If food is captured
           x_next = ind_x + eta*g_forage(res_next) - c_forage(res_next);
         } else {
           //If food is not captured
           x_next = ind_x - c_forage(res_next);
         }
         
         //Turn x_next into nearest integer
         x_next_rd(n) = round(x_next);
         
         //Boundary conditions and determine whether individual dies
         if (x_next_rd < xc_state) {
           x_next_rd(n) = 0;
           alive(n) = 0;
         }
         if (x_next_rd > xmax) {
           x_next_rd(n) = xmax;
           alive(n) = 1;
         }
         
         //If individual suffers stochastic mortality
       } else { 
         
         x_next_rd(n) = 0;
         alive(n) = 0;
         
       } //end stochmort if statement
       
       //Next individual time integer
       //Shouldn't matter if we advance time for dead individuals. Culling occurs later
       t_next(n) = time_v(n) + 1;
       
       
     } //end individual iterations
     
     int num_alive = sum(alive);
     
     //Density-dependent Reproduction :: TODO
     //Must add on new individuals to the current state vectors
     
     //This is the sum of biomass that is reproducing in this timestep
     //Regardless of mortality
     double SSB = sum(SSB_v); 
     double recruitB = (alpha*SSB) / (1 + beta*pow(SSB,1/comp));
     
     
     //The number of new individuals: Recruit biomass / mass of individual
     int num_recruits = (int) round(recruitB/xmax);
     
     //Total individual count = surviving individuals + recruits
     new_N = num_alive + num_recruits;
     //Add on new individuals
     IntegerVector new_alive(new_N);
     IntegerVector new_state_v(new_N);
     IntegerVector new_res_v(new_N);
     IntegerVector new_dec_v(new_N);
     IntegerVector new_time_v(new_N);
     
     //Transcribe over values for states at NEXT time interval
     int tic = 0;
     for (int i=0; i<num_alive; i++) {
       //only do this for living individuals
       if (alive(i) == 1) {
         new_alive(tic) = alive(i);
         new_state_v(tic) = x_next_rd(i);
         new_res_v(tic) = res_next(i);
         IntegerMatrix istar_t = istar_nr(new_res_v(i)); //Grab the istar for the focal resource
         new_dec_v(tic) = istar_t(new_state_v(i),time_v(i));
         new_time_v(tic) = time_v(i) + 1;
         tic = tic + 1;
       }
     }
     //Initiate new values for recruits
     NumericVector rdraw_v = runif(num_recruits);
     int tic = 0;
     for (int i=N; i<new_N; i++) {
       new_alive(i) = 1; //Recruits are all alive
       new_state_v(i) = xmax; //Initial energetic state is xmax
       new_time_v(i) = 0; //Initial time state is 0;
       //Randomly draw initial resource
       double rdraw = as<double>(rdraw_v(tic)); 
       new_res_v(i) = (int) floor(rdraw*(num_res));
       IntegerMatrix istar_t = istar_nr(new_res_v(i)); //Grab the istar for the focal resource
       new_dec_v(i) = istar_t(new_state_v(i),new_time_v(i)); //Grab the decision for a given state and time t=0;
       int tic = tic + 1;
     }
     
     //How many individuals survive this timestep?
     int num_alive = sum(new_alive);
     
     //Redefine states only for alive individuals
     IntegerVector state_v = new_state_v;
     IntegerVector time_v = new_time_v;
     IntegerVector res_v = new_res_v;
     IntegerVector dec_v = new_dec_v;
     
   } //end time iterations
   
}
