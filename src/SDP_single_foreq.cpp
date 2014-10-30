#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
List SDP_single_foreq(int tmax, NumericVector res_bs, int cons_bs, int xc, NumericVector rep_gain, 
NumericMatrix f_m, NumericVector mort, List dec_ls, NumericVector rho_vec, NumericMatrix c_learn, 
NumericVector g_forage, NumericVector c_forage, List W_nr, List istar_nr, int N, int tsim) {
   
   //Establish iniital variables required for the SDP
   int num_res = res_bs.size();
   //How many decisions? Count columns on first decision matrix
   NumericMatrix dec_ex = dec_ls(0);
   int num_dec = dec_ex.ncol();
   int max_enc = f_m.nrow(); //Rows = encounters; Columns = resources
   
   double xmax = (double) cons_bs;
   
   //Initialize the record-keeping of states
   //Record the state vector for each individual
   //This records the state of each individual... values will change over time
   IntegerVector state_v(N);
   for (i=0; i<N; i++) {
     state_v(i) = xmax; //The initial state is full health
   }
   //Record the time vector for each individual
   IntegerVector time_v(N);
   for (i=0; i<N; i++) {
     time_v(i) = 1; //The intitial time is 1 of course
   }
   //Record the current resource for each individual
   IntegerVector res_v(N);
   NumericVector rdraw_v = runif(N);
   for (i=0; i<N; i++) {
    double rdraw = as<double>(rdraw_v(i));
    res_v(i) = (int) floor(rdraw*(num_res)); //The initial resource is randomly drawn for each individual
   }
   
   //What are the optimal decisions for the current individual-level states?
   IntegerVector dec_v(N);
   for (i=0; i<N; i++) {
     istar_t = istar_nr(res_v(i)); //Grab the istar for the focal resource
     dec_v(i) = istar_t(state_v(i),0); //Grab the decision for a given state and time t=0;
   }
   
   
   
   //Begin time iteration (this is the simulation time)
   for (t=0; t<tsim; t++) {
     
     //Begin Individual iterations
     for (n=0; n<N; n++) {
       
       //Define states for the individual
       int ind_x = state_v(n);
       int ind_t = time_v(n);
       int ind_r = res_v(n);
       int ind_d = dec_v(n);
       
       //Choose next resource as a function of the preference probability distribution
       NumericMatrix pref_prob_m = dec_ls(ind_r);
       NumericVector pref_prob(num_res);
       IntegerVector num_prob(num_res);
       for (i=0;i<num_res; i++) {
         pref_prob(i) = pref_prob_m(i,ind_d);
         num_prob(i) = floor(pref_prob(i)*1000);
       } 
       int tot_num = sum(num_prob);
       IntegerVector k_bin(tot_num);
       int tic = 0;
       for (i=0; i<num_res; i++) {
         int num_k = num_prob(i);
         for (j=0; j<num_k; j++) {
           k_bin(tic) = i;
           tic = tic + 1;
         }
       }
       NumericVector rdraw_v = runif(1);
       double rdraw = as<double>(rdraw_v;
       int draw = (int) floor(rdraw*(tot_num));
       int k = k_bin(draw);
       
       
       
       
     } //end individual iterations
     
   } //end time iterations
   
}
