#include <Rcpp.h>
using namespace Rcpp;

//Function for a single-consumer diet-switching SDP.

// [[Rcpp::export]]
List SDP_single(int tmax, NumericVector res_bs, int cons_bs, int xc, NumericVector rep_gain, 
NumericMatrix f_m, NumericVector mort, List dec_ls, NumericVector rho_vec, NumericMatrix c_learn, NumericVector g_forage, NumericVector c_forage) {
   
   //Establish iniital variables required for the SDP
   int num_res = res_bs.size();
   //How many decisions? Count columns on first decision matrix
   NumericMatrix dec_ex = dec_ls(0);
   int num_dec = dec_ex.ncol();
   int max_enc = f_m.nrow(); //Rows = encounters; Columns = resources
   
   double cons_bs_state = (double) cons_bs;
   //Define xc_state as a double. This variable will be used for energetic calculations
   double xc_state = (double) xc;
   //Define xc, which will be used to locate the critical value in a matrix. 
   //Subtract one, because Cpp indices start at zero.
   xc = xc - 1;
   
   //Build Fitness List, Istar List, and terminal fitness values
   List W_nr(num_res);
   List istar_nr(num_res);
   
   NumericMatrix W_xt(cons_bs,tmax);
   //NumericMatrix istar_xt(cons_bs,tmax-1);
   
   //Initiate terminal fitness values in W_xt matrix
   //Don't forget that the matrix index BEGINS at 0.
   for (int i=xc; i<cons_bs; i++) {
     W_xt(i,tmax-1) = rep_gain(i); //tmax-1 because index starts at zero. So tmax-1 is the terminal time, as t=0 is initial time
   } 
   
   //Place fitness matrix into each list element for future updating
   for (int i=0; i<num_res; i++) {
     W_nr(i) = W_xt;
     //istar_nr(i) = istar_xt;
   }
   
   //Begin Backwards Equation
   for (int r=0; r<num_res; r++) {
     
     Rcpp::Rcout << "r = " << r << std::endl;
     
     NumericMatrix dec_m = dec_ls(r);
     
     NumericMatrix W_r = W_nr(r);
     IntegerMatrix istar_r(cons_bs,tmax-1);
     
     //Begin backwards calculations... start at tmax-1
     for (int t=(tmax-1); t --> 0;) {
       
       //Rcpp::Rcout << "t = " << t << std::endl;
       //int t_state = t - 1;
       
       for (int x=xc; x<cons_bs; x++) {
         
         //Define the energetic state to be used in calculations
         //Added one because the index starts at zero... but we care about the true energetic state, not the index
         double x_state = (double) x;
         x_state = x_state + 1;
         
         NumericVector value(num_dec);
         
         for (int i=0; i<num_dec; i++) {
           
           //Build preference vector for ith decision
           NumericVector pref_vec(num_res);
           for (int j=0; j<num_res; j++) {
             pref_vec(j) = dec_m(j,i);
           }
           
           //Loop over the changes in energetic reserves via consumption of different (t+1) resources
           NumericMatrix xp(num_res,max_enc);
           
           //xp high and low must be integers, because they are used as indices of the fitness matrices
           IntegerMatrix xp_high(num_res,max_enc); 
           IntegerMatrix xp_low(num_res,max_enc);
           
           //But these are not integer matrices
           NumericMatrix q(num_res,max_enc);
           NumericMatrix W(num_res,max_enc);
           
           NumericVector Wk(num_res);
           
           double Fx = 0; //Initialized at zero. Will be added/multiplied as dot product below
           
           for (int rr=0; rr<num_res; rr++) {
             for (int k=0; k<max_enc; k++) {
               
               //Change in Energetic state
               xp(rr,k) = x_state + rho_vec(k)*(g_forage(rr) - (c_forage(rr) + c_learn(r,rr))) - (1-rho_vec(k))*(c_forage(rr) + c_learn(r,rr));
               
               //Establish Boundary conditions en route
               //Lower boundary condition at x-critical
               if (xp(rr,k) < xc_state) {xp(rr,k) = xc_state;}
               //Higher boundary condition at xmax = cons_bs
               if (xp(rr,k) > cons_bs_state) {xp(rr,k) = cons_bs_state;}
               
               //Build xp_high, xp_low, and q matrices en route
               //Estabish Fitness en route
               if ((xp(rr,k) < cons_bs_state) && (xp(rr,k) > xc_state)) {
                 
                 xp_low(rr,k) = (int) floor(xp(rr,k));
                 xp_high(rr,k) = (int) xp_low(rr,k) + 1;
                 
                 //Make sure that we do not have an (integer - double)
                 double xp_h = (double) xp_high(rr,k);
                 q(rr,k) = xp_h - xp(rr,k);
                 
                 //Prep for indexing
                 //Adjust xp_low and xp_high to signify indices rather than energetic values
                 xp_low(rr,k) = xp_low(rr,k) - 1;
                 xp_high(rr,k) = xp_high(rr,k) - 1;
                 
                 W(rr,k) = q(rr,k)*W_r(xp_low(rr,k),t+1) + (1-q(rr,k))*W_r(xp_high(rr,k),t+1);
               }
               
               if (xp(rr,k) == xc_state) {
                 
                 xp_low(rr,k) = (int) xc_state;
                 xp_high(rr,k) = (int) xc_state+1;
                 q(rr,k) = 1;
                 
                 //Prep for indexing
                 //Adjust xp_low and xp_high to signify indices rather than energetic values
                 xp_low(rr,k) = xp_low(rr,k) - 1;
                 xp_high(rr,k) = xp_high(rr,k) - 1;
                 
                 W(rr,k) = q(rr,k)*W_r(xp_low(rr,k),t+1) + (1-q(rr,k))*W_r(xp_high(rr,k),t+1);
               }
               
               if (xp(rr,k) == cons_bs_state) {
                 
                 xp_low(rr,k) = (int) cons_bs_state - 1;
                 xp_high(rr,k) = (int) cons_bs_state;
                 q(rr,k) = 0;
                 
                 //Prep for indexing
                 //Adjust xp_low and xp_high to signify indices rather than energetic values
                 xp_low(rr,k) = xp_low(rr,k) - 1;
                 xp_high(rr,k) = xp_high(rr,k) - 1;
                 
                 W(rr,k) = q(rr,k)*W_r(xp_low(rr,k),t+1) + (1-q(rr,k))*W_r(xp_high(rr,k),t+1);
               }
               
               //Vector multiplication: fitness vector (over k) and probability of finding k resources (over k)
               //Note that the values for the W and f.m matrices are transposed relative to each other
               double vec = W(rr,k) * f_m(k,rr);
               Wk(rr) += vec; // The same as Wk(rr) = Wk(rr) + vec
               
             } // end k
             
             
             //Vector multiplication over preference probabilities and Wk ** over rr **
             double vec2 = pref_vec(rr) * (rep_gain(x) + (1-mort(rr))*Wk(rr));
             Fx += vec2;
             
             //Rcpp::Rcout << "Fx = " << Fx << std::endl;
             
           } // end rr
           
           //Record fitness value for decision i
           value(i) = Fx;
           
         } // end i (decision loop)
         
         //Find maximum value
         int istar = which_max(value);
         
         istar_r(x,t) = istar + 1; //The +1 is there because we are going from 'index' to 'decision number'
         W_r(x,t) = value(istar);
         
       } //End x iterations
       
     } //End t iterations
     
     //Replace updated W_nr
     //Update istar_nr list
     W_nr(r) = W_r;
     istar_nr(r) = istar_r;
     
   } //End r iterations
   
   List output(2);
   output(0) = W_nr;
   output(1) = istar_nr;
   
   return (output);
   
} //End cpp function
