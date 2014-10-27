#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
List SDP_single(int tmax, NumericVector res_bs, double cons_bs, int xc, NumericVector rep_gain, 
NumericVector f_m, NumericVector mort, List dec_ls, NumericVector rho_vec, NumericMatrix c_learn) {
   
   //Establish iniital variables required for the SDP
   int num_res = res_bs.size();
   int num_dec = dec_ls.size();
   int max_enc = f_m.size(2); //Make sure this is right :)
   
   //Adjust indices because the count now starts at zero
   double xc_state = as<double>(xc);
   xc = xc - 1;
   
   //Build Fitness List, Istar List, and terminal fitness values
   List W_nr(num_res);
   List istar_nr(num_res);
   
   NumericMatrix W_xt(cons_bs,tmax);
   NumericMatrix istar_xt(cons_bs,tmax-1);
   
   //Initiate terminal fitness values in W_xt matrix
   //Don't forget that the matrix index BEGINS at 0.
   for (int i=xc; i<cons_bs, i++) {
     W_xt(i,tmax) = rep.gain(i);
   } 
   
   for (int i=0; i<num_res; i++) {
     W_nr(i) = W_xt;
     istar_nr(i) = istar_xt;
   }
   
   
   
   //Begin Backwards Equation
   for (int r=0; r<num_res; r++) {
     
     NumericMatrix dec_m = dec_ls(r);
     
     NumericMatrix W_r = W_nr[r];
     NumericMatrix istar_r = istar_nr[r];
     
     for (int t=(tmax-1); t --> 0;) {
       
       for (int x=xc; x<cons_bs; x++) {
         
         NumericVector value(num_dec);
         
         for (int i=0; i<num_dec; i++) {
           
           //Build preference vector for ith decision
           NumericVector pref_vec(num_res);
           
           for (int j=0; j<num_res; j++) {
             pref_vec(j) = dec_m(j,i);
           }
           
           //Loop over the changes in energetic reserves via consumption of different (t+1) resources
           NumericMatrix xp(num_res,max_enc);
           NumericMatrix xp_high(num_res,max_enc);
           NumericMatrix xp_low(num_res,max_enc);
           NumericMatrix q(num_res,max_enc);
           
           NumericVector Wk(num_res);
           double Fx;
           
           for (int rr=0; rr<num_res; rr++) {
             for (int k=0; k<max_enc; k++) {
               
               //Change in Energetic state
               xp(rr,k) = x + rho_vec(k)*(g_forage(rr) - (c_forage(rr) + c_learn(r,rr))) - (1-rho_vec(k))*(c_forage(rr) + c_learn(r,rr));
               
               //Establish Boundary conditions en route
               if (xp(rr,k) < xc+1) {xp(rr,k) = xc_state;} //needs to be xc+1 because we modified xc for indexing
               if (xp(rr,k) > cons_bs) {xp(rr,k) = cons_bs;}
               
               //Build xp_high, xp_low, and q matrices en route
               //Estabish Fitness en route
               if ((xp(rr,k) < cons_bs) && (xp(rr,k) > xc_state)) {
                 
                 xp_low(rr,k) =floor(xp(rr,k));
                 xp_high(rr,k) = xp_low + 1;
                 q(rr,k) = xp_high(rr,k) - xp(rr,k);
                 
                 W(rr,k) = q(rr,k)*W_r(xp_low(rr,k),t+1) + (1-q)*W_r(xp_high,t+1);
               }
               
               if (xp(rr,k) == xc_state) {
                 
                 xp_low(rr,k) = xc_state;
                 xp_high(rr,k) = xc_state+1;
                 q(rr,k) = 1;
                 
                 W(rr,k) = q(rr,k)*W_r(xp_low(rr,k),t+1) + (1-q)*W_r(xp_high,t+1);
               }
               if (xp(rr,k) == cons_bs) {
                 
                 xp_low(rr,k) = cons_bs -1;
                 xp_high(rr,k) = cons_bs;
                 q(rr,k) = 0;
                 
                 W(rr,k) = q(rr,k)*W_r(xp_low(rr,k),t+1) + (1-q)*W_r(xp_high,t+1);
               }
               
               //Vector multiplication: fitness vector (over k) and probability of finding k resources (over k)
               double vec = W(rr,k) * f.m(k,rr);
               Wk(rr) += vec; // The same as Wk(rr) = Wk(rr) + vec
               
             } // end k
             
             //Vector multiplication over preference probabilities and Wk over rr
             double vec2 = (rep_gain(x) + (1-mort(rr))*Wk(rr)) * pref_vec(rr);
             Fx += vec2;
           
           } // end rr
           
           value(i) = Fx;
           
         } // end i (decision loop)
         
         //Find maximum value
         int istar = which_max(value);
         
         istar_r(x,t) = istar;
         W_r(x,t) = value(istar);
         
       } //End x iterations
       
     } //End t iterations
     
     //Update W_nr and istar_nr list
     W_nr[r] = W_r;
     istar_nr[r] = istar_r;
     
   } //End r iterations
   
   
}
