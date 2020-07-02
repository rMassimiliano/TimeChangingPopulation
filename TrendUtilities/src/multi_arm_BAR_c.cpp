 #include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
vec posterior_with_control(umat Outcome)
{
 vec grid = linspace(0.05, 0.95, 20);
 int A = Outcome.n_cols - 1; 
 mat posterior_param(2,A +1); // first is controll
 vec Pr0(grid.n_elem);
 mat Pra(A,grid.n_elem);

 vec F0(grid.n_elem);
 vec tmp_F0(grid.n_elem);

 // compute posterior parameters
 for(unsigned int a = 0; a< Outcome.n_cols; a++)
 {
   posterior_param(0,a) = 1 + Outcome(1,a);    
   posterior_param(1,a) = 1 + (Outcome(0,a) - Outcome(1,a));
 }
 // compute discrete beta approximation over a grid
 for(unsigned int x=0; x< grid.n_elem; x++)
 {
   Pr0[x]  = R::dbeta(grid[x],posterior_param(0,0), posterior_param(1,0),0);
   for(int a = 0; a < A; a++ )
   {
   Pra(a,x) = R::dbeta(grid[x],posterior_param(0,a + 1), posterior_param(1,a + 1),0);
   }
 }

 
   // standardize
   Pr0 = Pr0/sum(Pr0);
   for(int a = 0; a<A; a++ )
   {
   Pra.row(a) = Pra.row(a)/sum(Pra.row(a));
   }
   // compute cumulative distribution function for controls
   tmp_F0 = cumsum(Pr0);

   F0[0] = 0.0;
   for(unsigned int i = 1; i < grid.n_elem; i++)
   {
     F0[i] = tmp_F0[i - 1];
   }

   return (Pra * F0);
}







mat multi_arm_BAR_c (mat potential,
                     vec rand,
                     double omega,
                     double N_min = 0, // set to double for computational reasons
                     double outcome_delay =0
                     )
{
 //indices
 int N = potential.n_rows;
 int A = potential.n_cols-2;

 //quantities
 vec arm_prob(A);
 uvec indx;
 mat  accrued_data;
 umat SuffStat(2,A);
 vec PrPTE(A-1);
 double h_i;
 double min_cor;
 double sum_prob_a;


 // initialize arm_prob to 1/A
 SuffStat.zeros();

 for(int a=0;a<A;a++)
 {
  arm_prob[a] = pow(A,-1);
 }

 // discrete generator for the cpp standard library
 std::random_device rd;
 std::mt19937 gen(rd());


 // generate the first A at random
  std::discrete_distribution<> sample_arm(arm_prob.begin(),arm_prob.end()); 
    potential(0,1) = sample_arm(gen) + 1;

    // start BAR
    for(int i=1; i<N; i++)
    {
      // determine indices of observed responses
      indx = find((potential(i,0) - potential(span(0,i-1),0))>= outcome_delay);

      if(indx.n_elem > 0)
      {
       SuffStat.zeros();
       accrued_data =  potential.rows(indx);

      // compute sufficient statistics for the ind
      // 
       
       for(uword j = 0; j < accrued_data.n_rows; j++)
       {
         SuffStat(0,accrued_data(j,1)-1) +=1;

         if(accrued_data(j, accrued_data(j,1) -1 + 2 ) ==1)
        {
         SuffStat(1,accrued_data(j,1) -1) +=1;
        }
       }
      
      //
      //
      
      // update probability 
       PrPTE = posterior_with_control(SuffStat);
       h_i  = rand[0] * pow((i *pow(N,-1)),rand[1]);
       arm_prob[0] = as_scalar(exp(rand[2] * (SuffStat.cols(1,A-1)).max() -SuffStat(0,0)));
       arm_prob(span(1,A-1)) = exp(log(PrPTE)*h_i);

       // multiply by exponential function to approximately guarantee N_min patinent for a>0
       for(int a=1;a<A;a++)
       {
          min_cor = N_min - SuffStat(0,a);
          if(min_cor <= 0)
          {
            min_cor = 1;
          }
          else
          {
           min_cor = exp(omega *min_cor);
          }  
         arm_prob[a]  = min_cor*arm_prob[a];
       }

       // normalize probability
       sum_prob_a = sum(arm_prob);
       for(int a=0;a<A;a++)
       {
        arm_prob[a] =  arm_prob[a]/sum_prob_a;
       }

     }
      //sample arm for patient i
      // generate the first A at random
      std::discrete_distribution<> sample_arm(arm_prob.begin(),arm_prob.end()); 
      potential(i,1) = sample_arm(gen) + 1;
    }
 
    return potential;
   }

