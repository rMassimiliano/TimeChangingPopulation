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

//allocation function g() used in Doubly adaptive biased coin
double allocationFunction(const double x, const double y, const double gamma)
{
  double res;
  if(x ==0)
  {
     res =  1.0;
  }
  if(x == 1)
  {
     res = 0.0;
  }
  else
  {
    res = (y*pow((y/x),gamma))/ ( (y*pow((y/x),gamma)) + (1-y)* (pow((1-y)/(1-x),gamma)) ) ;
  }
  return res;
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

/*
// implement 2 arm randomization as in Thall et al. 2015
// [[Rcpp::export]]
mat Thall2015(mat potential, double c = 1.0, double outcome_delay =0, bool ret_prob = 0)
{
 //indices
 int N = potential.n_rows;
 int A = potential.n_cols-2; 
 if(A > 2)
 {
    stop("A>2; This procedure is for two-arm trials");
 }
 //quantities
 vec arm_prob(A);
 uvec indx;
 mat  accrued_data;
 mat newPotential;
 umat SuffStat(2,A);
 vec PrPTE(A-1);
 double sum_prob_a;



 SuffStat.zeros();

 // initialize arm_prob to 1/A
 for(int a=0;a<A;a++)
 {
  arm_prob[a] = pow(A,-1);
 }

 // discrete generator for the cpp standard library
 std::random_device rd;
 std::mt19937 gen(rd());


 // assign the first patients 
  std::discrete_distribution<> sample_arm(arm_prob.begin(),arm_prob.end()); 
    potential(0,1) = sample_arm(gen) + 1;



 if(ret_prob)
 {
  newPotential.set_size(potential.n_rows,potential.n_cols + 1);

  for(int i=0;i < potential.n_rows; i++)
  {
	  for(int j=0; j < potential.n_cols; j++)
	  {
            newPotential(i,j) = potential(i,j);
	  }
  }

   newPotential(0,newPotential.n_cols - 1) = 0.5;

 }
 else
 {
  newPotential.set_size(potential.n_rows,potential.n_cols);

  for(int i=0;i < potential.n_rows; i++)
  {
	  for(int j=0; j < potential.n_cols; j++)
	  {
            newPotential(i,j) = potential(i,j);
	  }
  }
 }



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
          
        // add 1 when observe a positive outcome.
         if(accrued_data(j, accrued_data(j,1) -1 + 2 ) ==1)
        {
         SuffStat(1,accrued_data(j,1) -1) +=1;
        }
       }
       // compute Pr(Treatment >  control| Accrued data )
       PrPTE = posterior_with_control(SuffStat);
       arm_prob[1] = pow(as_scalar(PrPTE),c);
       arm_prob[0] = pow(1.0 - as_scalar(PrPTE),c);
       // normalize probability
       sum_prob_a = sum(arm_prob);
       arm_prob[0] = arm_prob[0]/sum_prob_a;
       arm_prob[1] = arm_prob[1]/sum_prob_a;
      }
      // assign patients to arm
      std::discrete_distribution<> sample_arm(arm_prob.begin(),arm_prob.end()); 

      potential(i,1) = sample_arm(gen) + 1;
      newPotential(i,1) = potential(i,1);

      if(ret_prob)
      {
      	newPotential(i,newPotential.n_cols-1) =  arm_prob[1];
      }

      }
      if(ret_prob) return newPotential;
      else return potential;
}

*/





//' Simulate a two-arm trial using Efron's biased coin design 
//'
//' 
//'
//' @param outcome
//' A matrix with one potential outcome data for each arm. The first column is the enrolment time, the second a place holder for the arm the third and fourth potential outcome for the SOC and the experimental arm
//  
//' 
//'
//' @param alpha and beta, this are design parameter
//' 
//' 
//'
//' @param ret_prob 
//' if =1 return the probability to be assigned to the experimental arm for each of the subject.
//'                   //' @return
//' The function returns a outcome matrix with the allocations. 
//'
//'
//'@details This function implement the Efron's biased coin design
//'

//includes minor changes form the one implemented in the R package AddArmToTrial; refer to \url{http://bcb.dfci.harvard.edu/∼steffen/software.html}.
//' @export
// implement 2 arm biased coin design
// [[Rcpp::export]]
mat adaptiveBiasedCoin(mat potential, double alpha  = 0.0, double beta = 1.0, bool ret_prob = 0)
{
 //indices
 int N = potential.n_rows;
 int A = potential.n_cols-2; 
 if(A > 2)
 {
    stop("A>2; This procedure is for two-arm trials");
 }
 //quantities
 vec arm_prob(A);
 uvec indx;
 mat  accrued_data;
 mat newPotential;
 umat SuffStat(2,A);
 vec PrPTE(A-1);
 double sum_prob_a;

 // initialize arm_prob to 1/A
 for(int a=0;a<A;a++)
 {
  arm_prob[a] = pow(A,-1);
 }




 // discrete generator for the cpp standard library
 std::random_device rd;
 std::mt19937 gen(rd());


 // assign the first patients 
  std::discrete_distribution<> sample_arm(arm_prob.begin(),arm_prob.end()); 
    potential(0,1) = sample_arm(gen) + 1;


 if(ret_prob)
 {
  newPotential.set_size(potential.n_rows,potential.n_cols + 1);

  for(int i=0;i < potential.n_rows; i++)
  {
	  for(int j=0; j < potential.n_cols; j++)
	  {
            newPotential(i,j) = potential(i,j);
	  }
  }

   newPotential(0,newPotential.n_cols - 1) = 0.5;

 }
 else
 {
  newPotential.set_size(potential.n_rows,potential.n_cols);

  for(int i=0;i < potential.n_rows; i++)
  {
	  for(int j=0; j < potential.n_cols; j++)
	  {
            newPotential(i,j) = potential(i,j);
	  }
  }
 }



    // start BAR
    for(int i=1; i<N; i++)
    {
      indx = linspace<uvec>(0,i-1 ,i);
      SuffStat.zeros();
      accrued_data =  potential.rows(indx);
      // compute sufficient statistics for the ind
       for(uword j = 0; j < accrued_data.n_rows; j++)
       {
         SuffStat(0,accrued_data(j,1)-1) +=1;
          
        // add 1 when observe a positive outcome.
         if(accrued_data(j, accrued_data(j,1) -1 + 2 ) ==1)
        {
         SuffStat(1,accrued_data(j,1) -1) +=1;
        }
       }
       arm_prob[1] = (alpha + beta* SuffStat(0,0))/(2*alpha + beta*(i-1));
       arm_prob[0] = 1 - arm_prob[1];
      // assign patients to arm
      std::discrete_distribution<> sample_arm(arm_prob.begin(),arm_prob.end()); 
      potential(i,1)    = sample_arm(gen) + 1;
      newPotential(i,1) = potential(i,1);

      if(ret_prob)
      {
      	newPotential(i,newPotential.n_cols-1) =  arm_prob[1];
      }

    }
    if(ret_prob) return newPotential;
    else return potential;
}

// implements 2 arm bar randomization similar to Thall et a. 2015
// we change randomization any M observed patients, default =1
// [[Rcpp::export]]
mat bar2arms(mat potential, double c = 1.0,  int M = 1, bool ret_prob = 0)
{
 //indices
 int N = potential.n_rows;
 int A = potential.n_cols-2; 
 if(A > 2)
 {
    stop("A>2; This procedure is for two-arm trials");
 }
 //quantities
 vec arm_prob(A);
 uvec indx;
 mat  accrued_data;
 mat  newPotential;
 umat SuffStat(2,A);
 vec  PrPTE(A-1);
 double sum_prob_a;


 SuffStat.zeros();

 // initialize arm_prob to 1/A
 for(int a=0;a<A;a++)
 {
  arm_prob[a] = pow(A,-1);
 }

 // discrete generator for the cpp standard library
 std::random_device rd;
 std::mt19937 gen(rd());


 // assign the first patients 
  std::discrete_distribution<> sample_arm(arm_prob.begin(),arm_prob.end()); 
    potential(0,1) = sample_arm(gen) + 1;

 if(ret_prob)
 {
  newPotential.set_size(potential.n_rows,potential.n_cols + 1);

  for(int i=0;i < potential.n_rows; i++)
  {
	  for(int j=0; j < potential.n_cols; j++)
	  {
            newPotential(i,j) = potential(i,j);
	  }
  }

   newPotential(0,newPotential.n_cols - 1) = 0.5;

 }
 else
 {
  newPotential.set_size(potential.n_rows,potential.n_cols);

  for(int i=0;i < potential.n_rows; i++)
  {
	  for(int j=0; j < potential.n_cols; j++)
	  {
            newPotential(i,j) = potential(i,j);
	  }
  }
 }





    // start BAR
    for(int i=1; i<N; i++)
    {
      // determine indices of observed responses
      if(i % M == 0)
      {
       indx = linspace<uvec>(0,i-1 ,i);
       SuffStat.zeros();
       accrued_data =  potential.rows(indx);
      // compute sufficient statistics for the ind
       for(uword j = 0; j < accrued_data.n_rows; j++)
       {
         SuffStat(0,accrued_data(j,1)-1) +=1;
          
        // add 1 when observe a positive outcome.
         if(accrued_data(j, accrued_data(j,1) -1 + 2 ) ==1)
        {
         SuffStat(1,accrued_data(j,1) -1) +=1;
        }
       }
       // compute Pr(Treatment >  control| Accrued data )
       PrPTE = posterior_with_control(SuffStat);
       arm_prob[1] = pow(as_scalar(PrPTE),c);
       arm_prob[0] = pow(1.0 - as_scalar(PrPTE),c);
       // normalize probability
       sum_prob_a = sum(arm_prob);
       arm_prob[0] = arm_prob[0]/sum_prob_a;
       arm_prob[1] = arm_prob[1]/sum_prob_a;

        //Rcout << "I am at i = " << i << " and I am changing the randomization to " << arm_prob << std::endl;
      }

      // assign patients to arm
      std::discrete_distribution<> sample_arm(arm_prob.begin(),arm_prob.end()); 
      potential(i,1)    = sample_arm(gen) + 1;
      newPotential(i,1) = potential(i,1); 
      if(ret_prob)
      {
      	newPotential(i,newPotential.n_cols-1) =  arm_prob[1];
      }
      }
      if(ret_prob) return newPotential;
      else return potential;
}


// implements 2 arm DABC bar randomization similar to Duan and Hu (2009)
// we change randomization any M observed patients, default =1
// [[Rcpp::export]]
mat doublyAdaptiveBiasedCoin(mat potential, double gamma, int M =1, bool ret_prob =0)
{

 //indices
 int N = potential.n_rows;
 int A = potential.n_cols-2; 
 if(A > 2)
 {
    stop("A>2; This procedure is for two-arm trials");
 }
 //quantities
 vec arm_prob(A);
 uvec indx;
 mat  accrued_data;
 mat newPotential;
 umat SuffStat(2,A);
 SuffStat.zeros();
 vec PrPTE(A-1);
 double sum_prob_a;
 double p0,p1,rho;

 // initialize arm_prob to 1/A
 for(int a=0;a<A;a++)
 {
  arm_prob[a] = pow(A,-1);
 }

 // discrete generator for the cpp standard library
 std::random_device rd;
 std::mt19937 gen(rd());

 // assign the first patients 
  std::discrete_distribution<> sample_arm(arm_prob.begin(),arm_prob.end()); 
  potential(0,1) = sample_arm(gen) + 1;

 if(ret_prob)
 {
  newPotential.set_size(potential.n_rows,potential.n_cols + 1);

  for(int i=0;i < potential.n_rows; i++)
  {
	  for(int j=0; j < potential.n_cols; j++)
	  {
            newPotential(i,j) = potential(i,j);
	  }
  }

   newPotential(0,newPotential.n_cols - 1) = 0.5;

 }
 else
 {
  newPotential.set_size(potential.n_rows,potential.n_cols);

  for(int i=0;i < potential.n_rows; i++)
  {
	  for(int j=0; j < potential.n_cols; j++)
	  {
            newPotential(i,j) = potential(i,j);
	  }
  }
 }


  //start response-adaptive allocation
    for(int i=1; i<N; i++)
    {
      // determine indices of observed responses
      if(i % M == 0)
      {
       indx = linspace<uvec>(0,i-1 ,i);
       SuffStat.zeros();
       accrued_data =  potential.rows(indx);
      // compute sufficient statistics for the ind
       for(uword j = 0; j < accrued_data.n_rows; j++)
       {
         SuffStat(0,accrued_data(j,1)-1) +=1;
          
        // add 1 when observe a positive outcome.
         if(accrued_data(j, accrued_data(j,1) -1 + 2 ) ==1)
        {
         SuffStat(1,accrued_data(j,1) -1) +=1;
        }
       }
       // compute the MLE for the two groups
       p0 = (SuffStat(1,0) + 0.5)/ (SuffStat(0,0) +1.0);
       p1 = (SuffStat(1,1) + 0.5)/ (SuffStat(0,1) + 1.0); 
       //Rcout << "(p0,p1) = ("<< p0 << ", " << p1 <<")\n";   

       //Neyman allocation 
       rho = pow(p1*(1-p1), 0.5)/(pow(p1*(1-p1), 0.5) + pow(p0*(1-p0), 0.5));
       //Rcout << "rho = " << rho << "\n";
       // probabilities
          
          //Rcout <<"SS = " << SuffStat(0,1)  << "i-1 " << i-1 << "\n";
          //Rcout <<"SS/(i-1) = " << (1.0*SuffStat(0,1))/(1.0*(i-1))<< "\n";

       arm_prob[1] = allocationFunction((1.0*SuffStat(0,1))/(1.0*(i-1)), rho, gamma);
       arm_prob[0] = 1- arm_prob[1];

        //Rcout << "I am at i = " << i << " and I am changing the randomization to " << arm_prob << std::endl;
      }

      // assign patients to arm
      std::discrete_distribution<> sample_arm(arm_prob.begin(),arm_prob.end()); 
      potential(i,1) = sample_arm(gen) + 1;
      newPotential(i,1) = potential(i,1);
      if(ret_prob)
      {
      	newPotential(i,newPotential.n_cols-1) =  arm_prob[1];
      }
      }
      if(ret_prob) return newPotential;
      else return potential;
}




// implements 2 arm RPW randomization similar to Wei and Durham (1978)
// we change randomization any M observed patients, default =1
// [[Rcpp::export]]
mat randomizedPlayTheWinner(mat potential, int beta = 1, int alpha = 1,  int M =1, bool ret_prob = 0)
{

 //indices
 int N = potential.n_rows;
 int A = potential.n_cols-2; 
 if(A > 2)
 {
    stop("A>2; This procedure is for two-arm trials");
 }
 //quantities
 vec arm_prob(A);
 uvec indx;
 mat  accrued_data;
 mat newPotential;
 umat SuffStat(2,A);
 SuffStat.zeros();
 vec urnComposition(A);
 double sum_prob_a;
 int outcome, arm;

 // initialize arm_prob to 1/A
 for(int a=0;a<A;a++)
 {
  arm_prob[a] = pow(A,-1);
  urnComposition(a) = 1;
 }

 // discrete generator for the cpp standard library
 std::random_device rd;
 std::mt19937 gen(rd());

 // assign the first patients 
  std::discrete_distribution<> sample_arm(arm_prob.begin(),arm_prob.end()); 
  potential(0,1) = sample_arm(gen) + 1;

 if(ret_prob)
 {
  newPotential.set_size(potential.n_rows,potential.n_cols + 1);

  for(int i=0;i < potential.n_rows; i++)
  {
	  for(int j=0; j < potential.n_cols; j++)
	  {
            newPotential(i,j) = potential(i,j);
	  }
  }

   newPotential(0,newPotential.n_cols - 1) = 0.5;

 }
 else
 {
  newPotential.set_size(potential.n_rows,potential.n_cols);

  for(int i=0;i < potential.n_rows; i++)
  {
	  for(int j=0; j < potential.n_cols; j++)
	  {
            newPotential(i,j) = potential(i,j);
	  }
  }
 }

  //start response-adaptive allocation
    for(int i=1; i<N; i++)
    {
      // determine indices of observed responses
      if(i % M == 0)
      {
       indx = linspace<uvec>(0,i-1 ,i);
       SuffStat.zeros();
       accrued_data =  potential.rows(indx);
      // compute sufficient statistics for the ind
       for(uword j = 0; j < accrued_data.n_rows; j++)
       {
         SuffStat(0,accrued_data(j,1)-1) +=1;
          
        // add 1 when observe a positive outcome.
         if(accrued_data(j, accrued_data(j,1) -1 + 2 ) ==1)
        {
         SuffStat(1,accrued_data(j,1) -1) +=1;
        }
       }
       // update the urn
        // is success add beta balls for the treatments
        outcome = potential(i-1, potential(i-1,1)-1+2 );
        arm     = potential(i-1,1) -1;
        // Rcout << "(outcome,arm, abs(arm -1)) = ("<< outcome << ", " << arm << ", " << abs(arm -1) << ") \n";  
        if(outcome == 1 )
        {
         urnComposition(arm) += beta; 
        }
        if(outcome == 0 )
        {
         urnComposition(abs( arm -1)) += beta; 
        } 
        // Rcout << "urnComposition = " <<  urnComposition << std::endl;
        
       //Rcout << "rho = " << rho << "\n";
       // probabilities
       arm_prob[1] = urnComposition(1)/sum(urnComposition);
       arm_prob[0] = 1- arm_prob[1];

        //Rcout << "I am at i = " << i << " and I am changing the randomization to " << arm_prob << std::endl;
      }

      // assign patients to arm
      std::discrete_distribution<> sample_arm(arm_prob.begin(),arm_prob.end()); 
      potential(i,1) = sample_arm(gen) + 1;
      newPotential(i,1) = potential(i,1);
      if(ret_prob)
      {
      	newPotential(i,newPotential.n_cols-1) =  arm_prob[1];
      }
      }
      if(ret_prob) return newPotential;
      else return potential;
}



// implements 2 arm bar randomization similar to the Battle trial Zhou et al 2017
// we change randomization any M observed patients, default =1
// [[Rcpp::export]]
mat barBattle2arms(mat potential, double alpha = 1.0, double beta = 1.0, int M = 1, double minProb = 0.1, bool ret_prob = 0)
{
 //indices
 int N = potential.n_rows;
 int A = potential.n_cols-2; 
 if(A > 2)
 {
    stop("A>2; This procedure is for two-arm trials");
 }
 //quantities
 vec arm_prob(A);
 uvec indx;
 mat  accrued_data;
 mat newPotential;
 umat SuffStat(2,A);
 vec PrPTE(A-1);
 double sum_prob_a;


 SuffStat.zeros();

 // initialize arm_prob to 1/A
 for(int a=0;a<A;a++)
 {
  arm_prob[a] = pow(A,-1);
 }

 // discrete generator for the cpp standard library
 std::random_device rd;
 std::mt19937 gen(rd());


 // assign the first patients 
  std::discrete_distribution<> sample_arm(arm_prob.begin(),arm_prob.end()); 
    potential(0,1) = sample_arm(gen) + 1;

 if(ret_prob)
 {
  newPotential.set_size(potential.n_rows,potential.n_cols + 1);

  for(int i=0;i < potential.n_rows; i++)
  {
	  for(int j=0; j < potential.n_cols; j++)
	  {
            newPotential(i,j) = potential(i,j);
	  }
  }

   newPotential(0,newPotential.n_cols - 1) = 0.5;

 }
 else
 {
  newPotential.set_size(potential.n_rows,potential.n_cols);

  for(int i=0;i < potential.n_rows; i++)
  {
	  for(int j=0; j < potential.n_cols; j++)
	  {
            newPotential(i,j) = potential(i,j);
	  }
  }
 }



    // start BAR
    for(int i=1; i<N; i++)
    {
      // determine indices of observed responses
      if(i % M == 0)
      {
       indx = linspace<uvec>(0,i-1 ,i);
       SuffStat.zeros();
       accrued_data =  potential.rows(indx);
      // compute sufficient statistics for the ind
       for(uword j = 0; j < accrued_data.n_rows; j++)
       {
         SuffStat(0,accrued_data(j,1)-1) +=1;
          
        // add 1 when observe a positive outcome.
         if(accrued_data(j, accrued_data(j,1) -1 + 2 ) ==1)
        {
         SuffStat(1,accrued_data(j,1) -1) +=1;
        }
       }
       // compute Pr(Treatment >  control| Accrued data )
       arm_prob[0] = ((1.0*SuffStat(1,0)) + alpha)/(alpha + beta + (1.0*SuffStat(0,0)) ); 
       arm_prob[1] = ((1.0*SuffStat(1,1)) + alpha)/(alpha + beta + (1.0*SuffStat(0,1)) ); 
      // Rcout << "posterior mean " << arm_prob << std::endl; 

       // normalize probability
       sum_prob_a = sum(arm_prob);
       arm_prob[0] = arm_prob[0]/sum_prob_a;
       arm_prob[1] = arm_prob[1]/sum_prob_a;
      // Rcout << "randomization (uncorrected) " << arm_prob <<  std::endl;
       if(arm_prob[0] < minProb)
       {
        arm_prob[1] =  (minProb - arm_prob[0]);
        arm_prob[0] = minProb; 
       }
       if(arm_prob[1] < minProb)
       {
        arm_prob[0] -= (minProb - arm_prob[1]);
        arm_prob[1] = minProb; 
       }

        //Rcout << "I am at i = " << i << " and I am changing the randomization to " << arm_prob << std::endl;
      }

      // assign patients to arm
      std::discrete_distribution<> sample_arm(arm_prob.begin(),arm_prob.end()); 
      potential(i,1) = sample_arm(gen) + 1;
      newPotential(i,1) = potential(i,1);
      if(ret_prob)
      {
      	newPotential(i,newPotential.n_cols-1) =  arm_prob[1];
      }

      }
      if(ret_prob) return newPotential;
      else return potential;
}
