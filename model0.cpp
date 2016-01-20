// [[Rcpp]]
#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp; 


//[[Rcpp::export]]
NumericMatrix F_matrix( double incidence
  , NumericVector sizes
  , NumericVector theta
  , CharacterVector demes
  , IntegerVector NH // length m indicators for each deme 
  , IntegerVector AGE
  , IntegerVector CARE
  , IntegerVector RISK
  , NumericVector nh_wtransm // associated weight for each category
  , NumericVector age_wtransm
  , NumericVector care_wtransm
  , NumericVector risk_wtransm
  , NumericMatrix prRecipMat // pstartstage & age mixing & prisklevel
)
{
	int age, care, risk, nh; 
	int m = demes.size();
	NumericVector w(m, 0.); 
	for (int i = 0; i < m-1 ; i++){ //note not incl src
		age = AGE(i); 
		care = CARE(i);
		risk = RISK(i);
		nh = NH(i); 
//~ std::cout << age << " " << care << " " << risk << " " << nh << std::endl; 
		w(i) = nh_wtransm(nh) * age_wtransm(age) * care_wtransm(care) * risk_wtransm (risk ); 
	}
	w = w / sum(w); 
	NumericVector transm = incidence * w ; 
	NumericMatrix F(m,m ); //, 0.
	std::fill( F.begin(), F.end(), 0. ) ;
	for (int i = 0 ; i < m-1; i++){ //note not incl src
		F(i,_) = transm(i) * prRecipMat(i,_); 
	}
	
	// src
	// br = growth rate - death rate
	F(m-1,m-1) = (theta["srcGrowthRate"] - (1./10./365.)) * sizes(m-1); 
	
	return F; 
}


//[[Rcpp::export]]
NumericMatrix G_matrix( NumericVector sizes
  , NumericVector theta
  , CharacterVector demes
  , IntegerVector NH // length m indicators for each deme 
  , IntegerVector AGE
  , IntegerVector CARE
  , IntegerVector RISK
  , IntegerVector stageprog_recip // destination for migration. NOTE R indices
  , IntegerVector age_recip
  , IntegerVector care_recip
  , NumericVector stageprog_rates //rates for each deme
  , NumericVector age_rates
  , NumericVector care_rates // note these depend on time
)
{
	int m = demes.size() ; 
	NumericMatrix G(m,m); //,0.
	std::fill( G.begin(), G.end(), 0. ) ;
	int recip; 
	int age, care, risk, nh; 
	for (int i = 0; i < (m - 1); i++){
		nh = NH(i) ; 
		care = CARE(i) ;
		age = AGE(i); 
		
		
		recip = stageprog_recip(i) - 1; 
		if (recip>=0){
			G(i,recip) = sizes(i) * stageprog_rates(nh) ;
			if (care==2){
				G(i, recip) *= theta["treatmentEffectiveness"]; 
			}
		}
		
		recip = age_recip(i) - 1;
		if (recip >= 0){
			G(i, recip) = sizes(i) * age_rates( age );
		}
		
		recip = care_recip(i) - 1;
		if (recip >=0){
			G(i, recip) = sizes(i) * care_rates(care ); 
		}
	}
	
	// src
	double srcMigrationRate = theta["srcMigrationRate"];
	for (int i = 0; i < (m-1); i++){
		G(m-1, i) = srcMigrationRate * sizes(i); 
	}
	
	return G;
}

