
#include "RooBinnedGauss1d.h"

#include "RooRandom.h"
#include "RooMath.h"
#include "RooHelpers.h"
#include "RooBatchCompute.h"

bool edebug(true) ;
int  evalcount(0) ;

ClassImp(RooBinnedGauss1d);

////////////////////////////////////////////////////////////////////////////////

RooBinnedGauss1d::RooBinnedGauss1d(const char *name, const char *title,
          RooAbsReal& _x,
          RooAbsReal& _mean_x,
          RooAbsReal& _sigma2_x,
          int nbins,
          const RooArgList& bin_edges
          ) :
  RooAbsPdf(name,title),
  x("x","Observable",this,_x),
  mean_x("mean_x","Mean x",this,_mean_x),
  sigma2_x("sigma2_x","Sigma2 x",this,_sigma2_x)
{
   nbins_ = nbins ;
   bin_edges_ = new float(nbins+1) ;
   for ( int bi=0; bi<=nbins; bi++ ) {
      printf(" %3d :  ", bi ) ;
      ((RooConstVar*)bin_edges.at(bi))->Print() ;
      bin_edges_[bi] = ((RooConstVar*)bin_edges.at(bi))->getVal() ;
   } // bi
}

////////////////////////////////////////////////////////////////////////////////

RooBinnedGauss1d::RooBinnedGauss1d(const RooBinnedGauss1d& other, const char* name) :
  RooAbsPdf(other,name),
  x("x",this,other.x),
  mean_x("mean_x",this,other.mean_x),
  sigma2_x("sigma2_x",this,other.sigma2_x)
{
   nbins_ = other.nbins_ ;
   bin_edges_ = new float(nbins_+1) ;
   for ( int bi=0; bi<=nbins_; bi++ ) {
      bin_edges_[bi] = other.bin_edges_[bi] ;
   } // bi

}

////////////////////////////////////////////////////////////////////////////////

Double_t RooBinnedGauss1d::evaluate() const
{
//if ( edebug ) {
//   double xv = x ;
//   double yv = y ;
//   printf("  %5d : x = %f  y = %f\n", evalcount, xv, yv ) ;
//   evalcount++ ;
//   if ( evalcount > 10 ) edebug = false ;
//}

  /////double norm = 1. ;
  /////if ( sigma2_x > 0 ) {
  /////   norm = 1./(sqrt(2*3.14159265*sigma2_x)) ;
  /////}
  ///////double rv = norm * std::exp( -0.5 * (1./det)*(sigma2_y*(x-mean_x)*(x-mean_x) + sigma2_x*(y-mean_y)*(y-mean_y) - 2.*cov_xy*(x-mean_x)*(y-mean_y) ) ) ;
  int bi = 0 ;
  for ( int i=0; i<nbins_; i++ ) {
     if ( x >= bin_edges_[i] && x < bin_edges_[i+1] ) {
        bi = i ;
        break ;
     }
  } // i
  double rv = 0.5 * (
       ( 1. + TMath::Erf( (bin_edges_[bi+1] - mean_x)/(sqrt( 2. * sigma2_x ) ) ) )
     - ( 1. + TMath::Erf( (bin_edges_[bi  ] - mean_x)/(sqrt( 2. * sigma2_x ) ) ) )
  ) ;
  return rv ;
}


