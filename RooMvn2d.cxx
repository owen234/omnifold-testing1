
#include "RooMvn2d.h"

#include "RooRandom.h"
#include "RooMath.h"
#include "RooHelpers.h"
#include "RooBatchCompute.h"


ClassImp(RooMvn2d);

////////////////////////////////////////////////////////////////////////////////

RooMvn2d::RooMvn2d(const char *name, const char *title,
          RooAbsReal& _x, RooAbsReal& _y,
          RooAbsReal& _mean_x, RooAbsReal& _mean_y,
          RooAbsReal& _sigma2_x, RooAbsReal& _sigma2_y,
          RooAbsReal& _cov_xy) :
  RooAbsPdf(name,title),
  x("x","Observable",this,_x),
  y("y","Observable",this,_y),
  mean_x("mean_x","Mean x",this,_mean_x),
  mean_y("mean_y","Mean y",this,_mean_y),
  sigma2_x("sigma2_x","Sigma2 x",this,_sigma2_x),
  sigma2_y("sigma2_y","Sigma2 y",this,_sigma2_y),
  cov_xy("cov_xy","cov(x,y)",this,_cov_xy)
{
}

////////////////////////////////////////////////////////////////////////////////

RooMvn2d::RooMvn2d(const RooMvn2d& other, const char* name) :
  RooAbsPdf(other,name),
  x("x",this,other.x),
  y("y",this,other.y),
  mean_x("mean_x",this,other.mean_x),
  mean_y("mean_y",this,other.mean_y),
  sigma2_x("sigma2_x",this,other.sigma2_x),
  sigma2_y("sigma2_y",this,other.sigma2_y),
  cov_xy("cov_xy",this,other.cov_xy)
{
}

////////////////////////////////////////////////////////////////////////////////

Double_t RooMvn2d::evaluate() const
{
  double det = sigma2_x * sigma2_y - cov_xy * cov_xy ;
  if ( det <= 0 ) det = 1. ; // what else should I do?
  double norm = 1./((2*3.14159265)*sqrt(det)) ;
  return norm * std::exp( -0.5 * (1./det)*(sigma2_y*(x-mean_x)*(x-mean_x) + sigma2_x*(y-mean_y)*(y-mean_y) - 2.*cov_xy*(x-mean_x)*(y-mean_y) ) ) ;
}


