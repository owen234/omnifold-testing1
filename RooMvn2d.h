#ifndef ROO_MVN2D
#define ROO_MVN2D

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooAbsReal;

class RooMvn2d : public RooAbsPdf {
public:
  RooMvn2d() { };
  RooMvn2d(const char *name, const char *title,
          RooAbsReal& _x, RooAbsReal& _y,
          RooAbsReal& _mean_x, RooAbsReal& _mean_y,
          RooAbsReal& _sigma2_x, RooAbsReal& _sigma2_y,
          RooAbsReal& _cov_xy
          );
  RooMvn2d(const RooMvn2d& other, const char* name=0);
  virtual TObject* clone(const char* newname) const override {
    return new RooMvn2d(*this,newname);
  }
  inline virtual ~RooMvn2d() { }


protected:

  RooRealProxy x ;
  RooRealProxy y ;
  RooRealProxy mean_x ;
  RooRealProxy mean_y ;
  RooRealProxy sigma2_x ;
  RooRealProxy sigma2_y ;
  RooRealProxy cov_xy ;

  Double_t evaluate() const override;

private:

  ClassDefOverride(RooMvn2d,1)
};

#endif
