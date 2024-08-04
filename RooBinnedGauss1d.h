#ifndef ROO_MVN2D
#define ROO_MVN2D

#include "RooAbsPdf.h"
#include "RooConstVar.h"
#include "RooRealProxy.h"

class RooAbsReal;

class RooBinnedGauss1d : public RooAbsPdf {
public:
  RooBinnedGauss1d() { };
  RooBinnedGauss1d(const char *name, const char *title,
          RooAbsReal& _x,
          RooAbsReal& _mean_x,
          RooAbsReal& _sigma2_x,
          int nbins,
          const RooArgList& bin_edges
          );
  RooBinnedGauss1d(const RooBinnedGauss1d& other, const char* name=0);
  virtual TObject* clone(const char* newname) const override {
    return new RooBinnedGauss1d(*this,newname);
  }
  inline virtual ~RooBinnedGauss1d() { }


protected:

  RooRealProxy x ;
  RooRealProxy mean_x ;
  RooRealProxy sigma2_x ;

  Double_t evaluate() const override;

  int nbins_ ;
  float* bin_edges_ ;

private:

  ClassDefOverride(RooBinnedGauss1d,1)
};

#endif
