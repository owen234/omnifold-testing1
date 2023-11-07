/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id$
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_MULTI_VAR_GAUSSIAN2e
#define ROO_MULTI_VAR_GAUSSIAN2e

#include "RooAbsPdf.h"
#include "RooListProxy.h"
#include "TMatrixDSym.h"
#include "TMatrixD.h"
#include "TMatrixF.h"
#include "TVectorD.h"

class RooRealVar;
class RooFitResult ;

#include <map>
#include <vector>

class RooMultiVarGaussian2e : public RooAbsPdf {
public:

  RooMultiVarGaussian2e() {} ;
  //RooMultiVarGaussian2e(const char *name, const char *title, const RooArgList& xvec, const RooArgList& mu, const TMatrixDSym& covMatrix) ;
  RooMultiVarGaussian2e(const char *name, const char *title, const RooArgList& xvec, const RooArgList& mu, const RooArgList& covMatrixRAL ) ;

  RooMultiVarGaussian2e(const RooMultiVarGaussian2e& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooMultiVarGaussian2e(*this,newname); }
  inline virtual ~RooMultiVarGaussian2e() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ; 
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ; 

  class AnaIntData {
  public:
    TMatrixD    S22bar ;
    Double_t    S22det ;
    std::vector<int> pmap ;
    Int_t       nint ;
  } ;

  class GenData {
  public:
    TMatrixD    UT ;
    std::vector<int> omap ;
    std::vector<int> pmap ;
    TVectorD    mu1 ;
    TVectorD    mu2 ;
    TMatrixD    S12S22I ;
  } ;

  class BitBlock {
  public:
    BitBlock() : b0(0), b1(0), b2(0), b3(0) {} ;

    void setBit(Int_t ibit) ;      
    Bool_t getBit(Int_t ibit) ;
    Bool_t operator==(const BitBlock& other) ;

    Int_t b0 ;
    Int_t b1 ;
    Int_t b2 ;
    Int_t b3 ;
  } ;

  static void blockDecompose(const TMatrixD& input, const std::vector<int>& map1, const std::vector<int>& map2, TMatrixDSym& S11, TMatrixD& S12, TMatrixD& S21, TMatrixDSym& S22) ;

protected:

  RooListProxy _x ;
  RooListProxy _mu ;
  RooListProxy _covRAL ;
  mutable TMatrixDSym _cov ;
  mutable TMatrixDSym _covI ;
  mutable Double_t    _det ; 
  Double_t    _z ;
  Int_t       _ndim ;
  mutable Int_t _ncall_since_last_update ;

  TMatrixF _covind ;

  void syncCovMat() const ;

  double* _prevCovVals ;

  void syncMuVec() const ;
  mutable TVectorD _muVec ; //! Do not persist

  Double_t evaluate() const ;

  void decodeCode(Int_t code, std::vector<int>& map1, std::vector<int>& map2) const;
  AnaIntData& anaIntData(Int_t code) const ;

  mutable std::map<int,AnaIntData> _anaIntCache ; //!
  mutable std::map<int,GenData> _genCache ; //!

  mutable std::vector<BitBlock> _aicMap ; //!


private:

  ClassDef(RooMultiVarGaussian2e,1) // Multivariate Gaussian PDF with correlations
};

#endif
