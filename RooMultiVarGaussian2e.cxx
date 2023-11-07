/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id$
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

/**
\file RooMultiVarGaussian2e.cxx
\class RooMultiVarGaussian2e
\ingroup Roofitcore

Multivariate Gaussian p.d.f. with correlations
**/

#include "RooFit.h"

#include "Riostream.h"
#include <math.h>

#include "RooMultiVarGaussian2e.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooRandom.h"
#include "RooMath.h"
#include "RooGlobalFunc.h"
#include "RooConstVar.h"
#include "TDecompChol.h"
#include "RooFitResult.h"

using namespace std;

ClassImp(RooMultiVarGaussian2e);
  ;

////////////////////////////////////////////////////////////////////////////////

RooMultiVarGaussian2e::RooMultiVarGaussian2e(const char *name, const char *title,
					 const RooArgList& xvec, const RooArgList& mu, const RooArgList& covRAL ) :
  RooAbsPdf(name,title),
  _x("x","Observables",this,kTRUE,kFALSE),
  _mu("mu","Offset vector",this,kTRUE,kFALSE),
  _covRAL("covRAL","cov mat",this,kTRUE,kFALSE),
  _z(4)
{
 _x.add(xvec) ;

 _mu.add(mu) ;

 _covRAL.add(covRAL) ;

 _det = _cov.Determinant() ;

 _ncall_since_last_update = 0 ;

  _prevCovVals = new double[ _covRAL.size() ] ;
  for ( int i=0; i<_covRAL.size(); i++ ) {
     _prevCovVals[i] = ((RooAbsReal*)_covRAL.at( i )) -> getVal() ;
     printf("  Setting _prevCovVals[%d] to %f\n", i, _prevCovVals[i] ) ;
  } // i


  _ndim = _x.size() ;
  _cov.ResizeTo(_ndim,_ndim) ;
  _covind.ResizeTo(_ndim,_ndim) ;
  int covralind(0) ;
   for ( int i=0; i<_ndim; i++ ) {
      for ( int j=0; j<_ndim; j++ ) {
         if ( j < i ) continue ;
         printf(" i,j = %d,%d\n", i, j ) ;
         _cov(i,j) = ((RooAbsReal*)_covRAL.at( covralind )) -> getVal() ;
         _cov(j,i) = ((RooAbsReal*)_covRAL.at( covralind )) -> getVal() ;
         _covind(i,j) = covralind ;
         _covind(j,i) = covralind ;
         covralind++ ;
      }
   }
  _cov.Print() ;
  _covind.Print() ;

 _covI.ResizeTo(_ndim,_ndim) ;
 _covI = _cov ;

 // Invert covariance matrix
 _covI.Invert() ;

 printf("\n\n cov inverse.\n") ;
 _covI.Print() ;

}


////////////////////////////////////////////////////////////////////////////////

RooMultiVarGaussian2e::RooMultiVarGaussian2e(const RooMultiVarGaussian2e& other, const char* name) : 
  RooAbsPdf(other,name), _x("x",this,other._x), _mu("mu",this,other._mu), _covRAL(other._covRAL),
  _cov(other._cov), _covI(other._covI), _det(other._det), _z(other._z),
  _covind(other._covind),
  _ndim(other._ndim)
{
  _prevCovVals = new double[ other._covRAL.size() ] ;
  for ( int i=0; i<other._covRAL.size(); i++ ) {
     _prevCovVals[i] = ((RooAbsReal*)other._covRAL.at( i )) -> getVal() ;
     printf("  copy constructor, Setting _prevCovVals[%d] to %f\n", i, _prevCovVals[i] ) ;
  } // i
  _ncall_since_last_update = 0 ;
}



////////////////////////////////////////////////////////////////////////////////

void RooMultiVarGaussian2e::syncMuVec() const 
{
  _muVec.ResizeTo(_mu.getSize()) ;
  for (Int_t i=0 ; i<_mu.getSize() ; i++) {
    _muVec[i] = ((RooAbsReal*)_mu.at(i))->getVal() ;
  }
}

////////////////////////////////////////////////////////////////////////////////

void RooMultiVarGaussian2e::syncCovMat() const
{
   bool needToRecalc(false) ;
   for ( int i=0; i<_covRAL.size(); i++ ) {
     double current_val = ((RooAbsReal*)_covRAL.at( i )) -> getVal() ;
     if ( fabs(current_val - _prevCovVals[i]) > 1e-8 ) {
        needToRecalc = true ;
        ///printf(" syncCovMat : Need to recalculate.  Par %d changed from %f to %f.  %d calls since last recalc.\n", i, _prevCovVals[i], current_val, _ncall_since_last_update ) ;
        break ;
     }
   } // i

   if ( needToRecalc ) {

      int cmralind(0) ;
      for ( int i=0; i<_ndim; i++ ) {
         for ( int j=0; j<_ndim; j++ ) {
            if ( j < i ) continue ;
            double valij = ((RooAbsReal*)_covRAL.at( TMath::Nint(_covind(i,j)) )) -> getVal() ;
            double valii = ((RooAbsReal*)_covRAL.at( TMath::Nint(_covind(i,i)) )) -> getVal() ;
            double valjj = ((RooAbsReal*)_covRAL.at( TMath::Nint(_covind(j,j)) )) -> getVal() ;
            double lim = sqrt(valii*valjj) ;
            if ( fabs(valij) > lim ) {
               if ( valij > 0 ) {
                  valij = 0.99 * lim ;
               } else {
                  valij = -0.99 * lim ;
               }
               ///////printf("*** warning\n") ;
            }
            _cov(i,j) = valij  ;
            if ( i != j ) _cov(j,i) = _cov(i,j) ;
         }
      }

      _det = _cov.Determinant() ;
      /////////if ( _det <= 0 ) { printf("\n\n *** bad determinant  %f\n\n", _det ) ; _cov.Print() ; }
      if ( _det <= 0 ) { return ; }
      _covI = _cov ;
      _covI.Invert() ;
      _ncall_since_last_update = 0 ;

      for ( int i=0; i<_covRAL.size(); i++ ) {
         _prevCovVals[i] = ((RooAbsReal*)_covRAL.at( i )) -> getVal() ;
      } // i

   } else {
      _ncall_since_last_update ++ ;
   }

}

////////////////////////////////////////////////////////////////////////////////

Double_t RooMultiVarGaussian2e::evaluate() const
{

  syncCovMat() ;

  TVectorD x(_mu.getSize()) ;
  for (int i=0 ; i<_x.getSize() ; i++) {
     x[i] =  ((RooAbsReal*)_x.at(i))->getVal() - ((RooAbsReal*)_mu.at(i))->getVal() ;
  }

  Double_t alpha = x * ( _covI * x ) ;

//double norm2 = pow( 2.*3.14159265, _ndim ) * _det ;
//double norm(1.) ;
//if ( norm2 > 0 ) norm = sqrt( norm2 ) ;
  double norm(1.) ;
  double rv = exp( -0.5 * alpha ) / norm ;

  return rv ;



}

////////////////////////////////////////////////////////////////////////////////



Int_t RooMultiVarGaussian2e::getAnalyticalIntegral(RooArgSet& allVarsIn, RooArgSet& analVars, const char* rangeName) const 
{
  RooArgSet allVars(allVarsIn) ;
  
  // If allVars contains x_i it cannot contain mu_i
  for (Int_t i=0 ; i<_x.getSize() ; i++) {
    if (allVars.contains(*_x.at(i))) {
      allVars.remove(*_mu.at(i),kTRUE,kTRUE) ;
    }
  }
  
  
  // Analytical integral known over all observables
  if (allVars.getSize()==_x.getSize() && !rangeName) {
    analVars.add(allVars) ;
    return -1 ;
  }

  Int_t code(0) ;

  Int_t nx = _x.getSize() ;
  if (nx>127) {
    // Warn that analytical integration is only provided for the first 127 observables
    coutW(Integration) << "RooMultiVarGaussian2e::getAnalyticalIntegral(" << GetName() << ") WARNING: p.d.f. has " << _x.getSize() 
		       << " observables, analytical integration is only implemented for the first 127 observables" << endl ;
    nx=127 ;
  }

  // Advertise partial analytical integral over all observables for which is wide enough to
  // use asymptotic integral calculation
  BitBlock bits ;
  Bool_t anyBits(kFALSE) ;
  syncMuVec() ;
  for (int i=0 ; i<_x.getSize() ; i++) {

    // Check if integration over observable #i is requested
    if (allVars.find(_x.at(i)->GetName())) {
      // Check if range is wider than Z sigma 
      RooRealVar* xi = (RooRealVar*)_x.at(i) ;
      if (xi->getMin(rangeName)<_muVec(i)-_z*sqrt(_cov(i,i)) && xi->getMax(rangeName) > _muVec(i)+_z*sqrt(_cov(i,i))) {
	cxcoutD(Integration) << "RooMultiVarGaussian2e::getAnalyticalIntegral(" << GetName() 
			     << ") Advertising analytical integral over " << xi->GetName() << " as range is >" << _z << " sigma" << endl ;
	bits.setBit(i) ;
	anyBits = kTRUE ;
	analVars.add(*allVars.find(_x.at(i)->GetName())) ;
      } else {
	cxcoutD(Integration) << "RooMultiVarGaussian2e::getAnalyticalIntegral(" << GetName() << ") Range of " << xi->GetName() << " is <" 
			     << _z << " sigma, relying on numeric integral" << endl ;	
      }
    }

    // Check if integration over parameter #i is requested
    if (allVars.find(_mu.at(i)->GetName())) {
      // Check if range is wider than Z sigma 
      RooRealVar* pi = (RooRealVar*)_mu.at(i) ;
      if (pi->getMin(rangeName)<_muVec(i)-_z*sqrt(_cov(i,i)) && pi->getMax(rangeName) > _muVec(i)+_z*sqrt(_cov(i,i))) {
	cxcoutD(Integration) << "RooMultiVarGaussian2e::getAnalyticalIntegral(" << GetName() 
			     << ") Advertising analytical integral over " << pi->GetName() << " as range is >" << _z << " sigma" << endl ;
	bits.setBit(i) ;
	anyBits = kTRUE ;
	analVars.add(*allVars.find(_mu.at(i)->GetName())) ;
      } else {
	cxcoutD(Integration) << "RooMultiVarGaussian2e::getAnalyticalIntegral(" << GetName() << ") Range of " << pi->GetName() << " is <" 
			     << _z << " sigma, relying on numeric integral" << endl ;	
      }
    }


  }

  // Full numeric integration over requested observables maps always to code zero
  if (!anyBits) {
    return 0 ;
  }
  
  // Map BitBlock into return code
  for (UInt_t i=0 ; i<_aicMap.size() ; i++) {
    if (_aicMap[i]==bits) {
      code = i+1 ;
    }
  }
  if (code==0) {
    _aicMap.push_back(bits) ;
    code = _aicMap.size() ;
  }

  return code ;
}



////////////////////////////////////////////////////////////////////////////////
/// Handle full integral here

Double_t RooMultiVarGaussian2e::analyticalIntegral(Int_t code, const char* /*rangeName*/) const 
{
  if (code==-1) {
    return pow(2*3.14159268,_x.getSize()/2.)*sqrt(fabs(_det)) ;
  }

  // Handle partial integrals here

  // Retrieve |S22|, S22bar from cache
  AnaIntData& aid = anaIntData(code) ;
 
  // Fill position vector for non-integrated observables
  syncMuVec() ;
  TVectorD u(aid.pmap.size()) ;
  for (UInt_t i=0 ; i<aid.pmap.size() ; i++) {
    u(i) = ((RooAbsReal*)_x.at(aid.pmap[i]))->getVal() - _muVec(aid.pmap[i]) ;
  }

  // Calculate partial integral
  Double_t ret = pow(2*3.14159268,aid.nint/2.)/sqrt(fabs(aid.S22det))*exp(-0.5*u*(aid.S22bar*u)) ;

  return ret ;
}



////////////////////////////////////////////////////////////////////////////////
/// Check if cache entry was previously created

RooMultiVarGaussian2e::AnaIntData& RooMultiVarGaussian2e::anaIntData(Int_t code) const 
{
  map<int,AnaIntData>::iterator iter =  _anaIntCache.find(code) ;
  if (iter != _anaIntCache.end()) {
    return iter->second ;
  }

  // Calculate cache contents  

  // Decode integration code
  vector<int> map1,map2 ;
  decodeCode(code,map1,map2) ;
  
  // Rearrage observables so that all non-integrated observables
  // go first (preserving relative order) and all integrated observables
  // go last (preserving relative order)
  TMatrixDSym S11, S22 ;
  TMatrixD S12, S21 ;
  blockDecompose(_covI,map1,map2,S11,S12,S21,S22) ;

  // Begin calculation of partial integrals
  //                                          ___
  //      sqrt(2pi)^(#intObs)     (-0.5 * u1T S22 u1 )
  // I =  ------------------- * e
  //        sqrt(|det(S22)|)
  //                                                                        ___
  // Where S22 is the sub-matrix of covI for the integrated observables and S22
  // is the Schur complement of S22
  // ___                   -1
  // S22  = S11 - S12 * S22   * S21
  //
  // and u1 is the vector of non-integrated observables

  // Calculate Schur complement S22bar
  TMatrixD S22inv(S22) ;
  S22inv.Invert() ;  
  TMatrixD S22bar = S11 - S12*S22inv*S21 ;  

  // Create new cache entry
  AnaIntData& cacheData = _anaIntCache[code] ;
  cacheData.S22bar.ResizeTo(S22bar) ;
  cacheData.S22bar=S22bar ;
  cacheData.S22det= S22.Determinant() ;
  cacheData.pmap = map1  ;
  cacheData.nint = map2.size() ;

  return cacheData ;
}



////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
/// Decode analytical integration/generation code into index map of integrated/generated (map2)
/// and non-integrated/generated observables (map1)

void RooMultiVarGaussian2e::decodeCode(Int_t code, vector<int>& map1, vector<int>& map2) const
{
  if (code<0 || code> (Int_t)_aicMap.size()) {
    cout << "RooMultiVarGaussian2e::decodeCode(" << GetName() << ") ERROR don't have bit pattern for code " << code << endl ;
    throw string("RooMultiVarGaussian2e::decodeCode() ERROR don't have bit pattern for code") ;
  }

  BitBlock b = _aicMap[code-1] ;  
  map1.clear() ;
  map2.clear() ;
  for (int i=0 ; i<_x.getSize() ; i++) {
    if (b.getBit(i)) {
      map2.push_back(i) ;
    } else {
      map1.push_back(i) ;
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Block decomposition of covI according to given maps of observables

void RooMultiVarGaussian2e::blockDecompose(const TMatrixD& input, const vector<int>& map1, const vector<int>& map2, TMatrixDSym& S11, TMatrixD& S12, TMatrixD& S21, TMatrixDSym& S22)
{
  // Allocate and fill reordered covI matrix in 2x2 block structure

  S11.ResizeTo(map1.size(),map1.size()) ; 
  S12.ResizeTo(map1.size(),map2.size()) ;
  S21.ResizeTo(map2.size(),map1.size()) ;
  S22.ResizeTo(map2.size(),map2.size()) ;

  for (UInt_t i=0 ; i<map1.size() ; i++) {
    for (UInt_t j=0 ; j<map1.size() ; j++) 
      S11(i,j) = input(map1[i],map1[j]) ;
    for (UInt_t j=0 ; j<map2.size() ; j++) 
      S12(i,j) = input(map1[i],map2[j]) ;
  }
  for (UInt_t i=0 ; i<map2.size() ; i++) {
    for (UInt_t j=0 ; j<map1.size() ; j++) 
      S21(i,j) = input(map2[i],map1[j]) ;
    for (UInt_t j=0 ; j<map2.size() ; j++) 
      S22(i,j) = input(map2[i],map2[j]) ;
  }
  
}


void RooMultiVarGaussian2e::BitBlock::setBit(Int_t ibit) 
{
  if (ibit<32) { b0 |= (1<<ibit) ; return ; }
  if (ibit<64) { b1 |= (1<<(ibit-32)) ; return ; }
  if (ibit<96) { b2 |= (1<<(ibit-64)) ; return ; }
  if (ibit<128) { b3 |= (1<<(ibit-96)) ; return ; }
}

Bool_t RooMultiVarGaussian2e::BitBlock::getBit(Int_t ibit) 
{
  if (ibit<32) return (b0 & (1<<ibit)) ; 
  if (ibit<64) return (b1 & (1<<(ibit-32))) ; 
  if (ibit<96) return (b2 & (1<<(ibit-64))) ; 
  if (ibit<128) return (b3 & (1<<(ibit-96))) ; 
  return kFALSE ;
}

Bool_t RooMultiVarGaussian2e::BitBlock::operator==(const BitBlock& other) 
{
  if (b0 != other.b0) return kFALSE ;
  if (b1 != other.b1) return kFALSE ;
  if (b2 != other.b2) return kFALSE ;
  if (b3 != other.b3) return kFALSE ;
  return kTRUE ;
}







