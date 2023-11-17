

#include "corgen.h"

  //============================================================

  corgen::corgen( TMatrixD& covMat, Int_t rnSeed, Bool_t verb ) {

     _verb = verb ;

     _rand = new TRandom3( rnSeed ) ;

     Int_t nCols = covMat.GetNcols() ;
     Int_t nRows = covMat.GetNrows() ;

     if ( nCols != nRows ) {
        cout << " **** ERROR : corgen : nCols != nRows "
             << nCols << "," << nRows << endl ;
        return ;
     }

     _dim = nRows ;

     _covMat = new TMatrixD( nRows, nCols ) ;

     for ( Int_t i=0 ; i<nRows; i++ ) {
        for ( Int_t j=0 ; j<nCols; j++ ) {
           (*_covMat)(i,j) = covMat(i,j) ;
        } // j
     } // i

     cout << endl << endl << "  corgen constructor :" << endl ;
     cout << "    _covMat set to " << endl ;
     _covMat->Print() ;

     TVectorD eVals( nRows ) ;
     TMatrixD eVecMat = covMat.EigenVectors( eVals ) ;

     _eVals   = new TVectorD( eVals.GetNrows() ) ;
     _eVecMat = new TMatrixD( eVecMat.GetNrows(), eVecMat.GetNcols() ) ;

     for ( Int_t i=0 ; i<nRows; i++ ) {
        (*_eVals)(i) = eVals(i) ;
     } // i

     for ( Int_t i=0 ; i<nRows; i++ ) {
        for ( Int_t j=0 ; j<nCols; j++ ) {
           (*_eVecMat)(i,j) = eVecMat(i,j) ;
        } // j
     } // i

     cout << "    _eVals set to " << endl ;
     _eVals->Print() ;
     cout << "    _eVecMat set to " << endl ;
     _eVecMat->Print() ;

     cout << endl << endl ;

     //////// TMatrixD eVecMatInverse( nRows, nCols ) ;
     //////// eVecMatInverse = eVecMat ;
     //////// eVecMatInverse.Invert() ;


  } // constructor

  //============================================================

  corgen::~corgen() {
     delete _covMat ;
     delete _eVecMat ;
     delete _eVals ;
     delete _rand ;
  } // destructor

  //============================================================

  void corgen::genSet( TVectorD& rnVec ) {

     if ( rnVec.GetNrows() != _dim ) {
        cout << " *** ERROR : corgen::genSet : rnVec dim = " << rnVec.GetNrows()
             << ".  Was expecting " << _dim << endl ;
        return ;
     }

     for ( Int_t i=0; i<_dim; i++ ) {
        rnVec(i) = _rand->Gaus(  0. ,   sqrt( (*_eVals)(i) )   ) ;
     } // i


     rnVec *= (*_eVecMat) ;


  } // genSet

  //============================================================




