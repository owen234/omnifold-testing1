//
//   Owen Long
//
//   U. C. Riverside
//
//

#ifndef corgen_h
#define corgen_h

#include <TRandom3.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <iostream>

  using std::cout ;
  using std::endl ;

   class corgen  {

       public:

          corgen( TMatrixD& covMat, Int_t rnSeed=12345, Bool_t verb=kFALSE ) ;

          virtual ~corgen() ;

          void genSet( TVectorD& rnVec ) ;

       private:

          TMatrixD* _covMat ;
          TMatrixD* _eVecMat ;
          TVectorD* _eVals ;

          TRandom3* _rand ;

          Bool_t    _verb ;
          Int_t     _dim ;


   } ;

#endif // ifndef for corgen_h

