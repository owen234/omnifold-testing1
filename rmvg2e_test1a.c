
#include "RooMultiVarGaussian2e.h"

using namespace RooFit;



void rmvg2e_test1a( int dim=2 ) {


   TMatrixDSym covSetvals( dim ) ;

   covSetvals(0,0) = 1.1 ;
   covSetvals(1,1) = 1.3 ;
   covSetvals(0,1) = -0.6 ;
   covSetvals(1,0) = covSetvals(0,1) ;
   if ( dim > 2 ) {
      covSetvals(2,2) = 0.7 ;
      covSetvals(0,2) = -0.3 ;
      covSetvals(1,2) =  0.4 ;
      covSetvals(2,0) = covSetvals(0,2) ;
      covSetvals(2,1) = covSetvals(1,2) ;
   }
   if ( dim > 3 ) {
      covSetvals(3,3) = 0.3 ;
      covSetvals(0,3) = 0.0 ;
      covSetvals(1,3) = 0.0 ;
      covSetvals(2,3) = 0.0 ;
      covSetvals(3,0) = covSetvals(0,3) ;
      covSetvals(3,1) = covSetvals(1,3) ;
      covSetvals(3,2) = covSetvals(2,3) ;
   }
   if ( dim > 4 ) {
      covSetvals(4,4) = 1.3 ;
      covSetvals(0,4) = 0.0 ;
      covSetvals(1,4) = 0.0 ;
      covSetvals(2,4) = 0.0 ;
      covSetvals(3,4) = 0.0 ;
      covSetvals(4,0) = covSetvals(0,4) ;
      covSetvals(4,1) = covSetvals(1,4) ;
      covSetvals(4,2) = covSetvals(2,4) ;
      covSetvals(4,3) = covSetvals(3,4) ;
   }


   TVectorD muSetvals( dim ) ;
   muSetvals(0) = 0.2 ;
   muSetvals(1) = -0.3 ;
   if ( dim > 2 ) {
      muSetvals(2) =  0.5 ;
   }
   if ( dim > 3 ) {
      muSetvals(3) =  0.0 ;
   }
   if ( dim > 4 ) {
      muSetvals(4) = -0.1 ;
   }


   RooArgList xVec;
   RooArgList muVec;

   // make the observable and means
   Int_t i, j;
   RooRealVar *x;
   RooRealVar *mu_x;
   for (i = 0; i < dim; i++) {
      char *name = Form("x%d", i);
      x = new RooRealVar(name, name, 0, -5, 5);
      xVec.add(*x);

      char *mu_name = Form("mu_x%d", i);
      //mu_x = new RooRealVar(mu_name, mu_name, 0, -2, 2);
      mu_x = new RooRealVar(mu_name, mu_name, muSetvals(i), -2, 2);
      muVec.add(*mu_x);
   }
   RooArgSet obs_ras(xVec) ;

   // make a covariance matrix that is all 1's









   RooArgList covRAL ;
   for (i=0; i<dim; i++ ) {
      for (j=0; j<dim; j++ ) {
         if ( j<i ) continue ;
         RooRealVar* cov_rrv(0x0) ;
         if ( i == j ) {
            //cov_rrv = new RooRealVar( Form("cov%d%d", i, j), Form("cov%d%d", i, j), 1.0, 0.1, 9. )  ;
            cov_rrv = new RooRealVar( Form("cov%d%d", i, j), Form("cov%d%d", i, j), covSetvals(i,j), 0.1, 9. )  ;
         } else {
            //cov_rrv = new RooRealVar( Form("cov%d%d", i, j), Form("cov%d%d", i, j), 0.2, -9., 9. )  ;
            cov_rrv = new RooRealVar( Form("cov%d%d", i, j), Form("cov%d%d", i, j), covSetvals(i,j), -9., 9. )  ;
         }
         /////////if ( i!=j ) cov_rrv -> setConstant() ; // try fixing some parameters
         covRAL.add( *cov_rrv ) ;
      } // j
   } // i





   // now make the multivariate Gaussian
   RooMultiVarGaussian2e mvg("mvg", "mvg", xVec, muVec, covRAL);



   //RooDataSet* data = mvg.generate( obs_ras, 10000 ) ;
   RooDataSet* data = mvg.generate( obs_ras, 10000 ) ;

    data->Print("V") ;
    data->get(0) -> Print("V") ;
    data->get(1) -> Print("V") ;
    data->get(2) -> Print("V") ;
    data->get(3) -> Print("V") ;
    data->get(4) -> Print("V") ;

    //RooFitResult* rfr = mvg.fitTo( *data, Save(true) ) ;
    RooFitResult* rfr = mvg.fitTo( *data, NumCPU(20), Save(true) ) ;



    TCanvas* can = new TCanvas("can","can", 50, 50, 1200, 800 ) ;
    can -> Divide(2,1) ;


    can -> cd(1) ;
    //RooPlot* xframe = x.frame(Title("data")) ;
    RooPlot* xframe = ((RooRealVar*)xVec.at(0))->frame(Title("x frame")) ;
    data -> plotOn(xframe) ;
    mvg.plotOn(xframe, NumCPU(20)) ;
    xframe -> Draw() ;

    can -> cd(2) ;
    RooPlot* yframe = ((RooRealVar*)xVec.at(1))->frame(Title("data")) ;
    data -> plotOn(yframe) ;
    mvg.plotOn(yframe, NumCPU(20)) ;
    yframe -> Draw() ;

    printf("\n\n cov setvals:\n" ) ;
    covSetvals.Print() ;

    printf("\n\n mu setvals:\n" ) ;
    muSetvals.Print() ;

    printf("\n\n\n") ;
    int list_ind(0) ;
    for (i=0; i<dim; i++ ) {
       for (j=0; j<dim; j++ ) {
          if ( j<i ) continue ;
          double fitval = ((RooRealVar*)covRAL.at(list_ind))->getVal() ;
          double fiterr = ((RooRealVar*)covRAL.at(list_ind))->getError() ;
          double genval = covSetvals(i,j) ;
          char parname[100] ;
          sprintf( parname, "%s", ((RooRealVar*)covRAL.at(list_ind))->GetName() ) ;
          printf("  %8s :  true  %8.4f,   fit  %8.4f +/- %8.4f     diff %8.4f     diff/err %7.2f\n", parname, genval, fitval, fiterr, (fitval-genval), (fitval-genval)/fiterr ) ;
          list_ind++ ;
       } // j
    } // i
    printf("\n\n\n") ;
    for (i=0; i<dim; i++ ) {
       double fitval = ((RooRealVar*)muVec.at(i))->getVal() ;
       double fiterr = ((RooRealVar*)muVec.at(i))->getError() ;
       double genval = muSetvals(i) ;
       char parname[100] ;
       sprintf( parname, "%s", ((RooRealVar*)muVec.at(i))->GetName() ) ;
       printf("  %8s :  true  %8.4f,   fit  %8.4f +/- %8.4f     diff %8.4f     diff/err %7.2f\n", parname, genval, fitval, fiterr, (fitval-genval), (fitval-genval)/fiterr ) ;
    } // i


    printf("\n\n\n") ;
}
