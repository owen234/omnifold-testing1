
#include "RooMultiVarGaussian2e.h"

using namespace RooFit;



void rmvg2e_fit_bootstrap1( int bsi=0, const char* infile = "bootstrap-9d-export1g-diagonal-cov.root" ) {

   int dim=9 ;

    TChain ch("boot_ttree") ;
    ch.Add( infile ) ;

    float x0 ;
    float x1 ;
    float x2 ;
    float x3 ;
    float x4 ;
    float x5 ;
    float x6 ;
    float x7 ;
    float x8 ;
    float tt_weight ;
    int   dset_index ;
    ch.SetBranchAddress("x0",&x0) ;
    ch.SetBranchAddress("x1",&x1) ;
    ch.SetBranchAddress("x2",&x2) ;
    ch.SetBranchAddress("x3",&x3) ;
    ch.SetBranchAddress("x4",&x4) ;
    ch.SetBranchAddress("x5",&x5) ;
    ch.SetBranchAddress("x6",&x6) ;
    ch.SetBranchAddress("x7",&x7) ;
    ch.SetBranchAddress("x8",&x8) ;
    ch.SetBranchAddress("weight",&tt_weight) ;
    ch.SetBranchAddress("dset_index",&dset_index) ;


// TMatrixDSym covSetvals( dim ) ;

// covSetvals(0,0) = 1.1 ;
// covSetvals(1,1) = 1.3 ;
// covSetvals(0,1) = -0.6 ;
// covSetvals(1,0) = covSetvals(0,1) ;
// if ( dim > 2 ) {
//    covSetvals(2,2) = 0.7 ;
//    covSetvals(0,2) = -0.3 ;
//    covSetvals(1,2) =  0.4 ;
//    covSetvals(2,0) = covSetvals(0,2) ;
//    covSetvals(2,1) = covSetvals(1,2) ;
// }
// if ( dim > 3 ) {
//    covSetvals(3,3) = 0.3 ;
//    covSetvals(0,3) = 0.0 ;
//    covSetvals(1,3) = 0.0 ;
//    covSetvals(2,3) = 0.0 ;
//    covSetvals(3,0) = covSetvals(0,3) ;
//    covSetvals(3,1) = covSetvals(1,3) ;
//    covSetvals(3,2) = covSetvals(2,3) ;
// }
// if ( dim > 4 ) {
//    covSetvals(4,4) = 1.3 ;
//    covSetvals(0,4) = 0.0 ;
//    covSetvals(1,4) = 0.0 ;
//    covSetvals(2,4) = 0.0 ;
//    covSetvals(3,4) = 0.0 ;
//    covSetvals(4,0) = covSetvals(0,4) ;
//    covSetvals(4,1) = covSetvals(1,4) ;
//    covSetvals(4,2) = covSetvals(2,4) ;
//    covSetvals(4,3) = covSetvals(3,4) ;
// }


// TVectorD muSetvals( dim ) ;
// muSetvals(0) = 0.2 ;
// muSetvals(1) = -0.3 ;
// if ( dim > 2 ) {
//    muSetvals(2) =  0.5 ;
// }
// if ( dim > 3 ) {
//    muSetvals(3) =  0.0 ;
// }
// if ( dim > 4 ) {
//    muSetvals(4) = -0.1 ;
// }


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
      //mu_x = new RooRealVar(mu_name, mu_name, muSetvals(i), -2, 2);
      mu_x = new RooRealVar(mu_name, mu_name, 0., -2, 2);
      muVec.add(*mu_x);
   }
   RooArgSet obs_ras(xVec) ;



    ((RooRealVar*)muVec.at(0))->setVal(0.8) ;
    ((RooRealVar*)muVec.at(1))->setVal(0.1) ;
    ((RooRealVar*)muVec.at(2))->setVal(-0.6) ;
    ((RooRealVar*)muVec.at(3))->setVal(0.3) ;
    ((RooRealVar*)muVec.at(4))->setVal(1.3) ;
    ((RooRealVar*)muVec.at(5))->setVal(0.2) ;
    ((RooRealVar*)muVec.at(6))->setVal(-0.5) ;
    ((RooRealVar*)muVec.at(7))->setVal(0.2) ;
    ((RooRealVar*)muVec.at(8))->setVal(0.6) ;









   RooArgList covRAL ;
   for (i=0; i<dim; i++ ) {
      for (j=0; j<dim; j++ ) {
         if ( j<i ) continue ;
         RooRealVar* cov_rrv(0x0) ;
         if ( i == j ) {
            /////////cov_rrv = new RooRealVar( Form("cov%d%d", i, j), Form("cov%d%d", i, j), covSetvals(i,j), 0.1, 9. )  ;
            cov_rrv = new RooRealVar( Form("cov%d%d", i, j), Form("cov%d%d", i, j), 0.6, 0.02, 2. )  ;
         } else {
            ////////cov_rrv = new RooRealVar( Form("cov%d%d", i, j), Form("cov%d%d", i, j), covSetvals(i,j), -9., 9. )  ;
            ////cov_rrv = new RooRealVar( Form("cov%d%d", i, j), Form("cov%d%d", i, j), 0., -9., 9. )  ;
            cov_rrv = new RooRealVar( Form("cov%d%d", i, j), Form("cov%d%d", i, j), 0., -1.1, 1.1 )  ;
         }
         if ( i!=j ) cov_rrv -> setConstant() ; // try fixing some parameters
         ////////if ( i==j ) cov_rrv -> setConstant() ; // try fixing some parameters
         covRAL.add( *cov_rrv ) ;
      } // j
   } // i


  //--- initialize to true values.

   ((RooRealVar*)covRAL.find("cov00"))->setVal(0.64) ;
   ((RooRealVar*)covRAL.find("cov11"))->setVal(0.36) ;
   ((RooRealVar*)covRAL.find("cov22"))->setVal(1.00) ;
   ((RooRealVar*)covRAL.find("cov33"))->setVal(0.36) ;
   ((RooRealVar*)covRAL.find("cov44"))->setVal(1.44) ;
   ((RooRealVar*)covRAL.find("cov55"))->setVal(0.09) ;
   ((RooRealVar*)covRAL.find("cov66"))->setVal(0.09) ;
   ((RooRealVar*)covRAL.find("cov77"))->setVal(0.25) ;
   ((RooRealVar*)covRAL.find("cov88"))->setVal(1.00) ;




// ((RooRealVar*)covRAL.find("cov01"))->setVal( 0.288) ;
// ((RooRealVar*)covRAL.find("cov02"))->setVal( 0.400) ;
// ((RooRealVar*)covRAL.find("cov03"))->setVal(-0.144) ;
// ((RooRealVar*)covRAL.find("cov04"))->setVal( 0.096) ;
// ((RooRealVar*)covRAL.find("cov05"))->setVal(-0.024) ;
// ((RooRealVar*)covRAL.find("cov06"))->setVal( 0.024) ;
// ((RooRealVar*)covRAL.find("cov07"))->setVal( 0.000) ;
// ((RooRealVar*)covRAL.find("cov08"))->setVal( 0.000) ;

// ((RooRealVar*)covRAL.find("cov12"))->setVal( 0.420) ;
// ((RooRealVar*)covRAL.find("cov13"))->setVal( 0.144) ;
// ((RooRealVar*)covRAL.find("cov14"))->setVal( 0.288) ;
// ((RooRealVar*)covRAL.find("cov15"))->setVal( 0.054) ;
// ((RooRealVar*)covRAL.find("cov16"))->setVal( 0.036) ;
// ((RooRealVar*)covRAL.find("cov17"))->setVal(-0.030) ;
// ((RooRealVar*)covRAL.find("cov18"))->setVal( 0.000) ;

// ((RooRealVar*)covRAL.find("cov23"))->setVal( 0.360) ;
// ((RooRealVar*)covRAL.find("cov24"))->setVal(-0.600) ;
// ((RooRealVar*)covRAL.find("cov25"))->setVal(-0.090) ;
// ((RooRealVar*)covRAL.find("cov26"))->setVal( 0.060) ;
// ((RooRealVar*)covRAL.find("cov27"))->setVal( 0.050) ;
// ((RooRealVar*)covRAL.find("cov28"))->setVal( 0.000) ;

// ((RooRealVar*)covRAL.find("cov34"))->setVal(-0.500) ;
// ((RooRealVar*)covRAL.find("cov35"))->setVal( 0.108) ;
// ((RooRealVar*)covRAL.find("cov36"))->setVal(-0.072) ;
// ((RooRealVar*)covRAL.find("cov37"))->setVal(-0.090) ;
// ((RooRealVar*)covRAL.find("cov38"))->setVal( 0.000) ;

// ((RooRealVar*)covRAL.find("cov45"))->setVal(-0.216) ;
// ((RooRealVar*)covRAL.find("cov46"))->setVal( 0.144) ;
// ((RooRealVar*)covRAL.find("cov47"))->setVal( 0.120) ;
// ((RooRealVar*)covRAL.find("cov48"))->setVal(-0.120) ;

// ((RooRealVar*)covRAL.find("cov56"))->setVal( 0.027) ;
// ((RooRealVar*)covRAL.find("cov57"))->setVal( 0.000) ;
// ((RooRealVar*)covRAL.find("cov58"))->setVal( 0.000) ;

// ((RooRealVar*)covRAL.find("cov67"))->setVal(-0.045) ;
// ((RooRealVar*)covRAL.find("cov68"))->setVal( 0.000) ;

// ((RooRealVar*)covRAL.find("cov78"))->setVal( 0.100) ;



   ((RooRealVar*)covRAL.find("cov01"))->setVal( 0.000) ;
   ((RooRealVar*)covRAL.find("cov02"))->setVal( 0.000) ;
   ((RooRealVar*)covRAL.find("cov03"))->setVal( 0.000) ;
   ((RooRealVar*)covRAL.find("cov04"))->setVal( 0.000) ;
   ((RooRealVar*)covRAL.find("cov05"))->setVal( 0.000) ;
   ((RooRealVar*)covRAL.find("cov06"))->setVal( 0.000) ;
   ((RooRealVar*)covRAL.find("cov07"))->setVal( 0.000) ;
   ((RooRealVar*)covRAL.find("cov08"))->setVal( 0.000) ;

   ((RooRealVar*)covRAL.find("cov12"))->setVal( 0.000) ;
   ((RooRealVar*)covRAL.find("cov13"))->setVal( 0.000) ;
   ((RooRealVar*)covRAL.find("cov14"))->setVal( 0.000) ;
   ((RooRealVar*)covRAL.find("cov15"))->setVal( 0.000) ;
   ((RooRealVar*)covRAL.find("cov16"))->setVal( 0.000) ;
   ((RooRealVar*)covRAL.find("cov17"))->setVal( 0.000) ;
   ((RooRealVar*)covRAL.find("cov18"))->setVal( 0.000) ;

   ((RooRealVar*)covRAL.find("cov23"))->setVal( 0.000) ;
   ((RooRealVar*)covRAL.find("cov24"))->setVal( 0.000) ;
   ((RooRealVar*)covRAL.find("cov25"))->setVal( 0.000) ;
   ((RooRealVar*)covRAL.find("cov26"))->setVal( 0.000) ;
   ((RooRealVar*)covRAL.find("cov27"))->setVal( 0.000) ;
   ((RooRealVar*)covRAL.find("cov28"))->setVal( 0.000) ;

   ((RooRealVar*)covRAL.find("cov34"))->setVal( 0.000) ;
   ((RooRealVar*)covRAL.find("cov35"))->setVal( 0.000) ;
   ((RooRealVar*)covRAL.find("cov36"))->setVal( 0.000) ;
   ((RooRealVar*)covRAL.find("cov37"))->setVal( 0.000) ;
   ((RooRealVar*)covRAL.find("cov38"))->setVal( 0.000) ;

   ((RooRealVar*)covRAL.find("cov45"))->setVal( 0.000) ;
   ((RooRealVar*)covRAL.find("cov46"))->setVal( 0.000) ;
   ((RooRealVar*)covRAL.find("cov47"))->setVal( 0.000) ;
   ((RooRealVar*)covRAL.find("cov48"))->setVal( 0.000) ;

   ((RooRealVar*)covRAL.find("cov56"))->setVal( 0.000) ;
   ((RooRealVar*)covRAL.find("cov57"))->setVal( 0.000) ;
   ((RooRealVar*)covRAL.find("cov58"))->setVal( 0.000) ;

   ((RooRealVar*)covRAL.find("cov67"))->setVal( 0.000) ;
   ((RooRealVar*)covRAL.find("cov68"))->setVal( 0.000) ;

   ((RooRealVar*)covRAL.find("cov78"))->setVal( 0.000) ;






    RooRealVar weight( "weight", "weight", 0., 100. ) ;






   // now make the multivariate Gaussian
   RooMultiVarGaussian2e mvg("mvg", "mvg", xVec, muVec, covRAL);


//RooDataSet* data = new RooDataSet( "ds", "ds", RooArgSet(xVec,weight), WeightVar("weight") ) ;

//  int nadd(0) ;
//  for ( int i=0; i<ch.GetEntries(); i++ ) {
//     ch.GetEntry(i) ;
//     if ( dset_index != bsi ) continue ;
//     ((RooRealVar*)xVec.at(0))->setVal(x0) ;
//     ((RooRealVar*)xVec.at(1))->setVal(x1) ;
//     ((RooRealVar*)xVec.at(2))->setVal(x2) ;
//     ((RooRealVar*)xVec.at(3))->setVal(x3) ;
//     ((RooRealVar*)xVec.at(4))->setVal(x4) ;
//     ((RooRealVar*)xVec.at(5))->setVal(x5) ;
//     ((RooRealVar*)xVec.at(6))->setVal(x6) ;
//     ((RooRealVar*)xVec.at(7))->setVal(x7) ;
//     ((RooRealVar*)xVec.at(8))->setVal(x8) ;
//     weight.setVal(tt_weight) ;
//     data -> add( RooArgSet( xVec, weight ), weight.getVal() ) ;
//     nadd++ ;
//  }
//  printf("\n\n Added %d events.\n\n", nadd ) ;

    RooDataSet* data = mvg.generate( obs_ras, 10000 ) ;

    data->Print("V") ;
    data->get(0) -> Print("V") ;
    data->get(1) -> Print("V") ;
    data->get(2) -> Print("V") ;
    data->get(3) -> Print("V") ;
    data->get(4) -> Print("V") ;


//  ((RooRealVar*)muVec.at(0))->setConstant() ;
//  ((RooRealVar*)muVec.at(1))->setConstant() ;
//  ((RooRealVar*)muVec.at(2))->setConstant() ;
//  ((RooRealVar*)muVec.at(3))->setConstant() ;
//  ((RooRealVar*)muVec.at(4))->setConstant() ;
//  ((RooRealVar*)muVec.at(5))->setConstant() ;
//  ((RooRealVar*)muVec.at(6))->setConstant() ;
//  ((RooRealVar*)muVec.at(7))->setConstant() ;
//  ((RooRealVar*)muVec.at(8))->setConstant() ;



   //-- first pass floating only mu parameters and diagonal cov elements.

    RooFitResult* rfr = mvg.fitTo( *data, NumCPU(28), Save(true) ) ;





//  printf("\n\n ============================================================================================================\n\n") ;
//  printf("               Second pass fit\n") ;
//  printf("\n\n ============================================================================================================\n\n") ;


    for (i=0; i<dim; i++ ) {
       for (j=0; j<dim; j++ ) {
          if ( j<i ) continue ;
          RooRealVar* cov_rrv = (RooRealVar*) covRAL.find( Form("cov%d%d", i, j) )  ;
          cov_rrv -> setConstant(false) ;
       } // j
    } // i

    RooFitResult* rfr2 = mvg.fitTo( *data, NumCPU(28), Save(true) ) ;



// //-- third pass with diagonal cov elements and mean pars floating, rest fixed

//  printf("\n\n ============================================================================================================\n\n") ;
//  printf("               Third pass fit\n") ;
//  printf("\n\n ============================================================================================================\n\n") ;


//  for (i=0; i<dim; i++ ) {
//     for (j=0; j<dim; j++ ) {
//        if ( j<i ) continue ;
//        RooRealVar* cov_rrv = (RooRealVar*) covRAL.find( Form("cov%d%d", i, j) )  ;
//        if ( i == j ) {
//           cov_rrv -> setConstant(false) ;
//        } else {
//           cov_rrv -> setConstant(true) ;
//        }
//     } // j
//  } // i

//  ((RooRealVar*)muVec.at(0))->setConstant(false) ;
//  ((RooRealVar*)muVec.at(1))->setConstant(false) ;
//  ((RooRealVar*)muVec.at(2))->setConstant(false) ;
//  ((RooRealVar*)muVec.at(3))->setConstant(false) ;
//  ((RooRealVar*)muVec.at(4))->setConstant(false) ;
//  ((RooRealVar*)muVec.at(5))->setConstant(false) ;
//  ((RooRealVar*)muVec.at(6))->setConstant(false) ;
//  ((RooRealVar*)muVec.at(7))->setConstant(false) ;
//  ((RooRealVar*)muVec.at(8))->setConstant(false) ;

//  RooFitResult* rfr3 = mvg.fitTo( *data, NumCPU(28), Save(true) ) ;



//  TCanvas* can = new TCanvas("can","can", 50, 50, 1200, 800 ) ;
//  can -> Divide(2,1) ;


//  can -> cd(1) ;
//  //RooPlot* xframe = x.frame(Title("data")) ;
//  RooPlot* xframe = ((RooRealVar*)xVec.at(0))->frame(Title("x frame")) ;
//  data -> plotOn(xframe) ;
//  mvg.plotOn(xframe, NumCPU(20)) ;
//  xframe -> Draw() ;

//  can -> cd(2) ;
//  RooPlot* yframe = ((RooRealVar*)xVec.at(1))->frame(Title("data")) ;
//  data -> plotOn(yframe) ;
//  mvg.plotOn(yframe, NumCPU(20)) ;
//  yframe -> Draw() ;

//  printf("\n\n cov setvals:\n" ) ;
//  covSetvals.Print() ;

//  printf("\n\n mu setvals:\n" ) ;
//  muSetvals.Print() ;

//  printf("\n\n\n") ;
//  int list_ind(0) ;
//  for (i=0; i<dim; i++ ) {
//     for (j=0; j<dim; j++ ) {
//        if ( j<i ) continue ;
//        double fitval = ((RooRealVar*)covRAL.at(list_ind))->getVal() ;
//        double fiterr = ((RooRealVar*)covRAL.at(list_ind))->getError() ;
//        double genval = covSetvals(i,j) ;
//        char parname[100] ;
//        sprintf( parname, "%s", ((RooRealVar*)covRAL.at(list_ind))->GetName() ) ;
//        printf("  %8s :  true  %8.4f,   fit  %8.4f +/- %8.4f     diff %8.4f     diff/err %7.2f\n", parname, genval, fitval, fiterr, (fitval-genval), (fitval-genval)/fiterr ) ;
//        list_ind++ ;
//     } // j
//  } // i
//  printf("\n\n\n") ;
//  for (i=0; i<dim; i++ ) {
//     double fitval = ((RooRealVar*)muVec.at(i))->getVal() ;
//     double fiterr = ((RooRealVar*)muVec.at(i))->getError() ;
//     double genval = muSetvals(i) ;
//     char parname[100] ;
//     sprintf( parname, "%s", ((RooRealVar*)muVec.at(i))->GetName() ) ;
//     printf("  %8s :  true  %8.4f,   fit  %8.4f +/- %8.4f     diff %8.4f     diff/err %7.2f\n", parname, genval, fitval, fiterr, (fitval-genval), (fitval-genval)/fiterr ) ;
//  } // i


    printf("\n\n\n") ;
}
