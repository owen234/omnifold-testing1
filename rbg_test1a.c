
#include "RooBinnedGauss1d.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "TRandom2.h"
#include "TCanvas.h"
#include "RooPlot.h"

using namespace RooFit;



void rbg_test1a( int rand_seed = 12345 ) {

     float true_mean_x = 0.2 ;
     float true_sig2_x = 0.9 ;

   //float true_mean_x = 0.0 ;
   //float true_sig2_x = 1.0 ;

   float true_sig_x = sqrt( true_sig2_x ) ;

   int ngen = 1000 ;


   int hist_nbins = 90 ;
   //int hist_nbins = 190 ;
   //int hist_nbins = 900 ;

   float hist_xlow = -9. ;
   float hist_xhigh = 9. ;

   float hist_binwidth = (hist_xhigh - hist_xlow) / hist_nbins ;

   float hist_bin_centers[190] ;

   for ( int bi=0; bi<hist_nbins; bi++ ) {
      hist_bin_centers[bi] = hist_xlow + (bi+0.5)*hist_binwidth ;
      printf( "  %3d : bin center  %10.3f\n", bi, hist_bin_centers[bi] ) ;
   } // bi


   RooRealVar x( "x", "x", 0., -10., 10. ) ;
   RooRealVar mu( "mu", "mu", 0., -10., 10. ) ;
   RooRealVar sig2( "sig2", "sig2", 0.7, 0.1, 2.0 ) ;

   RooArgList bin_edges ;
   for ( int i=0; i<=hist_nbins; i++ ) {
      RooConstVar* be = new RooConstVar( Form("be%03d", i), Form("be%03d", i), hist_xlow + i*hist_binwidth ) ;
      bin_edges.add( *be ) ;
   }


   RooBinnedGauss1d rbg( "rbg", "rbg", x, mu, sig2, hist_nbins, bin_edges ) ;


   //--------------

   TRandom2 tran( rand_seed ) ;

   printf("\n\n Generating data.\n\n") ;
   fflush(stdout) ;

   RooDataSet* data = new RooDataSet( "ds", "ds", x ) ;
   for (int ei=0; ei<ngen; ei++ ) {
      float grn = tran.Gaus( true_mean_x, true_sig_x ) ;
      int bin_index = int( grn - hist_xlow ) / hist_binwidth ;
      if ( bin_index < 0 ) bin_index = 0 ;
      if ( bin_index >= hist_nbins ) bin_index = hist_nbins-1 ;
      float dgrn = hist_bin_centers[ bin_index ] ;
      x.setVal( dgrn ) ;
      data -> add( x ) ;
   } // ei


    data->Print("V") ;
    data->get(0) -> Print("V") ;
    data->get(1) -> Print("V") ;
    data->get(2) -> Print("V") ;
    data->get(3) -> Print("V") ;
    data->get(4) -> Print("V") ;

    printf("\n\n Fitting test data.\n\n") ;
    fflush(stdout) ;

    //RooFitResult* rfr = rbg.fitTo( *data, NumCPU(20), Save(true) ) ;
    RooFitResult* rfr = rbg.fitTo( *data, Save(true) ) ;

  return ;

    TCanvas* can = new TCanvas("can","can", 50, 50, 1200, 800 ) ;
    can -> Divide(2,1) ;


    can -> cd(1) ;
    RooPlot* xframe = x.frame(Title("x frame")) ;
    data -> plotOn(xframe) ;
    //rbg.plotOn(xframe, NumCPU(20)) ;
    xframe -> Draw() ;

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
