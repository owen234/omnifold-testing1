
#include "RooMvn2d.h"

using namespace RooFit;



void mvn2d_fit_bootstrap( int bsi=0, const char* infile = "bootstrap-2d-export1c.root" ) {

    TChain ch("boot_ttree") ;
    ch.Add( infile ) ;


    float x0 ;
    float x1 ;
    float tt_weight ;
    int   dset_index ;
    ch.SetBranchAddress("x0",&x0) ;
    ch.SetBranchAddress("x1",&x1) ;
    ch.SetBranchAddress("weight",&tt_weight) ;
    ch.SetBranchAddress("dset_index",&dset_index) ;

    //ch.GetEntry(0) ;
    //printf(" event 0 :  x0 = %8.3f, x1 = %8.3f, weight = %8.3f, dset_index = %d\n", x0, x1, weight, dset_index ) ;




    RooRealVar x( "x", "x", -10., 10. ) ;
    RooRealVar y( "y", "y", -10., 10. ) ;
    RooRealVar mean_x( "mean_x", "mean_x", 0.2, -3., 3. ) ;
    RooRealVar mean_y( "mean_y", "mean_y", 0.8, -3., 3. ) ;
    RooRealVar sigma2_x( "sigma2_x", "sigma2_x", 1.0*1.0, 0.1, 9. ) ;
    RooRealVar sigma2_y( "sigma2_y", "sigma2_y", 1.5*1.5, 0.1, 9. ) ;
    RooRealVar cov_xy( "cov_xy", "cov_xy", -0.6*1.0*1.5, -9., 9. ) ;
    RooRealVar weight( "weight", "weight", 0., 100. ) ;

    RooMvn2d rmvn( "rmvn", "rmvn", x, y, mean_x, mean_y, sigma2_x, sigma2_y, cov_xy ) ;

    //////////RooDataSet* data = rmvn.generate( RooArgSet(x,y), 10000 ) ;
    /////////RooDataSet* data = new RooDataSet( "ds", "ds", RooArgSet(x,y), Import(ch), Cut(Form("dset_index==%d", bsi)) ) ; //// NFG

    RooDataSet* data = new RooDataSet( "ds", "ds", RooArgSet(x,y,weight), WeightVar("weight") ) ;
    int nadd(0) ;
    for ( int i=0; i<ch.GetEntries(); i++ ) {
       ch.GetEntry(i) ;
       if ( dset_index != bsi ) continue ;
       x.setVal(x0) ;
       y.setVal(x1) ;
       weight.setVal(tt_weight) ;
       data -> add( RooArgSet( x, y, weight ), weight.getVal() ) ;
       nadd++ ;
    }
    printf("\n\n Added %d events.\n\n", nadd ) ;


    data->Print("V") ;
    data->get(0) -> Print("V") ;
    data->get(1) -> Print("V") ;
    data->get(2) -> Print("V") ;
    data->get(3) -> Print("V") ;
    data->get(4) -> Print("V") ;

    RooFitResult* rfr = rmvn.fitTo( *data, Save(true) ) ;

    mean_x.Print() ;
    mean_y.Print() ;
    sigma2_x.Print() ;
    sigma2_y.Print() ;
    cov_xy.Print() ;


    TCanvas* can = new TCanvas("can","can", 50, 50, 1200, 800 ) ;
    can -> Divide(2,1) ;


    can -> cd(1) ;
    RooPlot* xframe = x.frame(Title("data")) ;
    data -> plotOn(xframe) ;
    rmvn.plotOn(xframe) ;
    xframe -> Draw() ;

    can -> cd(2) ;
    RooPlot* yframe = y.frame(Title("data")) ;
    data -> plotOn(yframe) ;
    rmvn.plotOn(yframe) ;
    yframe -> Draw() ;



}
