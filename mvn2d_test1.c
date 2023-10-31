
#include "RooMvn2d.h"

using namespace RooFit;



void mvn2d_test1() {

    RooRealVar x( "x", "x", -10., 10. ) ;
    RooRealVar y( "y", "y", -10., 10. ) ;
    RooRealVar mean_x( "mean_x", "mean_x", 0.2, -3., 3. ) ;
    RooRealVar mean_y( "mean_y", "mean_y", 0.8, -3., 3. ) ;
    RooRealVar sigma2_x( "sigma2_x", "sigma2_x", 1.0*1.0, 0.1, 9. ) ;
    RooRealVar sigma2_y( "sigma2_y", "sigma2_y", 1.5*1.5, 0.1, 9. ) ;
    RooRealVar cov_xy( "cov_xy", "cov_xy", -0.6*1.0*1.5, -9., 9. ) ;

    RooMvn2d rmvn( "rmvn", "rmvn", x, y, mean_x, mean_y, sigma2_x, sigma2_y, cov_xy ) ;

    RooDataSet* data = rmvn.generate( RooArgSet(x,y), 10000 ) ;

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
