

#include "corgen.c"


   void toy_chi2_test3( double a0=0. , double b0=0.4, int nshow=25,  int ngen=10000   ) {

      TMatrixD covmat(4,4) ;

      covmat(0,0) = 1.0 ;
      covmat(1,1) = 1.0 ;
      covmat(2,2) = 1.0 ;
      covmat(3,3) = 1.0 ;

      covmat(0,3) = -0.8 ;
      covmat(3,0) = -0.8 ;

      //covmat(0,3) = 0.8 ;
      //covmat(3,0) = 0.8 ;

      printf("\n\n Covariance matrix:\n") ;
      covmat.Print() ;

      TMatrixD covInv = covmat ;
      covInv.Invert() ;
      printf("\n\n Inverse of Cov mat:\n") ;
      covInv.Print() ;

      TMatrixD invTest = covmat * covInv ;
      printf("\n\n Test of inverse:\n") ;
      invTest.Print() ;




      corgen cg( covmat ) ;

      TTree* tt = new TTree("tt","tt") ;

      double x0, x1, x2, x3 ;
      double ab, af, bb, bf, chib, chif ;


      tt->Branch( "x0", &x0, "x0/D" ) ;
      tt->Branch( "x1", &x1, "x1/D" ) ;
      tt->Branch( "x2", &x2, "x2/D" ) ;
      tt->Branch( "x3", &x3, "x3/D" ) ;

      tt->Branch( "ab", &ab, "ab/D" ) ;
      tt->Branch( "af", &af, "af/D" ) ;
      tt->Branch( "bb", &bb, "bb/D" ) ;
      tt->Branch( "bf", &bf, "bf/D" ) ;
      tt->Branch( "chib", &chib, "chib/D" ) ;
      tt->Branch( "chif", &chif, "chif/D" ) ;


      int nscan(80) ;
      double a_min(-3.) ;
      double a_max( 3.) ;
      double b_min(-4.) ;
      double b_max( 4.) ;
      double a_bw = (a_max-a_min)/(nscan-1) ;
      double b_bw = (b_max-b_min)/(nscan-1) ;

      TCanvas* can = new TCanvas( "can", "can", 50, 50, 800, 800 ) ;
      can -> Divide(2,2) ;
      can -> Draw() ;

      TH2F* h_chi2_basic_scan = new TH2F( "h_chi2_basic_scan", "chi2 basic scan", nscan, a_min-0.5*a_bw, a_max+0.5*a_bw,  nscan, b_min-0.5*b_bw, b_max+0.5*b_bw ) ;
      TH2F* h_chi2_full_scan = new TH2F( "h_chi2_full_scan", "chi2 full scan", nscan, a_min-0.5*a_bw, a_max+0.5*a_bw,  nscan, b_min-0.5*b_bw, b_max+0.5*b_bw ) ;

      TH1F* h_data = new TH1F( "h_data", "data", 4, -0.5, 3.5  ) ;
      h_data -> SetMarkerStyle(20) ;

      TLine* tl_basic = new TLine() ;
      tl_basic->SetLineColor(2) ;
      tl_basic->SetLineWidth(3) ;

      TLine* tl_full = new TLine() ;
      tl_full->SetLineColor(4) ;
      tl_full->SetLineWidth(3) ;


      for ( int i=0; i<ngen; i++ ) {

         TVectorD rnvec(4) ;

         cg.genSet( rnvec ) ;


         x0 = a0 + b0*0. + rnvec(0) ;
         x1 = a0 + b0*1. + rnvec(1) ;
         x2 = a0 + b0*2. + rnvec(2) ;
         x3 = a0 + b0*3. + rnvec(3) ;


         if ( i < nshow ) printf("  %5d :  (%7.3f, %7.3f, %7.3f)\n", i, x0, x1, x2 ) ;

         h_data -> SetBinContent(1, x0 ) ;
         h_data -> SetBinContent(2, x1 ) ;
         h_data -> SetBinContent(3, x2 ) ;
         h_data -> SetBinContent(4, x3 ) ;

         h_data -> SetBinError( 1, sqrt(covmat(0,0)) ) ;
         h_data -> SetBinError( 2, sqrt(covmat(1,1)) ) ;
         h_data -> SetBinError( 3, sqrt(covmat(2,2)) ) ;
         h_data -> SetBinError( 4, sqrt(covmat(3,3)) ) ;


         //-- scan the simple chi2
         double best_chi2_basic(1e9) ;
         double best_chi2_basic_aval(-9) ;
         double best_chi2_basic_bval(-9) ;
         for ( int ai=0; ai<nscan; ai++ ) {
            double ascan = a_min + ai*(a_max-a_min)/(nscan-1) ;
            for ( int bi=0; bi<nscan; bi++ ) {
               double bscan = b_min + bi*(b_max-b_min)/(nscan-1) ;
               double chi2_basic(0.) ;
               chi2_basic += pow( (x0-(ascan+bscan*0.)), 2 ) / covmat(0,0) ;
               chi2_basic += pow( (x1-(ascan+bscan*1.)), 2 ) / covmat(1,1) ;
               chi2_basic += pow( (x2-(ascan+bscan*2.)), 2 ) / covmat(2,2) ;
               chi2_basic += pow( (x3-(ascan+bscan*3.)), 2 ) / covmat(3,3) ;
               if ( i < nshow ) h_chi2_basic_scan -> SetBinContent( ai+1, bi+1, chi2_basic ) ;
               if ( chi2_basic < best_chi2_basic ) {
                  best_chi2_basic = chi2_basic ;
                  best_chi2_basic_aval = ascan ;
                  best_chi2_basic_bval = bscan ;
               }
            } // bi
         } // ai

         printf(" %5d : Best point:  a = %7.3f, b = %7.3f,  chi2_basic = %7.2f\n", i, best_chi2_basic_aval, best_chi2_basic_bval, best_chi2_basic ) ;


         ///   TVectorD deltaVec(3) ;
         ///   deltaVec(0) = x0 - (a0-b0) ;
         ///   deltaVec(1) = x1 - (a0) ;
         ///   deltaVec(2) = x2 - (a0+b0) ;
         ///   printf("delta vector:\n") ;
         ///   deltaVec.Print() ;
         ///   TVectorT<double> prod1 = covInv * deltaVec ;
         ///   prod1.Print() ;
         ///   double prod2 = deltaVec * prod1 ;
         ///   printf(" prod2 = %f\n", prod2 ) ;

         //-- scan the full chi2
         double best_chi2_full(1e9) ;
         double best_chi2_full_aval(-9) ;
         double best_chi2_full_bval(-9) ;
         for ( int ai=0; ai<nscan; ai++ ) {
            double ascan = a_min + ai*(a_max-a_min)/(nscan-1) ;
            for ( int bi=0; bi<nscan; bi++ ) {
               double bscan = b_min + bi*(b_max-b_min)/(nscan-1) ;
               TVectorD deltaVec(4) ;
               deltaVec(0) = x0 - (ascan+bscan*0.) ;
               deltaVec(1) = x1 - (ascan+bscan*1.) ;
               deltaVec(2) = x2 - (ascan+bscan*2.) ;
               deltaVec(3) = x3 - (ascan+bscan*3.) ;
               TVectorT<double> prod1 = covInv * deltaVec ;
               double chi2_full = deltaVec * prod1 ;
               if ( i < nshow ) h_chi2_full_scan -> SetBinContent( ai+1, bi+1, chi2_full ) ;
               if ( chi2_full < best_chi2_full ) {
                  best_chi2_full = chi2_full ;
                  best_chi2_full_aval = ascan ;
                  best_chi2_full_bval = bscan ;
               }
            } // bi
         } // ai

         printf(" %5d : Best point:  a = %7.3f, b = %7.3f,  chi2_full = %7.2f\n", i, best_chi2_full_aval, best_chi2_full_bval, best_chi2_full ) ;

         ab = best_chi2_basic_aval ;
         bb = best_chi2_basic_bval ;
         chib = best_chi2_basic ;

         af = best_chi2_full_aval ;
         bf = best_chi2_full_bval ;
         chif = best_chi2_full ;



         if ( i < nshow ) {
            can -> cd(1) ;
            h_chi2_basic_scan -> Draw("colz") ;
            tl_basic -> DrawLine( best_chi2_basic_aval, b_min, best_chi2_basic_aval, b_max ) ;
            tl_basic -> DrawLine( a_min, best_chi2_basic_bval, a_max, best_chi2_basic_bval ) ;

            can -> cd(2) ;
            h_data -> Draw("pe") ;
            tl_basic -> DrawLine( 0., best_chi2_basic_aval,  3.,  best_chi2_basic_aval + best_chi2_basic_bval*3. ) ;
            tl_full -> DrawLine( 0., best_chi2_full_aval,  3.,  best_chi2_full_aval + best_chi2_full_bval*3. ) ;

            can -> cd(3) ;
            h_chi2_full_scan -> Draw("colz") ;
            tl_full -> DrawLine( best_chi2_full_aval, b_min, best_chi2_full_aval, b_max ) ;
            tl_full -> DrawLine( a_min, best_chi2_full_bval, a_max, best_chi2_full_bval ) ;

            can -> Update() ;
            can -> Draw() ;
            gSystem -> ProcessEvents() ;
            char answ = getchar() ;
            if ( answ == 'q' ) return ;
         }

         tt->Fill() ;

      } // i

      TFile f("toy-chi2.root", "recreate") ;
      tt->Write() ;



   }







