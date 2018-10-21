
#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"
#include "TFile.h"
#include "TROOT.h"

using namespace RooStats;
using namespace HistFactory;


/*

 A ROOT script demonstrating
 an example of writing a HistFactory
 model using c++ only.

 This example was written to match
 the example.xml analysis in
 $ROOTSYS/tutorials/histfactory/

 Written by George Lewis

 */


void example() {


  std::string InputFile1 = "./example1.root";
  std::string InputFile2= "./example2.root";
  // in case the file is not found
//  bool bfile = gSystem->AccessPathName(InputFile.c_str());
//  if (bfile) {
//     std::cout << "Input file is not found - run prepareHistFactory script " << std::endl;
//     gROOT->ProcessLine(".! prepareHistFactory .");
//     bfile = gSystem->AccessPathName(InputFile.c_str());
//     if (bfile) {
//        std::cout << "Still no " << InputFile << ", giving up.\n";
//        exit(1);
//     }
//  }

  // Create the measurement
  Measurement meas("meas", "meas");

  meas.SetOutputFilePrefix( "./results/example_UsingC_twochannel" );
  meas.SetPOI( "SigXsecOverSM" );
//  meas.AddConstantParam("alpha_syst1");
  meas.AddConstantParam("Lumi");

  meas.SetLumi( 1.0 );//Scale the histogram, if the histogram is already normalized to the luminosity, then it should be set to 1
  meas.SetLumiRelErr( 0.10 );
  meas.SetExportOnly( false );
  meas.SetBinHigh( 2 );

  // Create a channel

  Channel chan1( "channel1" );
  chan1.SetData( "data", InputFile1 );
  chan1.SetStatErrorConfig( 0.05, "Poisson" );


  // Now, create some samples


  // Create the signal sample
  Sample signal1( "zprime", "zprime", InputFile1 );
  signal1.AddOverallSys( "syst1",  0.95, 1.05 );
  signal1.AddNormFactor( "SigXsecOverSM", 1, 0, 3 );
  chan1.AddSample( signal1 );

  // Background 1
  Sample background1_1( "background1", "background1", InputFile1 );
  background1_1.ActivateStatError( "background1_statUncert", InputFile1 );
  background1_1.AddOverallSys( "syst2", 0.95, 1.05  );
  chan1.AddSample( background1_1 );


  // Background 1
  Sample background1_2( "background2", "background2", InputFile1 );
  background1_2.ActivateStatError();
  background1_2.AddOverallSys( "syst3", 0.95, 1.05  );
  chan1.AddSample( background1_2 );


  // Done with this channel
  // Add it to the measurement:
  meas.AddChannel( chan1 );

  Channel chan2( "channel2" );
  chan2.SetData( "data", InputFile2 );
  chan2.SetStatErrorConfig( 0.05, "Poisson" );

  Sample signal2( "zprime", "zprime", InputFile2 );
  signal2.AddOverallSys( "syst4",  0.95, 1.05 );
  signal2.AddNormFactor( "SigXsecOverSM", 1, 0, 3 );
  chan2.AddSample( signal2 );

  Sample background2_1( "background1", "background1", InputFile2 );
  background2_1.ActivateStatError( "background1_statUncert", InputFile2 );
  background2_1.AddOverallSys( "syst5", 0.95, 1.05  );
  chan2.AddSample( background2_1 );

  Sample background2_2( "background2", "background2", InputFile2 );
  background2_2.ActivateStatError();
  background2_2.AddOverallSys( "syst6", 0.95, 1.05  );
  chan2.AddSample( background2_2 );

  meas.AddChannel(chan2);

  // Collect the histograms from their files,
  // print some output,
  meas.CollectHistograms();
  meas.PrintTree();

  // One can print XML code to an
  // output directory:
  meas.PrintXML( "xmlFromCCode", meas.GetOutputFilePrefix() );

  // Now, do the measurement
  MakeModelAndMeasurementFast( meas );


}
