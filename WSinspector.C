/*
Author: Valerio Dao

Based on the many features of FitCrossChecksForLimits from:
Email:  romain.madar@cern.ch, gabriel.facini@cern.ch


Description : 
- it dumps important information out of a workspace (
basically it allows to check a workspace in standalone mode without relying on underlying histograms)

*/


// C++
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <map>

// Root
#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TList.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGaxis.h"

// RooFit
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooAbsData.h"
#include "RooHist.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooFitResult.h"
#include "RooAbsData.h"
#include "RooRealSumPdf.h"
#include "Roo1DTable.h"
#include "RooConstVar.h"
#include "RooProduct.h"
#include "RooRandom.h"
#include "TStopwatch.h"
#include "RooNLLVar.h"
#include "RooMsgService.h"
#include "RooMinimizer.h"
///#include "HistFactory/FlexibleInterpVar.h"

#include <vector>

// RooStat
#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileInspector.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SamplingDistribution.h"
#include "RooStats/SamplingDistPlot.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/RooStatsUtils.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;

struct NPContainer{ 
  TString NPname; 
  double  NPvalue; 
  double  NPerrorHi;
  double  NPerrorLo;
  TString WhichFit;
};

static bool comp_second_abs_decend( const pair< RooRealVar*, float >& i, const pair< RooRealVar*, float >& j ) {
  return fabs(i.second) > fabs(j.second);
}

namespace LimitCrossCheck{  
  
  bool verbose=true;
  // Global variables;
  // User configuration one
  bool isBinned=true;                       // speeds the computation of binned likelihoods
  bool isGG;
  bool isComb;
   
  bool drawPlots(true);                    // create eps & png files and creat a webpage
  bool plotRelative(false);                 // plot % shift of systematic
  int isBlind(0);                           // 0: Use observed Data 1: use Asimov data 2: use toydata
  double mu_asimov(1.0);                    // mu value used to generate Asimov dataset (not used if isBlind==0)
  TString xAxisLabel("Final Distribution"); // set what the x-axis of the distribution is


  ////////////////////////////////////////////////////////////////////////////////////
  TH2F* dataStat;
  TH2F* mcStat;
  TFile*    outputfile; 

  // tmp infos
  int nCategories=0;
  int nNP=0;
  vector<int> nBins;
  
  // not switches
  RooWorkspace *w         ;
  ModelConfig  *mc        ;
  RooAbsData   *data      ;
  RooSimultaneous* pdf    ;
  RooCategory* channelCat ;
  RooRealVar * firstPOI   ;


  double        LumiRelError;
  TDirectory   *MainDirSyst;
  TDirectory   *MainDirFitGlobal;

  map <string,double>         MapNuisanceParamNom;
  map <TString,RooFitResult*> AllFitResults_map;
  map <TString,int>           AllFitStatus_map;
  vector<NPContainer>         AllNPafterEachFit_vec;
  TString OutputDir;
  
  //better debugging
  vector<string> sysNames;
  vector<string> sampleNames;
  vector<string> regionNames;
  map< string, vector< vector<float> > > sysEffect_up;
  map< string, vector< vector<float> > > sysEffect_do;  

  //Global functions
  void     PrintModelObservables();
  void     PrintNuisanceParameters();
  void     PrintAllParametersAndValues(RooArgSet para);
  void     PrintNumberOfEvents(RooAbsPdf *Tpdf);
  void     PrintSubChannels();

  void FillSysInfo();
  RooRealSumPdf* getModelPDF(RooAbsPdf* thePDF, TString catName);

  void PrintSysPerSample(string samName);

  string   decodeRegionName(string regName);
 
  bool     IsSimultaneousPdfOK();
  bool     IsChannelNameOK();

  void     SetAllNuisanceParaToSigma(double Nsigma);
  void     SetAllStatErrorToSigma(double Nsigma);
  void     SetNuisanceParaToSigma(RooRealVar *var, double Nsigma);
  void     GetNominalValueNuisancePara();
  void     SetNominalValueNuisancePara();
  void     SetPOI(double mu);
  void     SetStyle();
 
  bool     IsAnormFactor(RooRealVar *var);
   
  void     Initialize(const char* infile , const char* outputdir, const char* workspaceName, const char* modelConfigName, const char* ObsDataName);
  
  //// new methods
  void     FixUnwantedNF();
  void     PrintYields(string type);
  void     PrintSystematics(string type);


  RooRealSumPdf* getModelPDF(RooAbsPdf* thePDF, TString catName) {
    RooArgSet* tmpSet= thePDF->getComponents(); 
    TIterator* iter  = tmpSet->createIterator() ;
    RooAbsArg* MyObs = NULL;
    while( (MyObs = (RooAbsArg*) iter->Next()) ) {
      TString obsNam=MyObs->GetName();
      if ( obsNam.Contains( catName+"_model") ) {
	RooRealSumPdf* pdfmodel= (RooRealSumPdf*)MyObs;
	return pdfmodel;
      }
    }
    return NULL;
  }
  




  //==============================================================================================  //============================= Main function =====================
  //==============================================================================================
  void PlotFitCrossChecks(const char* infile          = "WorkspaceForTest1.root",
			  const char* outputdir       = "./results/",
			  const char* workspaceName   = "combined",
			  const char* modelConfigName = "ModelConfig",
			  const char* ObsDataName     = "obsData"){
    
    TString tmpFolder=outputdir;
    if (tmpFolder.Contains("Asimov") ) isBlind=1;
    
    // are these still needed????
    isGG=false;
    TString tmpName=infile;
    if (tmpName.Contains("yy") ||tmpName.Contains("GG") || tmpName.Contains("FULL") ) isGG=true;
    isComb=false;
    if (tmpName.Contains("combined") ) isComb=true;

    Initialize(infile, outputdir, workspaceName, modelConfigName, ObsDataName);

    //PlotHistosBeforeFit(0,1.0);

    return; 
  }
  

  string  decodeRegionName(string regName) {
    string reg="";
    if ( regName.find("Fat1")!=string::npos ) reg+="Mer_";
    else                                      reg+="Res_";
    if ( regName.find("BMin500")!=string::npos )      reg+="HiPtV_";
    else if ( regName.find("BMax500")!=string::npos ) reg+="LoPtV_";
    else                                              reg+="InPtV_";
    if ( regName.find("T3_")!=string::npos )       reg+="3T_";
    else if ( regName.find("T2_")!=string::npos )  reg+="2T_";
    else if ( regName.find("T1_")!=string::npos )  reg+="1T_";
    else                                           reg+="12T_";

    if ( regName.find("SR")!=string::npos )       reg+="SR";
    else if (regName.find("topem")!=string::npos) reg+="CRtop";
    else if (regName.find("mBB")!=string::npos)   reg+="CRmbb";

    return reg;
    TString newReg=regName;
    newReg=newReg.ReplaceAll("Region_","");
    newReg=newReg.ReplaceAll("_Y2015","");
    newReg=newReg.ReplaceAll("_distmva_DSR","");
    return newReg.Data();
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////
  void FixUnwantedNF() {
    if (isGG || isComb) {
      string listToFix[]={"mu_7TeV","mu_8TeV","mu_BR_gg","mu_VBFVH","mu_VBF_7TeV","mu_VBF_8TeV","mu_VBF","mu_VH","mu_WH_7TeV","mu_WH_8TeV","mu_WH","mu_ZH_7TeV","mu_ZH_8TeV","mu_ZH","mu_bbH_7TeV","mu_bbH_8TeV","mu_bbH","mu_ggH_7TeV","mu_ggH_8TeV","mu_ggH","mu_ggHttH","mu_tH_7TeV","mu_tH_8TeV","mu_tH","mu_ttH_7TeV","mu_ttH_8TeV","mu_ttH_C11toC12","mu_tHjb","mu_tHjb_8TeV","mu_tH","mu_tH_8TeV","mu_WtH","mu_tHjb_7TeV","mu_tH_7TeV","mu_WtH_7TeV"};
      int nFix=sizeof(listToFix)/sizeof(listToFix[0]);
      w->var("m_H")->setVal(125.4);
      w->var("m_H")->setConstant(true);
      w->var("mu_ttH")->setVal(1);
      w->var("mu_ttH")->setConstant(false);
      if (isGG) {
	w->var("mu")->setVal(1);
	w->var("mu")->setConstant(true);
      } else {
	w->var("mu_Incl")->setVal(1);
	w->var("mu_Incl")->setConstant(true);
      }  
      for (int j=0; j<nFix; j++) w->var(listToFix[j].c_str())->setConstant(true);
    } else {
      w->var("massCalib")->setVal(0.98);
      //w->var("massCalib")->setVal(1.0);
      w->var("massCalib")->setConstant(true);
      w->var("massCalib")->Print();
    }
  }
 
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void GetNominalValueNuisancePara(){
    TIterator *it = mc->GetNuisanceParameters()->createIterator();
    RooRealVar *var = NULL;
    if (MapNuisanceParamNom.size() > 0) MapNuisanceParamNom.clear();
    while ((var = (RooRealVar*)it->Next()) != NULL){
      const double val = var->getVal();
      MapNuisanceParamNom[(string)var->GetName()] = val;
    }
    return;
  }
  

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void SetNominalValueNuisancePara(){
    TIterator *it = mc->GetNuisanceParameters()->createIterator();
    RooRealVar *var = NULL;
    while ((var = (RooRealVar*)it->Next()) != NULL){
      string varname = (string)var->GetName();
      const double val =  MapNuisanceParamNom[varname];
      var->setVal(val);
    }
    return;
  }
  

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void SetAllStatErrorToSigma(double Nsigma){

    TIterator* it = mc->GetNuisanceParameters()->createIterator();
    RooRealVar* var = NULL;
    while( (var = (RooRealVar*) it->Next()) ){
      string varname = (string) var->GetName();
      if ( varname.find("gamma_stat")!=string::npos ){
	RooAbsReal* nom_gamma = (RooConstVar*) w->obj( ("nom_" + varname).c_str() );
	double nom_gamma_val = nom_gamma->getVal();
	//cout << " name: " << nom_gamma->GetName() << "  VALUE: " << nom_gamma_val << endl;
	double sigma = 1/TMath::Sqrt( nom_gamma_val );
	var->setVal(1 + Nsigma*sigma);
        }
    }
    
    return;
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void SetAllNuisanceParaToSigma(double Nsigma){
    
    TIterator* it = mc->GetNuisanceParameters()->createIterator();
    RooRealVar* var = NULL;
    while( (var = (RooRealVar*) it->Next()) ){
      string varname = (string) var->GetName();
      if ( varname.find("gamma_stat")!=string::npos ) continue;
      if(strcmp(var->GetName(),"Lumi")==0){
	var->setVal(w->var("nominalLumi")->getVal()*(1+Nsigma*LumiRelError));
      }
      else if (IsAnormFactor(var) && Nsigma==0) var->setVal(1.); // nominal for normfactor is 1.0
      else if (IsAnormFactor(var)) continue; // don't touch it if this is not nominal
      else  var->setVal(Nsigma);
    }
    
    return;
  }
  

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void SetNuisanceParaToSigma(RooRealVar *var, double Nsigma){
    
    string varname = (string) var->GetName();
    if ( varname.find("gamma_stat")!=string::npos ) return;
    
    if(strcmp(var->GetName(),"Lumi")==0){
      var->setVal(w->var("nominalLumi")->getVal()*(1+Nsigma*LumiRelError));
    } 
    else if (IsAnormFactor(var)) return;
    else var->setVal(Nsigma);

    return;
  }

  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  bool IsAnormFactor(RooRealVar *var){
    // potentially faulty method!!!!!!!!!
    bool result=false;
    string varname = (string) var->GetName();
    const double val =  MapNuisanceParamNom[varname];
    if (!((TString)varname).Contains("_stat_") && val==1.0) result=true;  
    return result;
  }

  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void SetPOI(double mu){
    firstPOI->setVal(mu);
    return  ;
  }


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  bool IsSimultaneousPdfOK(){
    bool IsSimultaneousPDF = strcmp(pdf->ClassName(),"RooSimultaneous")==0;
    if (!IsSimultaneousPDF){
      cout << " ERROR : no Simultaneous PDF was found, will stop here." << endl;
      cout << " You need to investigate your input histogramms." << endl;
      return false;
    }
    return true;
  }

  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  bool IsChannelNameOK(){
    if( !IsSimultaneousPdfOK() ) return false;
    
    TIterator* iter = channelCat->typeIterator() ;
    RooCatType* tt = NULL;    
    while((tt=(RooCatType*) iter->Next()) ){
      string channelName =  tt->GetName();
      if (channelName.find("/")!=string::npos){
	cout << endl;
	cout << "One of the channel name contain a caracter \"/\" : " << endl;
	cout << "  - "  << channelName << endl;
	cout << "This is mis-intrepreted by roofit in the reading of the workspace. " << endl;
	cout << "Please change the channel name in the xml file to run this code." << endl;
	cout << endl;
	return false;
      }
    }
   return true;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void PrintSysPerSample(string samName) {
    if ( sysEffect_up.find(samName)== sysEffect_up.end() ) {
      cout << " could not find sample: " << samName << " in the default sample map" << endl;
      return;
    }
    Float_t tot_up[regionNames.size()], tot_dn[regionNames.size()];
    vector< vector<float> > tmpV_up=sysEffect_up[samName];
    vector< vector<float> > tmpV_do=sysEffect_do[samName];
    cout << "=========================================================================" << endl;
    cout << "  Printing sys for sample:     ' " << samName << " ' " << endl;
    cout << "=========================================================================" << endl;
    string empty=" ";
    cout << Form(" %-40s |"," sys , region ");
    for (unsigned int iS=0; iS<regionNames.size(); iS++) {
      cout << Form(" %-18s |", (regionNames.at(iS)).c_str() );
    }   
    cout << endl;
    cout << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
   
    for (unsigned int iSys=0; iSys<sysNames.size(); iSys++) {
      TString tmpName=sysNames.at(iSys);
      tmpName=tmpName.ReplaceAll("alpha_","");
      tmpName=tmpName.ReplaceAll("Intercalibration","Inter");
      string newString=tmpName.Data();
      cout << Form(" %-40s |", newString.c_str() );
      //cout << Form(" %-40s |", sysNames.at(iSys).c_str() );
      for (unsigned int iS=0; iS<regionNames.size(); iS++) {
	if ( fabs(tmpV_up[iS][iSys])<1e-5 && fabs(tmpV_up[iS][iSys])<1e-5 )      cout << setw(21) << " -- / --  |";
	else if ( fabs(tmpV_up[iS][iSys])>1e-5 && fabs(tmpV_up[iS][iSys])<1e-5 ) cout << setw(21) << Form("%2.1f / -- |", tmpV_up[iS][iSys] );
	else if ( fabs(tmpV_up[iS][iSys])<1e-5 && fabs(tmpV_up[iS][iSys])>1e-5 ) cout << setw(21) << Form("-- / %2.1f |", tmpV_do[iS][iSys] );
	else                                                                     cout << setw(21) << Form("%2.1f / %2.1f |",  tmpV_up[iS][iSys], tmpV_do[iS][iSys] );
        //if ( std::find(newString.begin(), newString.end(), "mu") != newString.end() ) {
	//  tot_up[iS] += std::pow(tmpV_up[iS][iSys],2);
	//  tot_dn[iS] += std::pow(tmpV_do[iS][iSys],2);
	//}
      }
      cout << endl;
    }
    
    //for (unsigned int iS=0; iS<regionNames.size(); iS++) {
    //  cout << "Tot_up [" << decodeRegionName( (regionNames.at(iS)) ).c_str() << "] = " << sqrt(tot_up[iS]) << std::endl;
    //  cout << "Tot_dn [" << decodeRegionName( (regionNames.at(iS)) ).c_str() << "] = " << sqrt(tot_dn[iS]) << std::endl;
    //}

    cout << endl;
  }


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void FillSysInfo() {
    sysEffect_up.clear();
    sysEffect_do.clear();
    for (unsigned int sam=0; sam<sampleNames.size(); sam++) {
      vector< vector<float> > tmpVV;
      for (unsigned int reg=0; reg<regionNames.size(); reg++) {
	vector<float> tmpV;
	tmpV.resize( sysNames.size() );
	tmpVV.push_back(tmpV);
      } //regions
      sysEffect_up[sampleNames.at(sam)]=tmpVV;
      sysEffect_do[sampleNames.at(sam)]=tmpVV;
    } //samples
    
    ////////////////////////////////////////////////////////////////////////////////////////
    TIterator* iter = channelCat->typeIterator() ;
    RooCatType* tt = NULL;   
    int count=-1;
    while((tt=(RooCatType*) iter->Next()) ) {
      count++;
      cout << " category: " << tt->GetName() << endl;
      RooAbsPdf  *pdfReg  = pdf->getPdf( tt->GetName() );
      RooArgSet  *obstmp  = pdfReg->getObservables( *mc->GetObservables() ) ;
      RooRealVar *obs     = ((RooRealVar*) obstmp->first());
      cout << "   OBS2: " << tt->GetName() << endl;
      SetPOI(0.0);
      float nominal= pdfReg->expectedEvents(*obs);
      SetPOI(1.0);
      float nominalSig= pdfReg->expectedEvents(*obs)-nominal;
      SetPOI(0.0);
      //cout << "   nominal sig: " << tt->GetName() << endl;
      for (unsigned int iSys=0; iSys<sysNames.size(); iSys++) {
	//cout << " ... ..... " << sysNames[iSys] << endl;
	float iniV=w->var( sysNames.at(iSys).c_str() )->getVal();
	if ( iniV==0 ) w->var( sysNames.at(iSys).c_str() )->setVal(+1);
	else           w->var( sysNames.at(iSys).c_str() )->setVal(+2);
	float sys_up= pdfReg->expectedEvents(*obs);
	SetPOI(1.0);
	float sysALL_up= pdfReg->expectedEvents(*obs)-sys_up;
	//cout << " ... up to this point: " << endl;
	(sysEffect_up["background"])[count][iSys]= ((sys_up/nominal)-1)*100;
	(sysEffect_up["signal"])[count][iSys]= ((sysALL_up/nominalSig)-1)*100;
	SetPOI(0.0);
	//cout << " ... up to here: " << endl;

	if ( iniV==0 ) w->var( sysNames.at(iSys).c_str() )->setVal(-1);
	else           w->var( sysNames.at(iSys).c_str() )->setVal(0);
	float sys_do= pdfReg->expectedEvents(*obs);
	SetPOI(1.0);
	float sysALL_do= pdfReg->expectedEvents(*obs)-sys_do;
	(sysEffect_do["background"])[count][iSys]= ((sys_do/nominal)-1)*100;
	(sysEffect_do["signal"])[count][iSys]= ((sysALL_do/nominalSig)-1)*100;
	SetPOI(0.0);
	w->var( sysNames.at(iSys).c_str() )->setVal(iniV);
      }
      //cout << "   DONE sys: " << tt->GetName() << endl;
      SetPOI(1.0);
      //now each samples
      TString catName=tt->GetName();
      RooRealSumPdf* pdfmodel=getModelPDF(pdfReg,tt->GetName());
      RooArgList funcList =  pdfmodel->funcList();
      RooLinkedListIter funcIter = funcList.iterator() ;
      RooProduct* comp = 0;
      while( (comp = (RooProduct*) funcIter.Next())) { 
	TString compName=comp->GetName();
	TString newName( compName(0,compName.Index("_"+catName) ) );
	newName=newName.ReplaceAll( "L_x_", "" );
	float nomin=(comp->createIntegral(*obstmp))->getVal();
	cout << "  testVale:  for " << newName << " : " << nomin << endl; 
	for (unsigned int iSys=0; iSys<sysNames.size(); iSys++) {
	  float iniV=w->var( sysNames.at(iSys).c_str() )->getVal();
	  if ( iniV==0 ) w->var( sysNames.at(iSys).c_str() )->setVal(+1);
	  else           w->var( sysNames.at(iSys).c_str() )->setVal(+2);
	  float sys_up= (comp->createIntegral(*obstmp))->getVal();
	  (sysEffect_up[newName.Data()])[count][iSys]= ((sys_up/nomin)-1)*100;
	  if ( iniV==0 ) w->var( sysNames.at(iSys).c_str() )->setVal(-1);
	  else           w->var( sysNames.at(iSys).c_str() )->setVal(0);
	  float sys_do= (comp->createIntegral(*obstmp))->getVal();
	  (sysEffect_do[newName.Data()])[count][iSys]= ((sys_do/nomin)-1)*100;
	  w->var( sysNames.at(iSys).c_str() )->setVal(iniV);
	}
      }
    }
  }


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void PrintModelObservables(){
    regionNames.clear();
    sampleNames.clear();
    
    RooArgSet* AllObservables = (RooArgSet*) mc->GetObservables();
    TIterator* iter = AllObservables->createIterator() ;
    RooAbsArg* MyObs = NULL;    
    int TotBin=0;
    int maxBin=-1;
    cout << endl;
    cout << "------------------------------------------------------------------------"  << endl;
    cout << "  List of model Observables : "  << endl;
    cout << "------------------------------------------------------------------------"  << endl;
    nCategories=0;
    while( (MyObs = (RooAbsArg*) iter->Next()) ) {
      string obsNam=MyObs->GetName();
      if (obsNam.find("obs")==string::npos) continue;
      TString tmpNam=obsNam;
      TString stripS=tmpNam.ReplaceAll("obs_x_","");
      string tmpString=stripS.Data();
      //cout << "FUCKING FILLING WITH: " << tmpString  << endl;
      //regionNames.push_back( tmpString );
      cout << setw(80) << tmpNam.ReplaceAll("obs_x_","") << " HAS: " << setw(6) << (((RooRealVar*)MyObs)->getBinning()).numBins() << " bins" << endl;
      int tmpBins=(((RooRealVar*)MyObs)->getBinning()).numBins();
      TotBin+=tmpBins;
      nCategories++;
      if (tmpBins>=maxBin) maxBin=tmpBins;
      nBins.push_back(tmpBins);
      
    }
    cout << "Total Number of categories: " << nCategories << "  total number of bins is: " << TotBin << endl << endl;

    cout << "------------------------------------------------------------------------"  << endl;
    cout << " Samples in each region : "  << endl;
    cout << "------------------------------------------------------------------------"  << endl;
    regionNames.clear();
    iter = channelCat->typeIterator() ;
    RooCatType* tt = NULL;   
    int cat=0;
    set<string> tmpSamples;
    while((tt=(RooCatType*) iter->Next()) ) {    
      TString catName=tt->GetName();      
      RooAbsPdf  *pdftmp  = pdf->getPdf(tt->GetName()) ;
      RooRealSumPdf* pdfmodel=getModelPDF(pdftmp,catName);
      cout << "REGION: " << catName << " has components: " << endl;
      regionNames.push_back( catName.Data() );
      string listName="";
      RooArgList funcList =  pdfmodel->funcList();
      RooLinkedListIter funcIter = funcList.iterator() ;
      RooProduct* comp = 0;
      while( (comp = (RooProduct*) funcIter.Next())) { 
	TString compName=comp->GetName();
	TString newName( compName(0,compName.Index("_"+catName)) );
	newName=newName.ReplaceAll( "L_x_", "" );
	cout << " " << newName << " , ";
	tmpSamples.insert( newName.Data() );
      }
      cout << endl;
    }
    set<string>::iterator Itr=tmpSamples.begin();
    for ( ; Itr!=tmpSamples.end(); ++Itr) {
      sampleNames.push_back(*Itr);
    }
    sampleNames.push_back("signal");
    sampleNames.push_back("background");

    dataStat=new TH2F("dataStat","dataStat",maxBin,-0.5,maxBin-0.5,nCategories,-0.5,nCategories-0.5);
    mcStat=new TH2F("mcStat","mcStat",maxBin,-0.5,maxBin-0.5,nCategories,-0.5,nCategories-0.5);
    // put axis labels ...		      
    cout << endl;
    return;
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void PrintNuisanceParameters(){
    RooRealVar* arg;
    sysNames.clear();
    //w->var("normNF_ttll")->Print("v");
    //w->var("normNF_Lumi")->Print("v");

    RooArgSet nuis = *mc->GetNuisanceParameters();
    TIterator* itr = nuis.createIterator();
    nNP=0;
    int nNF=0;
    int totGamma=0;
    
    cout << endl;
    cout << "-------------------------------------------------------------------"  << endl;
    cout << "List of nuisance parameters : "  << endl;
    cout << "-------------------------------------------------------------------"  << endl;
    while ((arg=(RooRealVar*)itr->Next())) {
      if (!arg) continue;
      string name=arg->GetName();
      TString obsNam=arg->GetName();
    
      if (name.find("gamma")!=string::npos) {
	if (verbose) cout << arg->GetName()  << " : " << arg->getVal() << "+/-" << arg->getError() << "     IsNormFactor="<< IsAnormFactor(arg) << endl;
	totGamma++;
	continue;
      }
      cout << setw(45) << arg->GetName()  << " : " << arg->getVal() << "   Err:" << arg->getError() << "  Constant: " << arg->isConstant() << endl;   
      //IsNormFactor="<< IsAnormFactor(arg) << endl;
      //arg->Print("v");
      nNP++;
      //if ( obsNam.Contains("JER") )
      sysNames.push_back( obsNam.Data() );
    }
    cout << "------------------------------------" << endl;
    cout << "Total Number of Gammas: " << totGamma << endl;
    cout << "Total Number of NP: " << nNP << endl;
    
    
    cout << endl << "-------------------------------------------------------------------"  << endl;
    cout << "  List of constant parameters : "  << endl;
    cout << "-------------------------------------------------------------------"  << endl;


    //other:
    RooArgSet other = w->allVars();;
    TIterator* itr2 = other.createIterator();
    vector<string> toBeFixed;
    while ((arg=(RooRealVar*)itr2->Next())) {
      if (!arg) continue;
      string name=arg->GetName();
      //if ( arg->isConstant() ) {
	if ( name.find("nom")!=string::npos )   continue;
	if ( name.find("Width")!=string::npos ) continue;
	if ( name.find("gamma")!=string::npos ) continue;
	if ( name.find("obs")!=string::npos )   continue;
	if ( name.find("alpha")!=string::npos ) continue;
	cout << setw(45) << arg->GetName()  << " : " << arg->getVal() << "   Err:" << arg->getError() << "  Constant: " << arg->isConstant() << endl;   
	//sysNames.push_back( name );
	//w->var(arg->GetName())->Print();
	//}
	if ( !arg->isConstant() )  toBeFixed.push_back(name);
    }
        
    sysNames.push_back( mc->GetParametersOfInterest()->first()->GetName() );
    cout << endl;
    for (unsigned int i=0; i<toBeFixed.size(); i++) {
      cout << " \"" << toBeFixed[i] << "\", " ;
    } 
    cout << endl;
  

    cout << endl << "-------------------------------------------------------------------"  << endl;
    cout << "  List of parameters of interests : "  << endl;
    cout << "-------------------------------------------------------------------"  << endl;
    RooArgSet nuis2 = *mc->GetParametersOfInterest();
    itr = nuis2.createIterator();
    while ((arg=(RooRealVar*)itr->Next())) {
      if (!arg) continue;
      cout << setw(15) << arg->GetName()  << " : " << arg->getVal() << "   Err:" << arg->getError() << "  Constant: " << arg->isConstant() << endl;
      ///IsNormFactor="<< IsAnormFactor(arg) << endl;
      arg->Print();
  
    }
    cout << endl;
    return;
  }
  
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void PrintAllParametersAndValues(RooArgSet para){
    TIterator* itr = para.createIterator();
    RooRealVar* arg;
    cout << endl;
    cout << "List of parameters : "  << endl;
    cout << "----------------------------"  << endl;    
    while ((arg=(RooRealVar*)itr->Next())) {
      if (!arg) continue;
      cout << arg->GetName() << " = " << arg->getVal() << endl;
    }
    return;
  }


  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void PrintSubChannels(){
    
    TIterator* iter = channelCat->typeIterator() ;
    RooCatType* tt = NULL;    
    
    while((tt=(RooCatType*) iter->Next()) ){

      RooAbsPdf  *pdftmp  = pdf->getPdf( tt->GetName() );
      cout << endl;
      cout << "Details on channel " << tt->GetName() << " : "  << endl;
      cout << "----------------------------------------------------------" << endl;   
      PrintNumberOfEvents(pdftmp);
    }
    
    return;
  }


  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void PrintNumberOfEvents(RooAbsPdf *Tpdf){
    
    double val_sym=1;
      cout 
        << Form(" %3s |","")
        << Form(" %-32s |","Nuisance Parameter") 
        << Form(" %18s |","Signal events") 
        << Form(" %18s |","% Change (+1sig)") 
        << Form(" %18s |","% Change (-1sig)") 
        << Form(" %18s |","Background events") 
        << Form(" %18s |","% Change (+1sig)") 
        << Form(" %18s |","% Change (-1sig)") 
        << endl;

      int inuis=-1;
      RooArgSet  *obstmp  = Tpdf->getObservables( *mc->GetObservables() ) ;
      RooRealVar *myobs   = ((RooRealVar*) obstmp->first());

      RooArgSet nuis = *mc->GetNuisanceParameters();
      TIterator* itr = nuis.createIterator();
      RooRealVar* arg;
      while ((arg=(RooRealVar*)itr->Next())) {
        if (!arg) continue;
        //
        ++inuis;
        //

        double val_hi = val_sym;
        double val_lo = -val_sym;
        double val_nom = arg->getVal();
        if (string(arg->GetName()) == "Lumi"){
          val_nom = w->var("nominalLumi")->getVal();
          val_hi  = w->var("nominalLumi")->getVal() * (1+LumiRelError);
          val_lo  = w->var("nominalLumi")->getVal() * (1-LumiRelError);
        }
	
	if (string(arg->GetName()).find("ttH")==string::npos) continue;

        //
        arg->setVal(val_hi);
        firstPOI->setVal(0);
        double b_hi = Tpdf->expectedEvents(*myobs);
        firstPOI->setVal(1);
        double s_hi = Tpdf->expectedEvents(*myobs)-b_hi;
        //
        arg->setVal(val_lo);
        firstPOI->setVal(0);
        double b_lo = Tpdf->expectedEvents(*myobs);
        firstPOI->setVal(1);
        double s_lo = Tpdf->expectedEvents(*myobs)-b_lo;
        //
        arg->setVal(val_nom);
        firstPOI->setVal(0);
        double b_nom = Tpdf->expectedEvents(*myobs);
        firstPOI->setVal(1);
        double s_nom = Tpdf->expectedEvents(*myobs)-b_nom;
        //
        double x_nom = s_nom ;
        double x_hi  = 0; if (s_nom) x_hi = (s_hi-s_nom)/s_nom; 
        double x_lo  = 0; if (s_nom) x_lo = (s_lo-s_nom)/s_nom; 
        double y_nom = b_nom ;
        double y_hi  = 0; if (b_nom) y_hi = (b_hi-b_nom)/b_nom; 
        double y_lo  = 0; if (b_nom) y_lo = (b_lo-b_nom)/b_nom; 

        cout 
          << Form(" %3d |",inuis)
          << Form(" %-32s |",arg->GetName()) 
          << Form(" %18.2f |",x_nom) 
          << Form(" %18.2f |",100*x_hi) 
          << Form(" %18.2f |",100*x_lo) 
          << Form(" %18.2f |",y_nom) 
          << Form(" %18.2f |",100*y_hi) 
          << Form(" %18.2f |",100*y_lo) 
          << endl;
      } 

      return;
    }

    void SetStyle(){
      gStyle->SetOptStat(0);

      return;
    }


  ///////////////////////////////////////////////////////////////////////////////////////////////
  void     PrintYields(string type) {
    vector<double> yieldsV;
    yieldsV.resize(nCategories);
    for (unsigned int i=0; i<yieldsV.size(); i++) yieldsV[i]=0.0;

    TIterator* iter = channelCat->typeIterator() ;
    RooCatType* tt = NULL;    
    
    int index=-1;
    while((tt=(RooCatType*) iter->Next()) ){
      index++;
      RooAbsPdf  *pdftmp  = pdf->getPdf( tt->GetName() );
      RooArgSet  *obstmp  = pdftmp->getObservables( *mc->GetObservables() ) ;
      RooRealVar *myobs   = ((RooRealVar*) obstmp->first());
      TString catName=tt->GetName();        

      double yields=0;
      if (type=="Background") {
	SetPOI(0.0);
	yields=pdftmp->expectedEvents(*myobs);
	SetPOI(1.0);
      } else if (type=="Signal") {
	SetPOI(1.0);
	yields=pdftmp->expectedEvents(*myobs);
	SetPOI(0.0);
	yields-=pdftmp->expectedEvents(*myobs);
	SetPOI(1.0);
      } else {
	RooRealSumPdf* pdfmodel=getModelPDF(pdftmp,catName);
	RooArgList funcList =  pdfmodel->funcList();
	RooLinkedListIter funcIter = funcList.iterator() ;
	RooProduct* comp = 0;
	yields=-999;
	while( (comp = (RooProduct*) funcIter.Next())) { 
	  TString compName=comp->GetName();
	  TString newName( compName(0,compName.Index("_"+catName)) );
	  if ( newName.Contains(type) ) {
	    yields=(comp->createIntegral(*myobs))->getVal()*(float)nBins[index];
	  }
	}
	
      }
      yieldsV[index]=yields;
    }

    cout << " --> Yields for " << type << endl;
    for (unsigned int i=0; i<yieldsV.size(); i++) {
      if ((float)yieldsV[i]>-400) {
	cout << Form("  %9.1f  |",(float)yieldsV[i]);
      } else {
	cout << "         --  |";
      }
    }
    cout << endl;
    return;
  }


   ///////////////////////////////////////////////////////////////////////////////////////////////
  void PrintSystematics(string type) {
    //vector<string> sysNames;
    sysNames.resize(nNP);
    vector< vector<double> > yieldsSys;
    for (int iS=0; iS<nNP; iS++) {
      vector<double> sysV;
      sysV.resize(nCategories);
      for (unsigned int i=0; i<sysV.size(); i++) sysV[i]=0.0;
      yieldsSys.push_back(sysV);
    }
    vector< double > yieldsV;
    yieldsV.resize(nCategories);
    for (unsigned int i=0; i< yieldsV.size(); i++)  yieldsV[i]=0.0;
    
    cout << "before first loop" << endl;

    // first cache the nominal yields
    TIterator* iter = channelCat->typeIterator() ;
    RooCatType* tt = NULL;    
    int index=-1;
    while((tt=(RooCatType*) iter->Next()) ){
      index++;
      RooAbsPdf  *pdftmp  = pdf->getPdf( tt->GetName() );
      RooArgSet  *obstmp  = pdftmp->getObservables( *mc->GetObservables() ) ;
      RooRealVar *myobs   = ((RooRealVar*) obstmp->first());
      double yields=0;
      if (type=="Background") {
	SetPOI(0.0);
	yields=pdftmp->expectedEvents(*myobs);
	SetPOI(1.0);
      } else if (type=="Signal") {
	SetPOI(1.0);
	yields=pdftmp->expectedEvents(*myobs);
	SetPOI(0.0);
	yields-=pdftmp->expectedEvents(*myobs);
	SetPOI(1.0);
      }
      yieldsV[index]=yields;
    }

    cout << "before sys loop" << endl;

    /// sys lists
    iter = channelCat->typeIterator() ;
    index=-1;
    while((tt=(RooCatType*) iter->Next()) ){
      index++;
      int inuis=-1;
      RooAbsPdf  *pdftmp  = pdf->getPdf( tt->GetName() );
      RooArgSet  *obstmp  = pdftmp->getObservables( *mc->GetObservables() ) ;
      RooRealVar *myobs   = ((RooRealVar*) obstmp->first());
      //cout << pdftmp << " " << obstmp << " " << myobs << endl;

      RooArgSet nuis = *mc->GetNuisanceParameters();
      TIterator* itr = nuis.createIterator();
      RooRealVar* arg;
      while ((arg=(RooRealVar*)itr->Next())) {
        if (!arg) continue;
        inuis++;
	//cout << "inuis: " << inuis << " --> " << arg->GetName() << endl;
	string tmpName=arg->GetName();
	if (tmpName.find("gamma")!=string::npos) break;
	sysNames[inuis]=arg->GetName();
	arg->setVal(1.0);
	double yields=0;
	if (type=="Background") {
	  SetPOI(0.0);
	  yields=pdftmp->expectedEvents(*myobs);
	  SetPOI(1.0);
	} else if (type=="Signal") {
	  SetPOI(1.0);
	  yields=pdftmp->expectedEvents(*myobs);
	  SetPOI(0.0);
	  yields-=pdftmp->expectedEvents(*myobs);
	  SetPOI(1.0);
	}
	arg->setVal(0.0);
	yieldsSys[inuis][index]=yields;
      }
    }

    cout << " --> Sys effects for " << type << endl;
    for (int is=0; is<nNP; is++) {
      bool doPrint=false;
      for (unsigned int i=0; i<yieldsV.size(); i++) {
	float val=(yieldsSys[is][i]-yieldsV[i])/yieldsV[i]*100;      
	if ( fabs(val)>0.1 ) {
	  doPrint=true;
	  break;
	}
      }

      if (!doPrint) continue;
      cout << Form(" %-40s | ", (sysNames[is]).c_str() );
      for (unsigned int i=0; i<yieldsV.size(); i++) {
	float val=(yieldsSys[is][i]-yieldsV[i])/yieldsV[i]*100;
	if (val==0)  cout << "     --   |";
	else         cout << Form(" %8.1f |",val); 
      }
      cout << endl;
    }
    return;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void Initialize(const char* infile , const char* outputdir, const char* workspaceName, const char* modelConfigName, const char* ObsDataName) {
    
    cout << endl << endl << "=====================================================================================================" << endl << endl;
    
    RooMsgService::instance().setGlobalKillBelow(ERROR);
    // Cosmetics
    SetStyle();
    
    // Container for the plots
    OutputDir = (TString) outputdir;
    gSystem->Exec("mkdir -p " + OutputDir+"Checks");
    outputfile = new TFile(OutputDir+"Checks/OutPutChecks.root","RECREATE");
    
    // Load workspace, model and data
    TFile *file = TFile::Open(infile);
    if (!file) {
      cout << "The file " << infile << " is not found/created, will stop here." << endl;
      return;
    } else if (verbose) 
      cout << endl << ">>>>>>>>>> SUCCESSFULLY opened file: " << infile << " <<<<<<<<<<" <<endl << endl;
    
    w      = (RooWorkspace*) file->Get(workspaceName);
    if(!w){
      cout <<"workspace: " << workspaceName << " not found" << endl;
      return;
    } else if (verbose) 
      cout << ">>>>>>>>>> SUCCESSFULLY retrieved WS: " <<  workspaceName << " <<<<<<<<<<" <<endl << endl;
    w->SetName("w");
    w->SetTitle("w");
    //w->Print();
    
    //FixUnwantedNF();

    // Activate binned likelihood calculation for binned models
    if(isBinned){
      RooFIter iter = w->components().fwdIterator() ;
      RooAbsArg* arg ;
      while((arg = iter.next())) {
	if (arg->IsA() == RooRealSumPdf::Class()) {
	  arg->setAttribute("BinnedLikelihood");
	}
      }
    }

    mc     = (ModelConfig*) w->obj(modelConfigName);
    if (!mc) {
      cout <<"model config: " << modelConfigName << " not found" << endl;
      return;
    } else if (verbose) 
      cout << ">>>>>>>>>> SUCCESSFULLY retrieved ModelConfig: " <<  modelConfigName << " <<<<<<<<<<" <<endl << endl;
    
    pdf = (RooSimultaneous*)(mc->GetPdf());
    if (!pdf) {
      cout <<"pdf not found" << endl;
      return;
    } else if (verbose) 
      cout << ">>>>>>>>>> SUCCESSFULLY retrieved PDF <<<<<<<<<<" <<endl << endl;
    channelCat = (RooCategory*) (&pdf->indexCat());

    data   = w->data(ObsDataName);
    if (!data) {
      cout <<"Data not found" << endl;
      //return;
    } else if (verbose) 
      cout << ">>>>>>>>>> SUCCESSFULLY retrieved Data <<<<<<<<<<" <<endl << endl;

    if (!data || isBlind>0){
      //if (isBlind == 1) data = makeAsimovData(mu_asimov);
      //else              data = makeAsimovData(mu_asimov, true); //fluctuated
    }
  
    // save snapshot before any fit has been done
    RooArgSet* params = (RooArgSet*) pdf->getParameters(*data) ;
    if(!w->loadSnapshot("NominalParamValues"))  
      w->saveSnapshot("NominalParamValues",*params);
    else 
      cout << " Snapshot 'NominalParamValues' already exists in  workspace, will not overwrite it" << endl;
    firstPOI= dynamic_cast<RooRealVar*>(mc->GetParametersOfInterest()->first());
    
    // Some sanity checks on the workspace
    if ( !IsChannelNameOK()   ) return;
    
    GetNominalValueNuisancePara();
    AllNPafterEachFit_vec.clear();
    AllFitResults_map.clear();
    //////////////////////////////////////////////////////////////////////////////////////////
    // Print some information
    if (!isGG) PrintModelObservables();
    
    PrintNuisanceParameters();

    cout << endl << endl << endl;
    //exit(-1);
    FillSysInfo();    
    cout << "DONE WITH PREPARESYSINFO" << endl << endl;


    PrintSysPerSample("Zprime");
    PrintSysPerSample("ttbar");
    PrintSysPerSample("multijet");

    return; /// VALERIO BREAK!!!!!!!!
    cout << endl << endl << endl;

    ///////////////////////////////////////////////////////////////////////
    TIterator* iter = channelCat->typeIterator() ;
    RooCatType* tt = NULL;   
    int cat=0;
    //w->var("normNF_LUMI")->setVal( 6.0/3.2 );
    while((tt=(RooCatType*) iter->Next()) ) {
      cout << " ---> ANALYSIS region: " << tt->GetName() << endl;
      RooAbsPdf  *pdfReg  = pdf->getPdf( tt->GetName() );
      //pdfReg->Print("v");
      RooArgSet  *obstmp  = pdfReg->getObservables( *mc->GetObservables() ) ;
      obstmp->Print("v");
      RooRealSumPdf* pdfmodel=getModelPDF(pdfReg,tt->GetName());
      pdfmodel->Print("v");
      RooArgList funcList =  pdfmodel->funcList();
      funcList.Print("v");
      RooLinkedListIter funcIter = funcList.iterator() ;
      RooProduct* comp = 0;
      while( (comp = (RooProduct*) funcIter.Next())) { 
	TString compName=comp->GetName();
	float yields=(comp->createIntegral(*obstmp))->getVal();
	cout << "  testVale:  for " << compName << " : " << yields << endl; 
      }
  
      break;

      /*
      RooAbsData *datatmp = data->reduce(Form("%s==%s::%s",channelCat->GetName(),channelCat->GetName(),tt->GetName()));
      float DataEvents=datatmp->sumEntries();
      cout << "Region: " << tt->GetName() << "  has: " 
	   << Form("%6.0f",(float)datatmp->sumEntries()) 
	   << "  with relative uncertainty: " 
	   << Form(" %7.3f ",1/sqrt(DataEvents)*100) << " % " << endl;
     
      RooAbsPdf  *pdftmp  = pdf->getPdf( tt->GetName() );
      RooArgSet  *obstmp  = pdftmp->getObservables( *mc->GetObservables() ) ;
      RooRealVar *obs     = ((RooRealVar*) obstmp->first());
      RooPlot* frame = obs->frame();
      datatmp->plotOn(frame,Name("Data"));
      for (int iB=0; iB<frame->GetNbinsX(); iB++) {
	double xP=0;
	double yP=0;
	frame->getHist()->GetPoint(iB,xP,yP);
	if (yP!=0) {
	  dataStat->SetBinContent(iB+1,cat+1,1/sqrt(yP));
	}
      }
      
      SetPOI(1.0);	    
      ///// central
      string newNameH="central"+string(tt->GetName());
      TH1F* centralH=new TH1F( newNameH.c_str(),
			       newNameH.c_str(),
			       frame->GetNbinsX(),0,1);
      RooPlot* frame2 = obs->frame();
      float postFitIntegral = pdftmp->expectedEvents(*obs);
      pdftmp->plotOn(frame2,LineWidth(2),Normalization(postFitIntegral,RooAbsReal::NumEvent),Name("Central_MC"));
      double prevVal=-999;
      int EffBin=0;
      for (int iB=0; iB<frame2->getCurve("Central_MC")->GetN(); iB++) {
	double xP=0;
	double yP=0;
	frame2->getCurve("Central_MC")->GetPoint(iB,xP,yP);
	if (xP<=0) continue;
	if (xP>=1) continue;
	if (yP!=prevVal) {
	  EffBin++;
	  centralH->SetBinContent(EffBin,yP);
	  prevVal=yP;
	}
      }


      // shifting all the gammas up
      SetAllStatErrorToSigma(1);
      postFitIntegral = pdftmp->expectedEvents(*obs);
      pdftmp->plotOn(frame2,LineWidth(2),Normalization(postFitIntegral,RooAbsReal::NumEvent),Name("StatUp"));
      prevVal=-999;
      EffBin=0;
      for (int iB=0; iB<frame2->getCurve("StatUp")->GetN(); iB++) {
	double xP=0;
	double yP=0;
	frame2->getCurve("StatUp")->GetPoint(iB,xP,yP);
	if (xP<=0) continue;
	if (xP>=1) continue;
	if (yP!=prevVal) {
	  EffBin++;
	  centralH->SetBinError(EffBin,yP-centralH->GetBinContent(EffBin));
	  if (centralH->GetBinContent(EffBin)!=0) mcStat->SetBinContent(EffBin, cat+1,centralH->GetBinError(EffBin)/centralH->GetBinContent(EffBin));
	  prevVal=yP;
	  
	}
      }
      SetAllStatErrorToSigma(0);
      cat++;
      */
    }
    return;
    cout << endl;

    outputfile->cd();
    dataStat->Write();
    mcStat->Write();

    outputfile->Close();
    

    
    //PrintYields("ttH125");
    //PrintYields("ttbar-light");
    //PrintYields("QCD");
    //PrintYields("Wjets");
    //PrintYields("Background");
    //PrintYields("Signal");
    
    cout << endl;
    //PrintSystematics("Background");
    //PrintSystematics("Signal");
    //PrintSubChannels();
    
    iter = channelCat->typeIterator() ;
    tt = NULL;    
    while((tt=(RooCatType*) iter->Next()) ){
      cout << " -- On category " << tt->GetName() << " " << endl;
      TString tmpName=tt->GetName();      
      
      RooAbsPdf  *pdftmp  = pdf->getPdf(tt->GetName()) ;
      RooArgSet  *obstmp  = pdftmp->getObservables( *mc->GetObservables() ) ;
      RooRealVar *obs     = ((RooRealVar*) obstmp->first());
      //pdftmp->getComponents()->Print();
      RooArgSet* tmpSet=pdftmp->getComponents(); 
      iter = tmpSet->createIterator() ;
      RooAbsArg* MyObs2 = NULL;
      while( (MyObs2 = (RooAbsArg*) iter->Next()) ) {
	TString obsNam=MyObs2->GetName();
	if ( obsNam.Contains("alpha") )      continue;
	if ( obsNam.Contains("gamma") )      continue;
	if ( obsNam.Contains("Constraint") ) continue;
	if ( !obsNam.Contains("model") )     continue;
	
	if ( !obsNam.Contains("epsilon") )     continue;
	cout << obsNam << " -------------------------------------------------------- " << endl;

	
	if ( !obsNam.Contains("ttbarV") )     continue;
	((RooAbsReal*)MyObs2)->Print("verbose");
	
	cout << endl; 
	cout << endl;
	ostringstream ValeTest;
	MyObs2->printArgs(ValeTest);
	//cout << ValeTest << endl;

	// this is crazy, I easily see the imformation I need but there is no way to get it out????


	//((FlexibleInterpVar*)MyObs2->Find(MyObs2))->Print("verbose");

      }
      exit(-1);
    }

    exit(-1);

    // Prepare the directory structure of the outputfile
    //MainDirSyst              = (TDirectory*) outputfile->mkdir("PlotsBeforeFit");
    //MainDirFitGlobal         = (TDirectory*) outputfile->mkdir("PlotsAfterGlobalFit");
   
    gROOT->cd(); 
  }
}    
