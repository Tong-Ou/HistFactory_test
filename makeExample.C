
void makeSimpleExample(){
  TFile* example = new TFile("example1.root","RECREATE");
  TH1F* data = new TH1F("data","data", 2,1,2);
  TH1F* signal = new TH1F("zprime","signal histogram (pb)", 2,1,2);
  TH1F* background1 = new TH1F("background1","background 1 histogram (pb)", 2,1,2);
  TH1F* background2 = new TH1F("background2","background 2 histogram (pb)", 2,1,2);
  TH1F* statUncert = new TH1F("background1_statUncert", "statUncert", 2,1,2);

  // run with 1 pb
  data->SetBinContent(1,122);
  data->SetBinContent(2,112);

  signal->SetBinContent(1,20);
  signal->SetBinContent(2,10);

  background1->SetBinContent(1,100);
  background1->SetBinContent(2,50);
  background2->SetBinContent(1,0);
  background2->SetBinContent(2,100);

  // A small statistical uncertainty
  statUncert->SetBinContent(1, .05);  // 5% uncertainty
  statUncert->SetBinContent(2, .05);  // 5% uncertainty


  example->Write();
  example->Close();
  //////////////////////


}

void makeExample(){
//  makeDataDriven();
//  makeShapeSys2DDataset();
//  makeShapeSysDataset();
//  makeStatErrorDataSet();
  makeSimpleExample();
}
