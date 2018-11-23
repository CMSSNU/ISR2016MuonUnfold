void SaveBinning(){
  TFile *fout=new TFile("binning.root","recreate");
  const int nmassbin_fine=38;
  double massbin_fine[nmassbin_fine+1]={40,44,48,52,56,60,64,68,72,76,80,84,88,92,96,100,104,108,112,116,120,130,140,150,160,170,180,190,200,220,240,250,260,270,280,290,300,320,350};
  const int nmassbin_wide=26;
  double massbin_wide[nmassbin_wide+1]={40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,136,152,168,184,200,225,250,275,300,350};
  const int nptbin_fine=56;
  double ptbin_fine[nptbin_fine+1]={0,0.8,1.6,2.4,3.2,4,4.8,5.6,6.4,7.2,8,8.8,9.6,10.4,11.2,12,12.8,13.6,14.4,15.2,16,16.8,17.6,18.4,19.2,20,21.5,23,24.5,26,27.5,29,30.5,32,33.5,35,36.5,38,39.5,41,42.5,45,48,51,54,57,60,64,68,72,76,80,84,88,92,96,100};
  const int nptbin_wide=43;
  double ptbin_wide[nptbin_wide+1]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,24,26,28,30,32,34,36,38,40,44,48,52,56,60,65,70,75,80,85,90,95,100};

  TUnfoldBinningV17 *binning_Rec=new TUnfoldBinningV17("Rec");
  binning_Rec->AddAxis("pt",nptbin_fine,ptbin_fine,false,true);
  binning_Rec->AddAxis("mass",nmassbin_fine,massbin_fine,true,true);
  binning_Rec->PrintStream(std::cout); 

  TUnfoldBinningV17 *binning_Gen=new TUnfoldBinningV17("Gen");
  binning_Gen->AddAxis("pt",nptbin_wide,ptbin_wide,false,true);
  binning_Gen->AddAxis("mass",nmassbin_wide,massbin_wide,true,true);
  binning_Gen->PrintStream(std::cout); 
  
  binning_Rec->Write();
  binning_Gen->Write();
  fout->Close();
  delete fout;
 }

void SaveHist(){
  TFile fbin("binning.root");
  TUnfoldBinningV17* binning_Rec=(TUnfoldBinningV17*)fbin.Get("Rec");
  TUnfoldBinningV17* binning_Gen=(TUnfoldBinningV17*)fbin.Get("Gen");
  
  TH2 *hmcGenRec=TUnfoldBinningV17::CreateHistogramOfMigrations(binning_Gen,binning_Rec,"hmcGenRec");
  TH1 *hdataRec=binning_Rec->CreateHistogram("hdataRec");
  TH2 *hdataRec_RecAxis=(TH2*)binning_Rec->CreateHistogram("hdataRec_RecAxis",true);
  TH2 *hmcRec_RecAxis=(TH2*)binning_Rec->CreateHistogram("hmcRec_RecAxis",true);
  TH2 *hmcGen_GenAxis=(TH2*)binning_Gen->CreateHistogram("hmcGen_GenAxis",true);
  TH2 *hmcGen_RecAxis=(TH2*)binning_Rec->CreateHistogram("hmcGen_RecAxis",true);

  Double_t weightGen,weightTotal;
  int ispassRec,isfiducialPreFSR,DYtautau;
  int nentries;

  vector<Double_t> *ptPreFSR = 0;
  vector<Double_t> *mPreFSR = 0;
  vector<Double_t> *ptRec = 0;
  vector<Double_t> *mRec = 0;

  
  TFile *fdata=new TFile("ISR2016MuonAnalyzer_data_DoubleMuon_cat_v8-0-7.root");
  TTree *tdata=(TTree *)fdata->Get("tree");
  tdata->SetBranchAddress("ptRec",&ptRec);
  tdata->SetBranchAddress("mRec",&mRec);
  tdata->SetBranchAddress("ispassRec",&ispassRec);
  tdata->SetBranchAddress("weightTotal",&weightTotal);
  nentries=tdata->GetEntries();
  cout<<"start data loop"<<endl;
  for(int i=0;i<nentries;i++){
    if(i%10000000==0) cout<<i<<endl;
    tdata->GetEntry(i);
    if(ispassRec){
      hdataRec->Fill(binning_Rec->GetGlobalBinNumber(ptRec->at(2),mRec->at(2)));
      hdataRec_RecAxis->Fill(ptRec->at(2),mRec->at(2));
    }
  }

  TFile fsignal("ISR2016MuonAnalyzer_signal.root");
  TTree* tsignal=(TTree*)fsignal.Get("tree");  
  tsignal->SetBranchAddress("ptRec",&ptRec);
  tsignal->SetBranchAddress("mRec",&mRec);
  tsignal->SetBranchAddress("ispassRec",&ispassRec);
  tsignal->SetBranchAddress("isfiducialPreFSR",&isfiducialPreFSR);
  tsignal->SetBranchAddress("ptPreFSR",&ptPreFSR);
  tsignal->SetBranchAddress("mPreFSR",&mPreFSR);
  tsignal->SetBranchAddress("weightGen",&weightGen);
  tsignal->SetBranchAddress("weightTotal",&weightTotal);
  tsignal->SetBranchAddress("DYtautau",&DYtautau);
  nentries=tsignal->GetEntries();
  cout<<"start signal loop"<<endl;
  TH1F* hmGen=new TH1F("hmGen","hmGen",400,0,400);
  TH1F* hptGen=new TH1F("hptGen","hptGen",200,0,200);
  TH1F* hmRec=new TH1F("hmRec","hmRec",400,0,400);
  TH1F* hptRec=new TH1F("hptRec","hptRec",200,0,200);
  for(int i=0;i<nentries;i++){
    tsignal->GetEntry(i);
    if(i%10000000==0) cout<<i<<endl;
    if(!isfiducialPreFSR) continue;
    if(DYtautau==0){
      hmcGenRec->Fill(binning_Gen->GetGlobalBinNumber(ptPreFSR->at(2),mPreFSR->at(2)),ispassRec?binning_Rec->GetGlobalBinNumber(ptRec->at(2),mRec->at(2)):0,weightGen);
      hmcGen_GenAxis->Fill(ptPreFSR->at(2),mPreFSR->at(2),weightGen);
      hmcGen_RecAxis->Fill(ptPreFSR->at(2),mPreFSR->at(2),weightGen);
      if(ispassRec){
	hmcRec_RecAxis->Fill(ptRec->at(2),mRec->at(2),weightTotal);
      }
    }else{
      if(ispassRec){
	hdataRec->Fill(binning_Rec->GetGlobalBinNumber(ptRec->at(2),mRec->at(2)),-1.*weightTotal);
	hdataRec_RecAxis->Fill(ptRec->at(2),mRec->at(2),-1.*weightTotal);
      }
    }
  }

  TFile fbackground("ISR2016MuonAnalyzer_background.root");
  TTree* tbackground=(TTree*)fbackground.Get("tree");  
  tbackground->SetBranchAddress("ptRec",&ptRec);
  tbackground->SetBranchAddress("mRec",&mRec);
  tbackground->SetBranchAddress("ispassRec",&ispassRec);
  tbackground->SetBranchAddress("isfiducialPreFSR",&isfiducialPreFSR);
  tbackground->SetBranchAddress("weightGen",&weightGen);
  tbackground->SetBranchAddress("weightTotal",&weightTotal);
  nentries=tbackground->GetEntries();
  cout<<"start background loop"<<endl;
  for(int i=0;i<nentries;i++){
    tbackground->GetEntry(i);
    if(i%10000000==0) cout<<i<<endl;
    if(!ispassRec) continue;
    hdataRec->Fill(binning_Rec->GetGlobalBinNumber(ptRec->at(2),mRec->at(2)),-1.*weightTotal);
    hdataRec_RecAxis->Fill(ptRec->at(2),mRec->at(2),-1.*weightTotal);
  }

  TFile fout("hist.root","recreate");
  hdataRec->Write();
  hmcGenRec->Write();
  binning_Gen->Write();
  binning_Rec->Write();
  hdataRec_RecAxis->Write();
  hmcGen_GenAxis->Write();
  hmcRec_RecAxis->Write();
  hmcGen_RecAxis->Write();
  fout.Close();
}

void Unfold(TUnfoldV17::ERegMode regmode=TUnfoldV17::kRegModeCurvature,double scalebias=0.0){
  TFile fhist("hist.root");
  TH1* hdataRec=(TH1*)fhist.Get("hdataRec");
  TH2* hmcGenRec=(TH2*)fhist.Get("hmcGenRec");
  TH2* hdataRec_RecAxis=(TH2*)fhist.Get("hdataRec_RecAxis");
  TH2* hmcGen_GenAxis=(TH2*)fhist.Get("hmcGen_GenAxis");
  TH2* hmcRec_RecAxis=(TH2*)fhist.Get("hmcRec_RecAxis");
  TH2* hmcGen_RecAxis=(TH2*)fhist.Get("hmcGen_RecAxis");
  TUnfoldBinningV17* binning_Rec=(TUnfoldBinningV17*)fhist.Get("Rec");
  TUnfoldBinningV17* binning_Gen=(TUnfoldBinningV17*)fhist.Get("Gen");

  TCanvas* canvas=new TCanvas;
  TUnfoldDensityV17 *unfold;
  if(regmode==TUnfoldV17::kRegModeNone){
    unfold=new TUnfoldDensityV17(hmcGenRec,TUnfoldV17::kHistMapOutputHoriz,regmode,TUnfoldV17::kEConstraintArea,TUnfoldDensityV17::kDensityModeBinWidth,binning_Gen,binning_Rec);
    cout<<&unfold<<endl;
    unfold->SetInput(hdataRec,scalebias);
    unfold->DoUnfold(0);
  }else{
    unfold=new TUnfoldDensityV17(hmcGenRec,TUnfoldV17::kHistMapOutputHoriz,regmode,TUnfoldV17::kEConstraintArea,TUnfoldDensityV17::kDensityModeBinWidth,binning_Gen,binning_Rec);
    cout<<&unfold<<endl;
    unfold->SetInput(hdataRec,scalebias);
    
    TH2 *histL= unfold->GetL("L");
    for(Int_t j=1;j<=histL->GetNbinsY();j++) {
      cout<<"L["<<unfold->GetLBinning()->GetBinName(j)<<"]";
      for(Int_t i=1;i<=histL->GetNbinsX();i++) {
	Double_t c=histL->GetBinContent(i,j);
	if(c!=0.0) cout<<" ["<<i<<"]="<<c;
      }
      cout<<"\n";
    }
    
    Int_t nScan=30;
    TSpline *rhoLogTau=0;
    TGraph *lCurve=0;
    // for determining tau, scan the correlation coefficients
    // correlation coefficients may be probed for all distributions
    // or only for selected distributions
    // underflow/overflow bins may be included/excluded
    //
    const char *SCAN_AXISSTEERING=0;
    Int_t iBest=unfold->ScanTau(nScan,0.,0.,&rhoLogTau,
				TUnfoldDensity::kEScanTauRhoMax,
				0,0,
				&lCurve);
    // create graphs with one point to visualize best choice of tau
    Double_t t[1],rho[1],x[1],y[1];
    rhoLogTau->GetKnot(iBest,t[0],rho[0]);
    lCurve->GetPoint(iBest,x[0],y[0]);
    TGraph *bestRhoLogTau=new TGraph(1,t,rho);
    TGraph *bestLCurve=new TGraph(1,x,y);
    Double_t *tAll=new Double_t[nScan],*rhoAll=new Double_t[nScan];
    for(Int_t i=0;i<nScan;i++) {
      rhoLogTau->GetKnot(i,tAll[i],rhoAll[i]);
    }
    TGraph *knots=new TGraph(nScan,tAll,rhoAll);
    cout<<"chi**2="<<unfold->GetChi2A()<<"+"<<unfold->GetChi2L()
	<<" / "<<unfold->GetNdf()<<"\n";
    
    canvas->Divide(1,2);
    // (4) scan of correlation vs tau
    canvas->cd(1);
    rhoLogTau->Draw();
    knots->Draw("*");
    bestRhoLogTau->SetMarkerColor(kRed);
    bestRhoLogTau->Draw("*");
    // (6) L-curve
    canvas->cd(2);
    lCurve->Draw("AL");
    bestLCurve->SetMarkerColor(kRed);
    bestLCurve->Draw("*");
  }
  TFile fout("unfold.root","recreate");
  TH1* hdataUnfold=unfold->GetOutput("hdataUnfold");
  hdataUnfold->Write();
  hdataRec->Write();
  hmcGenRec->Write();
  hdataRec_RecAxis->Write();
  hmcGen_GenAxis->Write();
  hmcRec_RecAxis->Write();
  hmcGen_RecAxis->Write();
  binning_Gen->Write();
  binning_Rec->Write();
  canvas->Write();
    
  fout.Close();
}

void RenormWithBinWidth(TH1* hist){
  for(int i=1;i<=hist->GetNbinsX();i++){
    double y=hist->GetBinContent(i);
    double ey=hist->GetBinError(i);
    double width=hist->GetBinWidth(i);
    hist->SetBinContent(i,y/width);
    hist->SetBinError(i,ey/width);
  }
}
void Compare(TH2* h1,TH2* h2,int iaxis,double min, double max){
  int firstbin,lastbin;
  TH1 *hist1,*hist2;
  if(iaxis==0){
    firstbin=h1->GetYaxis()->FindBin(min);
    lastbin=h1->GetYaxis()->FindBin(max-0.000001);
    cout<<"ibin "<<firstbin<<" "<<lastbin<<endl;
    hist1=h1->ProjectionX(h1->GetName()+TString("_x"),firstbin,lastbin);
    firstbin=h2->GetYaxis()->FindBin(min);
    lastbin=h2->GetYaxis()->FindBin(max-0.000001);
    cout<<"ibin "<<firstbin<<" "<<lastbin<<endl;
    hist2=h2->ProjectionX(h2->GetName()+TString("_x"),firstbin,lastbin);
  }else if(iaxis==1){
    firstbin=h1->GetXaxis()->FindBin(min);
    lastbin=h1->GetXaxis()->FindBin(max-0.000001);
    cout<<"ibin "<<firstbin<<" "<<lastbin<<endl;
    hist1=h1->ProjectionY(h1->GetName()+TString("_y"),firstbin,lastbin); 
    firstbin=h2->GetXaxis()->FindBin(min);
    lastbin=h2->GetXaxis()->FindBin(max-0.000001);
    cout<<"ibin "<<firstbin<<" "<<lastbin<<endl;
    hist2=h2->ProjectionY(h2->GetName()+TString("_y"),firstbin,lastbin);
  }
  RenormWithBinWidth(hist1);
  RenormWithBinWidth(hist2);
  TCanvas* c1=new TCanvas;  
  hist1->Draw();
  hist2->Draw("sames");
  hist2->SetLineColor(2);
  hist1->GetYaxis()->SetRangeUser(0,hist1->GetMaximum()>hist2->GetMaximum()?1.1*hist1->GetMaximum():1.1*hist2->GetMaximum());
  cout<<hist1->Integral("width")<<" "<<hist2->Integral("width")<<endl;
}
void Compare(TString fin,TString hist1name,TString hist2name,int iaxis,double min,double max){
  TFile funfold(fin);
  TH2* hist1=(TH2*)funfold.Get(hist1name);
  TH2* hist2=(TH2*)funfold.Get(hist2name);
  hist1->AddDirectory(0);hist2->AddDirectory(0);
  Compare(hist1,hist2,iaxis,min,max);
}
void Compare(TString hist1name,TString hist2name,int iaxis,double min,double max){
  TFile funfold("unfold.root");
  TH2* hist1=(TH2*)funfold.Get(hist1name);
  TH2* hist2=(TH2*)funfold.Get(hist2name);
  hist1->AddDirectory(0);hist2->AddDirectory(0);
  Compare(hist1,hist2,iaxis,min,max);
}

void Compare(int iaxis,double min,double max){
  Compare("hdataUnfold","hmcGen_GenAxis",iaxis,min,max);
}

void CheckCondition(){
  TFile fhist("hist.root");
  TH2* hmcGenRec=(TH2*)fhist.Get("hmcGenRec");
  TUnfoldBinningV17* binning_Rec=(TUnfoldBinningV17*)fhist.Get("Rec");
  TUnfoldBinningV17* binning_Gen=(TUnfoldBinningV17*)fhist.Get("Gen");

  TUnfoldDensityV17 *unfold;
  unfold=new TUnfoldDensityV17(hmcGenRec,TUnfoldV17::kHistMapOutputHoriz,TUnfoldV17::kRegModeNone,TUnfoldV17::kEConstraintArea,TUnfoldDensityV17::kDensityModeBinWidth,binning_Gen,binning_Rec);
  TH2D* hprob=(TH2D*)unfold->GetProbabilityMatrix("hprob");
  TMatrixD* mprob=new TMatrixD(hprob->GetNbinsX()+2,hprob->GetNbinsY()+2,hprob->GetArray(),"F");
  mprob->Print();
  TMatrixD* mtprob=new TMatrixD(hprob->GetNbinsY()+2,hprob->GetNbinsX()+2);
  mtprob->Transpose(*mprob);
  mtprob->Print();
  TDecompSVD* svdprob=new TDecompSVD(*mtprob);
  cout<<svdprob->Condition()<<endl;
  TFile f("matrix.root","recreate");
  svdprob->Write();
  svdprob->Print();
}
