void SaveHist(TString foutstring){
  const int nmassbin_fine=10;
  double massbin_fine[nmassbin_fine+1]={40,50,60,70,80,90,100,120,160,200,350};
  const int nmassbin_wide=10;
  double massbin_wide[nmassbin_wide+1]={40,50,60,70,80,90,100,120,160,200,350};
  const int nptbin_fine=43;
  double ptbin_fine[nptbin_fine+1]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,24,26,28,30,32,34,36,38,40,44,48,52,56,60,65,70,75,80,85,90,95,100};
  const int nptbin_wide=28;
  double ptbin_wide[nptbin_wide+1]={0,2,4,6,8,10,12,14,17,20,23,26,29,32,35,38,42,45,48,52,56,60,64,68,72,76,84,92,100};

  TUnfoldBinningV17 *binning_Rec=new TUnfoldBinningV17("Rec");
  binning_Rec->AddAxis("pt",nptbin_fine,ptbin_fine,false,true);
  binning_Rec->AddAxis("mass",nmassbin_fine,massbin_fine,true,true);
  binning_Rec->PrintStream(std::cout); 

  TUnfoldBinningV17 *binning_Gen=new TUnfoldBinningV17("Gen");
  binning_Gen->AddAxis("pt",nptbin_wide,ptbin_wide,false,true);
  binning_Gen->AddAxis("mass",nmassbin_wide,massbin_wide,true,true);
  binning_Gen->PrintStream(std::cout);
 
  
  TH2 *hmcGenRec=TUnfoldBinningV17::CreateHistogramOfMigrations(binning_Gen,binning_Rec,"hmcGenRec");
  TH1 *hdataRec=binning_Rec->CreateHistogram("hdataRec");
  TH2 *hdataRec_RecAxis=(TH2*)binning_Rec->CreateHistogram("hdataRec_RecAxis",true);
  TH2 *hmcRec_RecAxis=(TH2*)binning_Rec->CreateHistogram("hmcRec_RecAxis",true);
  TH2 *hmcGen_GenAxis=(TH2*)binning_Gen->CreateHistogram("hmcGen_GenAxis",true);
  TH2 *hmcGen_RecAxis=(TH2*)binning_Rec->CreateHistogram("hmcGen_RecAxis",true);
  TH1D *hpt=new TH1D("hpt","hpt",nptbin_wide,ptbin_wide);
  
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
      hpt->Fill(ptPreFSR->at(2),weightGen);
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

  TFile fout(foutstring,"recreate");
  hdataRec->Write();
  hmcGenRec->Write();
  binning_Gen->Write();
  binning_Rec->Write();
  hdataRec_RecAxis->Write();
  hmcGen_GenAxis->Write();
  hmcRec_RecAxis->Write();
  hmcGen_RecAxis->Write();
  hpt->Write();
  fout.Close();
}
TUnfoldDensityV17* GetUnfold(TString fhiststring,int regmode,double f){
  TFile fhist(fhiststring);
  TH1* hdataRec=(TH1*)fhist.Get("hdataRec");
  TH2* hmcGenRec=(TH2*)fhist.Get("hmcGenRec");
  TUnfoldBinningV17* binning_Rec=(TUnfoldBinningV17*)fhist.Get("Rec");
  TUnfoldBinningV17* binning_Gen=(TUnfoldBinningV17*)fhist.Get("Gen");

  TUnfoldDensityV17 *unfold=new TUnfoldDensityV17(hmcGenRec,TUnfoldV17::kHistMapOutputHoriz,(TUnfoldV17::ERegMode)regmode,TUnfoldV17::kEConstraintArea,TUnfoldDensityV17::kDensityModeBinWidth,binning_Gen,binning_Rec);
  unfold->SetInput(hdataRec,f);
  return unfold;
}
void SaveUnfold(TString fhiststring,TString foutstring,TUnfoldDensityV17* unfold){
  TFile fhist(fhiststring);
  TH1* hdataRec=(TH1*)fhist.Get("hdataRec");
  TH2* hmcGenRec=(TH2*)fhist.Get("hmcGenRec");
  TH2* hdataRec_RecAxis=(TH2*)fhist.Get("hdataRec_RecAxis");
  TH2* hmcGen_GenAxis=(TH2*)fhist.Get("hmcGen_GenAxis");
  TH2* hmcRec_RecAxis=(TH2*)fhist.Get("hmcRec_RecAxis");
  TH2* hmcGen_RecAxis=(TH2*)fhist.Get("hmcGen_RecAxis");

  TFile f(foutstring,"recreate");
  unfold->GetInputBinning()->Write();
  unfold->GetOutputBinning()->Write();
  unfold->GetOutput("hdataUnfold")->Write();
  hdataRec->Write();
  hmcGenRec->Write();
  hdataRec_RecAxis->Write();
  hmcGen_GenAxis->Write();
  hmcGen_RecAxis->Write();
  hmcRec_RecAxis->Write();
}
void Unfold(TString fhiststring,TString foutstring,int mode,int scanmethod,double f){
  int regmode=0;
  if(mode==-1){
    regmode=0;
  }else{
    regmode=mode;
  } 
  cout<<"[ISR2016MuonUnfold::Unfold] mode="<<mode<<" regmode="<<regmode<<endl;

  TUnfoldDensityV17 *unfold=GetUnfold(fhiststring,regmode,f);
  if(mode==0){
    cout<<"[ISR2016MuonUnfold::Unfold] No regularization"<<endl;
    unfold->DoUnfold(0);
    SaveUnfold(fhiststring,foutstring,unfold);
    return;
  }else if(mode==-1){
    cout<<"[ISR2016MuonUnfold::Unfold] private regularization"<<endl;
    const TUnfoldBinningV17* bin=unfold->GetOutputBinning();
    int istart=bin->GetGlobalBinNumber(0,200);
    int iend=bin->GetGlobalBinNumber(100-0.01,200);
    cout<<istart<<" "<<iend<<endl;
    unfold->RegularizeBins(istart,1,iend-istart+1,TUnfoldV17::kRegModeCurvature);
  }
			   
  TH2 *histL= unfold->GetL("L");
  if(histL){
    for(Int_t j=1;j<=histL->GetNbinsY();j++) {
      cout<<"L["<<unfold->GetLBinning()->GetBinName(j)<<"]";
      for(Int_t i=1;i<=histL->GetNbinsX();i++) {
	Double_t c=histL->GetBinContent(i,j);
	if(c!=0.0) cout<<" ["<<i<<"]="<<c;
      }
      cout<<"\n";
    }
  }else{
    cout<<"[ISR2016MuonUnfold::Unfold] No L"<<endl;
  }
  TCanvas* c1;
  if(scanmethod==0){
    cout<<"[ISR2016MuonUnfold::Unfold] L curve scan"<<endl;
    TGraph *lcurve;
    TSpline *logtaux,*logtauy,*logtaucurvature;
    int ibest=unfold->ScanLcurve(20,0,0,&lcurve,&logtaux,&logtauy,&logtaucurvature);
    c1=new TCanvas;
    c1->Divide(2,2);
    c1->cd(1);
    lcurve->Draw();
    TGraph* best1=new TGraph;
    best1->SetPoint(0,lcurve->GetX()[ibest],lcurve->GetY()[ibest]);
    best1->SetMarkerStyle(22);
    best1->Draw("p");

    c1->cd(2);
    logtaux->Draw();
    TGraph* best2=new TGraph;
    double x,y;
    logtaux->GetKnot(ibest,x,y);
    best2->SetPoint(0,x,y);
    best2->SetMarkerStyle(22);
    best2->Draw("p");

    c1->cd(3);
    logtauy->Draw();
    TGraph* best3=new TGraph;
    logtauy->GetKnot(ibest,x,y);
    best3->SetPoint(0,x,y);
    best3->SetMarkerStyle(22);
    best3->Draw("p");

    c1->cd(4);
    logtaucurvature->Draw();
    TGraph* best4=new TGraph;
    logtaucurvature->GetKnot(ibest,x,y);
    best4->SetPoint(0,x,y);
    best4->SetMarkerStyle(22);
    best4->Draw("p");
    cout<<"[ISR2016MuonUnfold::Unfold] ibest="<<ibest<<endl;
  }else if(scanmethod==1){
    cout<<"[ISR2016MuonUnfold::Unfold] No code for ScanTau. Please add it."<<ibest<<endl;
    return;
  }else{
    cout<<"[ISR2016MuonUnfold::Unfold] Unkown scanmethod"<<ibest<<endl;
    return;
  }    
  SaveUnfold(fhiststring,foutstring,unfold);
}
/*				      
void Unfold(TString finstring="hist.root",TString foutstring="unfold.root",TUnfoldV17::ERegMode regmode=TUnfoldV17::kRegModeCurvature,double scalebias=0.0){
  TFile fhist(finstring);
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
    Int_t iBest=unfold->ScanTau(100,0.01,10,&rhoLogTau,
				TUnfoldDensityV17::kEScanTauRhoAvg,
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
  TFile fout(foutstring,"recreate");
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
*/

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
    cout<<"[ISR2016MuonUnfold::Compare] ibin "<<firstbin<<" "<<lastbin<<endl;
    hist1=(TH1D*)h1->ProjectionX(h1->GetName()+TString("_x"),firstbin,lastbin)->Clone(h1->GetName());
    firstbin=h2->GetYaxis()->FindBin(min);
    lastbin=h2->GetYaxis()->FindBin(max-0.000001);
    cout<<"[ISR2016MuonUnfold::Compare] ibin "<<firstbin<<" "<<lastbin<<endl;
    hist2=(TH1D*)h2->ProjectionX(h2->GetName()+TString("_x"),firstbin,lastbin)->Clone(h2->GetName());
  }else if(iaxis==1){
    firstbin=h1->GetXaxis()->FindBin(min);
    lastbin=h1->GetXaxis()->FindBin(max-0.000001);
    cout<<"[ISR2016MuonUnfold::Compare] ibin "<<firstbin<<" "<<lastbin<<endl;
    hist1=(TH1D*)h1->ProjectionY(h1->GetName()+TString("_y"),firstbin,lastbin)->Clone(h1->GetName()); 
    firstbin=h2->GetXaxis()->FindBin(min);
    lastbin=h2->GetXaxis()->FindBin(max-0.000001);
    cout<<"[ISR2016MuonUnfold::Compare] ibin "<<firstbin<<" "<<lastbin<<endl;
    hist2=(TH1D*)h2->ProjectionY(h2->GetName()+TString("_y"),firstbin,lastbin)->Clone(h2->GetName());
  }
  RenormWithBinWidth(hist1);
  RenormWithBinWidth(hist2);
  hist1->SetDirectory(0);hist2->SetDirectory(0);
  TCanvas* c1=new TCanvas;  
  hist1->Draw();
  hist2->Draw("sames");
  hist2->SetLineColor(2);
  hist1->GetYaxis()->SetRangeUser(0,hist1->GetMaximum()>hist2->GetMaximum()?1.1*hist1->GetMaximum():1.1*hist2->GetMaximum());
  cout<<"[ISR2016MuonUnfold::Compare] Integral:"<<hist1->Integral("width")<<" "<<hist2->Integral("width")<<endl;
}
void Compare(TString fin,TString hist1name,TString hist2name,int iaxis,double min,double max){
  TFile funfold(fin);
  TH2* hist1=(TH2*)funfold.Get(hist1name);
  TH2* hist2=(TH2*)funfold.Get(hist2name);
  Compare(hist1,hist2,iaxis,min,max);
}

void Compare(TString hist1name,TString hist2name,int iaxis,double min,double max){
  Compare("unfold.root",hist1name,hist2name,iaxis,min,max);
}

void Compare(int iaxis,double min,double max){
  Compare("hdataUnfold","hmcGen_GenAxis",iaxis,min,max);
}

void CheckCondition(TString histfilename){
  TFile fhist(histfilename);
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
