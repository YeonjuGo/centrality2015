#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TChain.h>
#include <iostream>

TGraphAsymmErrors* getEfficiency(TTree *t,char *var,char *cut, char *preselection, int nBin=100, double binL=0, double binH=100)
{
   static int i = 1;
   TH1D *hPass = new TH1D(Form("hPass%d",i),"",nBin,binL,binH);
   TH1D *hAll = new TH1D(Form("hAll%d",i),"",nBin,binL,binH);
   t->Draw(Form("%s>>hAll%d",var,i),Form("(%s)",preselection));
   t->Draw(Form("%s>>hPass%d",var,i),Form("(%s)&&(%s)",preselection,cut));
   TGraphAsymmErrors *g = new TGraphAsymmErrors;
   g->BayesDivide(hPass,hAll);
   g->SetMarkerStyle(20);
   g->SetLineStyle(2);
   i++;
   cout <<cut<<" given "<<preselection<<" "<<hPass->GetEntries()/hAll->GetEntries()<<" "<<hAll->GetEntries()<<endl;
   return g;  
}
void plotEff2(char *infname)
{
   TFile *inf = new TFile(infname);
   TTree *tEvt = (TTree*)inf->Get("hiEvtAnalyzer/HiTree");
   TTree *tHLT = (TTree*)inf->Get("hltanalysis/HltTree");
   TTree *tSkim = (TTree*)inf->Get("skimanalysis/HltTree");
    
   tEvt->AddFriend(tHLT);
   tEvt->AddFriend(tSkim);

//   TGraphAsymmErrors *g = getEfficiency(tEvt,"hiBin","HLT_HIL1MinimumBiasHF1ANDPixel_SingleTrack_v1","HLT_HIL1MinimumBiasHF1ANDExpress_v1",100,0,200);
   TGraphAsymmErrors *g = getEfficiency(tEvt,"hiBin","pcollisionEventSelection","HLT_HIL1MinimumBiasHF1ANDExpress_v1",50,0,200);
   TGraphAsymmErrors *gHFHFplusANDminusTH0 = getEfficiency(tEvt,"hiBin","pcollisionEventSelection","HLT_HIL1HFplusANDminusTH0Express_v1",50,0,200);
   TGraphAsymmErrors *gHF2And = getEfficiency(tEvt,"hiBin","pcollisionEventSelection","HLT_HIL1MinimumBiasHF2ANDExpress_v1",50,0,200);
   
   TH1D *h= new TH1D("h","",200,0,200);
   h->SetXTitle("Centrality Bin");
   h->SetYTitle("Trigger Efficiency");
   
   TCanvas *c = new TCanvas("c","",800,800);
   h->Draw();
   g->Draw("p same");
   gHF2And->Draw("p same");
   gHF2And->SetLineColor(2);
   gHF2And->SetMarkerColor(2);
   gHF2And->SetMarkerStyle(24);
   //gHFHFplusANDminusTH0->Draw("p same");
   gHFHFplusANDminusTH0->SetLineColor(4);
   gHFHFplusANDminusTH0->SetMarkerColor(4);
   gHFHFplusANDminusTH0->SetMarkerStyle(24);
   
   TLegend *leg = new TLegend (0.44,0.26,0.9,0.53);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   leg->AddEntry(g,"HIL1MinimumBiasHF1ANDExpress_v1","pl");  
   leg->AddEntry(gHF2And,"HIL1MinimumBiasHF2ANDExpress_v1","pl");
   //leg->AddEntry(gHFHFplusANDminusTH0,"L1_HFplusANDminusTH0","pl");
   leg->Draw();
}
