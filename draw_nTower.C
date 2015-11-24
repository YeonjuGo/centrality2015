// Author Yeonju Go
// last modification : 2015/01/27 
// 
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TProfile.h"
#include "TPad.h"
#include "TLatex.h"
#include "stdio.h"
#include "../../yjUtility.h"
const int rebinN = 5;
float normHistHere(TH1* hNom=0, TH1* hDen=0, TH1* hRatio=0, double cut_i=700, double cut_f=900);

void draw_nTower(){
    SetHistTitleStyle();
    SetyjPadStyle();
    gStyle -> SetOptStat(0);

    double thr[] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
    //double thr[] = {0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0};
    int nbin = sizeof(thr)/sizeof(double);
    double thrArr[nbin+1];
    for(int i=0;i<nbin+1;i++){
        double binwidth = (thr[1]-thr[0])/2.;
        if(i!=nbin) thrArr[i] = thr[i]-binwidth;
        else thrArr[i] = thr[nbin-1]+binwidth;
    }
#if 1    
    TH1D* h1D_e = new TH1D("h1D_e_ratio",";(# of HF tower) E Threshold (GeV);DATA/MC Ratio", nbin,thrArr);
    for(int i=0;i<nbin;i++){
        cout << "START transverse energy threshold : " << thr[i] << endl;
        TFile* f = new TFile(Form("histfiles/nTower_eThr%.1f_etThr0.0.root",thr[i]));
        TH1F* hdata = (TH1F*) f->Get("h1F_sample0");
        TH1F* hmc = (TH1F*) f->Get("h1F_sample1");
        hdata->GetYaxis()->SetRangeUser(0.1,100000);
        hmc->GetYaxis()->SetRangeUser(0.1,100000);
        TH1F* hratio = (TH1F*) hdata->Clone("h1F_ratio");
        float norm_i = 300;
        float norm_f = 700;
        float ratio = normHistHere(hdata, hmc, hratio, norm_i, norm_f);
        h1D_e->SetBinContent(i+1,ratio);

        //cout << ratio << endl;

        TCanvas *c = new TCanvas("c","c", 500,1000);
        c->Divide(1,2);
        c->cd(1);
        gPad->SetLogy();
        hdata->DrawCopy();
        hmc->DrawCopy("hist same");
        //drawText(Form("tower et > 0.0"), 0.46,0.80);
        drawText(Form("tower e > %.1f",thr[i]), 0.46,0.88);
        Double_t range = cleverRange(hdata,hmc);
        jumSun(norm_i, 0.1, norm_i, range);
        jumSun(norm_f, 0.1, norm_f, range);
        c->cd(2);    
        hratio->Rebin(rebinN);
        hratio->Scale(1./rebinN); 
        hratio->GetYaxis()->SetTitle("DATA/MC");
        hratio->GetYaxis()->SetRangeUser(0,2);
        hratio->DrawCopy();
        jumSun(0, 1, 1000, 1);
        drawText(Form("DATA/MC = %.3f", ratio), 0.26,0.88);
        c->SaveAs(Form("png/nTower_eThr%.1f_etThr0.0_normRange%dto%d.png",thr[i],(int)norm_i,(int)norm_f));
    }
    TCanvas* c_e = new TCanvas("c_e", "",300,300);
    h1D_e->Draw("hist");
    c_e->SaveAs("png/eThr_vs_ratio_distribution.png");

#endif
    TH1D* h1D_et = new TH1D("h1D_et_ratio",";(# of HF tower) E_{T} Threshold (GeV);DATA/MC Ratio", nbin,thrArr);
   for(int i=0;i<nbin;i++){
        cout << "START transverse energy threshold : " << thr[i] << endl;
        TFile* f = new TFile(Form("histfiles/recTower_n_eThr0.0_etThr%.1f.root",thr[i]));
        TH1F* hdata = (TH1F*) f->Get("h1F_sample0");
        TH1F* hmc = (TH1F*) f->Get("h1F_sample1");
        hdata->GetYaxis()->SetRangeUser(0.1,100000);
        hmc->GetYaxis()->SetRangeUser(0.1,100000);
        TH1F* hratio = (TH1F*) hdata->Clone("h1F_ratio");
        float norm_i = 100;
        float norm_f = 700;
        //if(i>5) { norm_i=200; norm_f=400; }
        float ratio = normHistHere(hdata, hmc, hratio, norm_i, norm_f);
        h1D_et->SetBinContent(i+1,ratio);

        TCanvas *c = new TCanvas("c","c", 500,1000);
        c->Divide(1,2);
        c->cd(1);
        gPad->SetLogy();
        hdata->DrawCopy();
        hmc->DrawCopy("hist same");
        drawText(Form("tower et > %.1f",thr[i]), 0.46,0.80);
        //drawText(Form("tower e > 0.0"), 0.46,0.88);
        Double_t range = cleverRange(hdata,hmc);
        jumSun(norm_i, 0.1, norm_i, range);
        jumSun(norm_f, 0.1, norm_f, range);
        c->cd(2);    
        hratio->Rebin(rebinN);
        hratio->Scale(1./rebinN);
        hratio->GetYaxis()->SetTitle("DATA/MC");
        hratio->GetYaxis()->SetRangeUser(0,2);
        hratio->DrawCopy();
        jumSun(0, 1, 1000, 1);
        drawText(Form("DATA/MC = %.4f", ratio), 0.26,0.88);
        c->SaveAs(Form("png/nTower_eThr0.0_etThr%.1f_normRange%dto%d.png",thr[i],(int)norm_i,(int)norm_f));
    }

    TCanvas* c_et = new TCanvas("c_e", "",300,300);
    h1D_et->Draw("hist");
    c_et->SaveAs("png/etThr_vs_ratio_distribution.png");

}


float normHistHere(TH1* hNom, TH1* hDen, TH1* hRatio, double cut_i, double cut_f){
    int cutBinFrom = hDen->FindBin(cut_i);
    int cutBinTo = hDen->FindBin(cut_f);
    int lastBin = hDen->FindBin(990);

    TH1::SetDefaultSumw2();
    hDen->Scale(hNom->Integral(cutBinFrom,cutBinTo)/hDen->Integral(cutBinFrom,cutBinTo));
    hRatio->Divide(hNom,hDen);

    //cout << "total integral of DATA in the whole range: " << hData->Integral()<< endl;
    //cout << "total integral of MC in the whole range: " << hMC->Integral()<< endl;
    //cout << "full efficiency in the whole range: " <<  hData->Integral()/hMC->Integral()<< endl;

    double nomInt = hNom->Integral(0,cutBinTo)+hNom->Integral(cutBinTo+1,lastBin);
    double denInt = hDen->Integral(0,cutBinTo)+hNom->Integral(cutBinTo+1,lastBin);
    return nomInt/denInt;
} 
