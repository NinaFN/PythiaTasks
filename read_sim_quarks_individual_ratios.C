#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

void read_sim_quarks_individual_ratios() {
    TFile *f = TFile::Open("sim_tuples_quarks_individual_save.root","READ");
    
    
    vector<double> *binLuminocity;
    vector<double> *binEvCounts;
    vector<double> *weightSums;
    vector<double> *sigmaGens;

    vector<TNtuple*> qeTuples(8); //should be count of bins but 7 is easier
    vector<TNtuple*> qmTuples(8);

    f->GetObject("luminocities",binLuminocity);
    f->GetObject("eventCounts",binEvCounts);
    f->GetObject("weightSums",weightSums);
    f->GetObject("sigmaGens",sigmaGens);

    int binCount = 0;

    //Create Hists
    TH1F *qePtTotal = new TH1F("qe_full","Events that produce c#bar{c} or b#bar{b} pairs in Pythia HardQCD processes;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 35, 0.0, 120.0);
    TH1F *qePtPart = new TH1F("qe_pt_part","", 35, 0.0, 120.0);
    TH1F *qePtCheckEta = new TH1F("qe_pt_check_eta",";#hat{p}_{T} (GeV/c);Ratio", 35, 0.0, 120.0);
    TH1F *qePtCheckDec = new TH1F("qe_pt_check_dec","", 35, 0.0, 120.0);
    TH1F *qePtCheckDecEta = new TH1F("qe_pt_check_dec_eta","", 35, 0.0, 120.0);
    
    TH1F *qmPtTotal = new TH1F("qm_full","Events that produce c#bar{c} or b#bar{b} pairs in Pythia HardQCD processes;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 35, 0.0, 120.0);
    TH1F *qmPtPart = new TH1F("qm_pt_part","", 35, 0.0, 120.0);
    TH1F *qmPtCheckEta = new TH1F("qm_pt_check_eta",";#hat{p}_{T} (GeV/c);Ratio", 35, 0.0, 120.0);
    TH1F *qmPtCheckDec = new TH1F("qm_pt_check_dec","", 35, 0.0, 120.0);
    TH1F *qmPtCheckDecEta = new TH1F("qm_pt_check_dec_eta","", 35, 0.0, 120.0);
    

    for(std::vector<double>::iterator it = binLuminocity->begin(); it != binLuminocity->end(); ++it){
        
        //Fill Tuples
        qeTuples[binCount] = (TNtuple*)f->Get(Form("qe%d", binCount));
        qmTuples[binCount] = (TNtuple*)f->Get(Form("qm%d", binCount));

        ////Fill Histograms

        //Electron Decays
      
        qePtPart->Reset();
        qeTuples[binCount]->Draw("ptHat>>qe_pt_part");
        qePtPart->Scale(1/(*it),"width");
        qePtTotal->Add(qePtPart);

        qePtPart->Reset();
        qeTuples[binCount]->Draw("ptHat>>qe_pt_part","etaQ<=1.5 && etaQ>=-1.5");
        qePtPart->Scale(1/(*it), "width");
        qePtCheckEta->Add(qePtPart);

        qePtPart->Reset();
        qeTuples[binCount]->Draw("ptHat>>qe_pt_part","decayFlag==1");
        qePtPart->Scale(1/(*it), "width");
        qePtCheckDec->Add(qePtPart);

        qePtPart->Reset();
        qeTuples[binCount]->Draw("ptHat>>qe_pt_part","decayFlag==1 && etaDecFlag==1");
        qePtPart->Scale(1/(*it), "width");
        qePtCheckDecEta->Add(qePtPart);

        //Muon Decays
      
        qmPtPart->Reset();
        qmTuples[binCount]->Draw("ptHat>>qm_pt_part");
        qmPtPart->Scale(1/(*it),"width");
        qmPtTotal->Add(qmPtPart);

        qmPtPart->Reset();
        qmTuples[binCount]->Draw("ptHat>>qm_pt_part","etaQ<=-2 && etaQ>=-5");
        qmPtPart->Scale(1/(*it), "width");
        qmPtCheckEta->Add(qmPtPart);

        qmPtPart->Reset();
        qmTuples[binCount]->Draw("ptHat>>qm_pt_part","decayFlag==1");
        qmPtPart->Scale(1/(*it), "width");
        qmPtCheckDec->Add(qmPtPart);

        qmPtPart->Reset();
        qmTuples[binCount]->Draw("ptHat>>qm_pt_part","decayFlag==1 && etaDecFlag==1");
        qmPtPart->Scale(1/(*it), "width");
        qmPtCheckDecEta->Add(qmPtPart);

        binCount++;
    }

    ////Plotting
    // qes
    TFile *outf =  new TFile("sim_hists_individual_ratios.root", "RECREATE");
    TCanvas *canvasQE = new TCanvas("qe_sigma","qe_sigma");
    
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetLogy();
    pad1->SetBottomMargin(0);

    qePtTotal->GetYaxis()->SetTitleOffset(1);

    pad1->Draw();
    ///pad1->SetGridx();
    pad1->cd(); 

    qePtTotal->SetLineColor(kBlack);
    qePtTotal->SetStats(0);
    //qePtTotal->GetYaxis()->SetLabelSize(0.);
    qePtTotal->Draw("SAME");

    qePtCheckEta->SetLineColor(kRed);
    qePtCheckEta->SetStats(0);
    qePtCheckEta->DrawCopy("SAME");

    qePtCheckDec->SetLineColor(kGreen);
    qePtCheckDec->SetStats(0);
    qePtCheckDec->DrawCopy("SAME");

    qePtCheckDecEta->SetLineColor(kBlue);
    qePtCheckDecEta->SetStats(0);
    qePtCheckDecEta->DrawCopy("SAME");

    auto legendqe = new TLegend();
    legendqe->AddEntry(qePtTotal,"All","l");
    legendqe->AddEntry(qePtCheckEta,"Quark produced with -1.5<#eta_{quark}<1.5","l");
    legendqe->AddEntry(qePtCheckDec,"Quark produced, decays to e","l");
    legendqe->AddEntry(qePtCheckDecEta,"Quark produced, decays to e with -0.9<#eta_{e}<0.9","l");
    
    legendqe->Draw("SAME");

    canvasQE->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);

    qePtCheckEta->GetXaxis()->SetTitleSize(13);
    qePtCheckEta->GetXaxis()->SetTitleOffset(1);
    qePtCheckEta->GetXaxis()->SetLabelSize(10);

    qePtCheckEta->GetYaxis()->SetTitleSize(13);
    qePtCheckEta->GetYaxis()->SetTitleOffset(1);
    qePtCheckEta->GetYaxis()->SetLabelSize(10);

    pad2->Draw("SAME");
    ///pad2->SetGridx();
    pad2->cd(); 

    qePtCheckEta->Divide(qePtTotal);
    qePtCheckEta->SetStats(0);
    qePtCheckEta->SetAxisRange(0,0.85,"Y");
    qePtCheckEta->DrawCopy();

    qePtCheckDec->Divide(qePtTotal);
    qePtCheckDec->SetStats(0);
    qePtCheckDec->DrawCopy("SAME");

    qePtCheckDecEta->Divide(qePtTotal);
    qePtCheckDecEta->SetStats(0);
    qePtCheckDecEta->DrawCopy("SAME");


    canvasQE->Write();

    ////

    // qms
    TCanvas *canvasQM = new TCanvas("qm_sigma","qm_sigma");
    
    TPad *pad1M = new TPad("pad1M", "pad1M", 0, 0.3, 1, 1.0);
    pad1M->SetLogy();
    pad1M->SetBottomMargin(0);
    qmPtTotal->GetYaxis()->SetTitleOffset(1);
    pad1M->Draw();
    ///pad1M->SetGridx();
    pad1M->cd(); 

    qmPtTotal->SetLineColor(kBlack);
    qmPtTotal->SetStats(0);
    qmPtTotal->Draw("SAME");

    qmPtCheckEta->SetLineColor(kRed);
    qmPtCheckEta->SetStats(0);
    qmPtCheckEta->DrawCopy("SAME");

    qmPtCheckDec->SetLineColor(kGreen);
    qmPtCheckDec->SetStats(0);
    qmPtCheckDec->DrawCopy("SAME");

    qmPtCheckDecEta->SetLineColor(kBlue);
    qmPtCheckDecEta->SetStats(0);
    qmPtCheckDecEta->DrawCopy("SAME");

    auto legendqm = new TLegend();
    legendqm->AddEntry(qmPtTotal,"All","l");
    legendqm->AddEntry(qmPtCheckEta,"Quark produced with -5<#eta_{quark}<-2","l");
    legendqm->AddEntry(qmPtCheckDec,"Quark produced, decays to #mu","l");
    legendqm->AddEntry(qmPtCheckDecEta,"Quark produced, decays to #mu with -4.5<#eta_{#mu}<-2.5","l");
    
    legendqm->Draw("SAME");

    canvasQM->cd();
    TPad *pad2M = new TPad("pad2M", "pad2M", 0, 0.05, 1, 0.3);

    pad2M->SetTopMargin(0);
    pad2M->SetBottomMargin(0.2);

    qmPtCheckEta->GetXaxis()->SetTitleSize(13);
    qmPtCheckEta->GetXaxis()->SetTitleOffset(1);
    qmPtCheckEta->GetXaxis()->SetLabelSize(10);

    qmPtCheckEta->GetYaxis()->SetTitleSize(13);
    qmPtCheckEta->GetYaxis()->SetTitleOffset(1);
    qmPtCheckEta->GetYaxis()->SetLabelSize(10);

    pad2M->Draw("SAME");
    ///pad2M->SetGridx();
    pad2M->cd(); 

    qmPtCheckEta->Divide(qmPtTotal);
    qmPtCheckEta->SetStats(0);
    qmPtCheckEta->SetAxisRange(0,0.35,"Y");
    qmPtCheckEta->DrawCopy();

    qmPtCheckDec->Divide(qmPtTotal);
    qmPtCheckDec->SetStats(0);
    qmPtCheckDec->DrawCopy("SAME");

    qmPtCheckDecEta->Divide(qmPtTotal);
    qmPtCheckDecEta->SetStats(0);
    qmPtCheckDecEta->DrawCopy("SAME");


    canvasQM->Write();
    
    delete outf;

}