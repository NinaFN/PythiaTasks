#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

void read_sim_quarks_individual() {
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
    TH1F *qePtTotal = new TH1F("qe_full","Heading;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} ???", 35, 0.0, 120.0);
    TH1F *qePtPart = new TH1F("qe_pt_part","", 35, 0.0, 120.0);
    TH1F *qePtCheckEta = new TH1F("qe_pt_check_eta","", 35, 0.0, 120.0);
    TH1F *qePtCheckDec = new TH1F("qe_pt_check_dec","", 35, 0.0, 120.0);
    TH1F *qePtCheckDecEta = new TH1F("qe_pt_check_dec_eta","", 35, 0.0, 120.0);
    
    TH1F *qmPtTotal = new TH1F("qm_full","Heading;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} ???", 35, 0.0, 120.0);
    TH1F *qmPtPart = new TH1F("qm_pt_part","", 35, 0.0, 120.0);
    TH1F *qmPtCheckEta = new TH1F("qm_pt_check_eta","", 35, 0.0, 120.0);
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
    TFile *outf =  new TFile("sim_hists_individual.root", "RECREATE");
    TCanvas *canvasQE = new TCanvas("qe_sigma","qe_sigma");
    
    gPad->SetLogy();
    qePtTotal->SetLineColor(kBlack);
    qePtTotal->Draw();

    qePtCheckEta->SetLineColor(kRed);
    qePtCheckEta->SetStats(0);
    qePtCheckEta->Draw("SAME");

    qePtCheckDec->SetLineColor(kGreen);
    qePtCheckDec->SetStats(0);
    qePtCheckDec->Draw("SAME");

    qePtCheckDecEta->SetLineColor(kBlue);
    qePtCheckDecEta->SetStats(0);
    qePtCheckDecEta->Draw("SAME");


    auto legendqe = new TLegend();
    legendqe->AddEntry(qePtTotal,"All","l");
    legendqe->AddEntry(qePtCheckEta,"-1.5<#eta_{quark}<1.5","l");
    legendqe->AddEntry(qePtCheckDec,"Decay to e","l");
    legendqe->AddEntry(qePtCheckDecEta,"Decay to e, -0.9<#eta_{e}<0.9","l");
    
    legendqe->Draw("SAME");

    canvasQE->Write();

    //qms
    TCanvas *canvasQM = new TCanvas("qm_sigma","qm_sigma");
    
    gPad->SetLogy();
    qmPtTotal->SetLineColor(kBlack);
    qmPtTotal->Draw();

    qmPtCheckEta->SetLineColor(kRed);
    qmPtCheckEta->SetStats(0);
    qmPtCheckEta->Draw("SAME");

    qmPtCheckDec->SetLineColor(kGreen);
    qmPtCheckDec->SetStats(0);
    qmPtCheckDec->Draw("SAME");

    qmPtCheckDecEta->SetLineColor(kBlue);
    qmPtCheckDecEta->SetStats(0);
    qmPtCheckDecEta->Draw("SAME");


    auto legendqm = new TLegend();
    legendqm->AddEntry(qmPtTotal,"All","l");
    legendqm->AddEntry(qmPtCheckEta,"-5<#eta_{quark}<-2","l");
    legendqm->AddEntry(qmPtCheckDec,"Decay to #mu","l");
    legendqm->AddEntry(qmPtCheckDecEta,"Decay to #mu, -0.9<#eta_{#mu}<0.9","l");
    
    legendqm->Draw("SAME");

    canvasQM->Write();

    delete outf;

}