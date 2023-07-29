#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

void read_sim_leps_individual_cuts() {
    TFile *f = TFile::Open("sim_tuples_leps_individual.root","READ");
    
    
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
    TH1F *qePtTotal = new TH1F("qe_full","Electrons produced in Pythia HardQCD process collisions;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 35, 0.0, 120.0);
    TH1F *qePtPart = new TH1F("qe_pt_part","", 35, 0.0, 120.0);
    TH1F *qePtCut1 = new TH1F("qe_pt_cut1",";#hat{p}_{T} (GeV/c);#frac{N_{e}}{N_{e paired}}", 35, 0.0, 120.0);
    TH1F *qePtCut3 = new TH1F("qe_pt_cut3","", 35, 0.0, 120.0);
    TH1F *qePtCut5 = new TH1F("qe_pt_cut5","", 35, 0.0, 120.0);
    TH1F *qePtCheckPair = new TH1F("qe_pt_check_pair","", 35, 0.0, 120.0);
    
    TH1F *qmPtTotal = new TH1F("qm_full","Muons produced in Pythia HardQCD process collisions;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 35, 0.0, 120.0);
    TH1F *qmPtPart = new TH1F("qm_pt_part","", 35, 0.0, 120.0);
    TH1F *qmPtCut1 = new TH1F("qm_pt_cut1",";#hat{p}_{T} (GeV/c);#frac{N_{#mu}}{N_{#mu paired}}", 35, 0.0, 120.0);
    TH1F *qmPtCut3 = new TH1F("qm_pt_cut3","", 35, 0.0, 120.0);
    TH1F *qmPtCut5 = new TH1F("qm_pt_cut5","", 35, 0.0, 120.0);
    TH1F *qmPtCheckPair = new TH1F("qm_pt_check_pair","", 35, 0.0, 120.0);


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
        qeTuples[binCount]->Draw("ptHat>>qe_pt_part","ptL>1");
        qePtPart->Scale(1/(*it), "width");
        qePtCut1->Add(qePtPart);

        qePtPart->Reset();
        qeTuples[binCount]->Draw("ptHat>>qe_pt_part","ptL>3");
        qePtPart->Scale(1/(*it), "width");
        qePtCut3->Add(qePtPart);

        qePtPart->Reset();
        qeTuples[binCount]->Draw("ptHat>>qe_pt_part","ptL>5");
        qePtPart->Scale(1/(*it), "width");
        qePtCut5->Add(qePtPart);

        qePtPart->Reset();
        qeTuples[binCount]->Draw("ptHat>>qe_pt_part","parentFlag==1 && pairPartner==1");
        qePtPart->Scale(1/(*it), "width");
        qePtCheckPair->Add(qePtPart);

        //Muon Decays
      
        qmPtPart->Reset();
        qmTuples[binCount]->Draw("ptHat>>qm_pt_part");
        qmPtPart->Scale(1/(*it),"width");
        qmPtTotal->Add(qmPtPart);

        qmPtPart->Reset();
        qmTuples[binCount]->Draw("ptHat>>qm_pt_part","ptL>1");
        qmPtPart->Scale(1/(*it), "width");
        qmPtCut1->Add(qmPtPart);

        qmPtPart->Reset();
        qmTuples[binCount]->Draw("ptHat>>qm_pt_part","ptL>4");
        qmPtPart->Scale(1/(*it), "width");
        qmPtCut3->Add(qmPtPart);

        qmPtPart->Reset();
        qmTuples[binCount]->Draw("ptHat>>qm_pt_part","ptL>7");
        qmPtPart->Scale(1/(*it), "width");
        qmPtCut5->Add(qmPtPart);

        qmPtPart->Reset();
        qmTuples[binCount]->Draw("ptHat>>qm_pt_part","parentFlag==1 && pairPartner==1");
        qmPtPart->Scale(1/(*it), "width");
        qmPtCheckPair->Add(qmPtPart);

        binCount++;
    }

    ////Plotting
    // qes
    TFile *outf =  new TFile("sim_hists_individual_ratios_cuts.root", "RECREATE");
    TCanvas *canvasQE = new TCanvas("qe_sigma","qe_sigma");
    
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetLogy();
    pad1->SetBottomMargin(0);

    qePtTotal->GetYaxis()->SetTitleOffset(1);

    pad1->Draw();
    ///pad1->SetGridx();
    pad1->cd(); 

    qePtTotal->SetLineColor(kGray);
    qePtTotal->SetStats(0);
    //qePtTotal->GetYaxis()->SetLabelSize(0.);
    qePtTotal->Draw("SAME");

    qePtCut1->SetLineColor(kRed);
    qePtCut1->SetStats(0);
    qePtCut1->DrawCopy("SAME");

    qePtCut3->SetLineColor(kGreen);
    qePtCut3->SetStats(0);
    qePtCut3->DrawCopy("SAME");

    qePtCut5->SetLineColor(kBlue);
    qePtCut5->SetStats(0);
    qePtCut5->DrawCopy("SAME");

    qePtCheckPair->SetLineColor(kBlack);
    qePtCheckPair->SetStats(0);
    qePtCheckPair->DrawCopy("SAME");

    auto legendqe = new TLegend();
    legendqe->AddEntry(qePtTotal,"All e","l");
    legendqe->AddEntry(qePtCut1,"p_{T}^{e}>1 GeV","l");
    legendqe->AddEntry(qePtCut3,"p_{T}^{e}>3 GeV","l");
    legendqe->AddEntry(qePtCut5,"p_{T}^{e}>5 GeV","l");
    legendqe->AddEntry(qePtCheckPair,"parent cut and pair cut","l");
    
    legendqe->Draw("SAME");

    canvasQE->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);

    qePtCut1->GetXaxis()->SetTitleSize(13);
    qePtCut1->GetXaxis()->SetTitleOffset(1);
    qePtCut1->GetXaxis()->SetLabelSize(10);

    qePtCut1->GetYaxis()->SetTitleSize(13);
    qePtCut1->GetYaxis()->SetTitleOffset(1);
    qePtCut1->GetYaxis()->SetLabelSize(10);

    pad2->Draw("SAME");
    ///pad2->SetGridx();
    pad2->cd(); 

    qePtCut1->Divide(qePtCheckPair);
    qePtCut1->SetStats(0);
    qePtCut1->SetAxisRange(0,2,"Y");
    qePtCut1->DrawCopy();

    qePtCut3->Divide(qePtCheckPair);
    qePtCut3->SetStats(0);
    qePtCut3->DrawCopy("SAME");

    qePtCut5->Divide(qePtCheckPair);
    qePtCut5->SetStats(0);
    qePtCut5->DrawCopy("SAME");

    TLine *line = new TLine(0,1,120,1);
    line->SetLineColor(kBlack);
    line->SetLineStyle(kDashed);
    line->Draw();

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

    qmPtTotal->SetLineColor(kGray);
    qmPtTotal->SetStats(0);
    //qmPtTotal->GetYaxis()->SetLabelSize(0.);
    qmPtTotal->Draw("SAME");

    qmPtCut1->SetLineColor(kRed);
    qmPtCut1->SetStats(0);
    qmPtCut1->DrawCopy("SAME");

    qmPtCut3->SetLineColor(kGreen);
    qmPtCut3->SetStats(0);
    qmPtCut3->DrawCopy("SAME");

    qmPtCut5->SetLineColor(kBlue);
    qmPtCut5->SetStats(0);
    qmPtCut5->DrawCopy("SAME");

    qmPtCheckPair->SetLineColor(kBlack);
    qmPtCheckPair->SetStats(0);
    qmPtCheckPair->DrawCopy("SAME");

    auto legendqm = new TLegend();
    legendqm->AddEntry(qmPtTotal,"All #mu","l");
    legendqm->AddEntry(qmPtCut1,"p_{T}^{#mu}>1 GeV","l");
    legendqm->AddEntry(qmPtCut3,"p_{T}^{#mu}>3 GeV","l");
    legendqm->AddEntry(qmPtCut5,"p_{T}^{#mu}>5 GeV","l");
    legendqm->AddEntry(qmPtCheckPair,"parent cut and pair cut","l");
    
    legendqm->Draw("SAME");

    canvasQM->cd();
    TPad *pad2M = new TPad("pad2M", "pad2M", 0, 0.05, 1, 0.3);

    pad2M->SetTopMargin(0);
    pad2M->SetBottomMargin(0.2);

    qmPtCut1->GetXaxis()->SetTitleSize(13);
    qmPtCut1->GetXaxis()->SetTitleOffset(1);
    qmPtCut1->GetXaxis()->SetLabelSize(10);

    qmPtCut1->GetYaxis()->SetTitleSize(13);
    qmPtCut1->GetYaxis()->SetTitleOffset(1);
    qmPtCut1->GetYaxis()->SetLabelSize(10);

    pad2M->Draw("SAME");
    ///pad2M->SetGridx();
    pad2M->cd(); 

    qmPtCut1->Divide(qmPtCheckPair);
    qmPtCut1->SetStats(0);
    qmPtCut1->SetAxisRange(0,3,"Y");
    qmPtCut1->DrawCopy();

    qmPtCut3->Divide(qmPtCheckPair);
    qmPtCut3->SetStats(0);
    qmPtCut3->DrawCopy("SAME");

    qmPtCut5->Divide(qmPtCheckPair);
    qmPtCut5->SetStats(0);
    qmPtCut5->DrawCopy("SAME");

    TLine *lineM = new TLine(0,1,120,1);
    lineM->SetLineColor(kBlack);
    lineM->SetLineStyle(kDashed);
    lineM->Draw();


    canvasQM->Write();
    
    delete outf;

}