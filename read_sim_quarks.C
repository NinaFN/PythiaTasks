#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

void read_sim_quarks() {
    TFile *f = TFile::Open("sim_tuples_quarks.root","READ");
    
    
    vector<double> *binLuminocity;
    vector<TNtuple*> quarkTuples(7); //should be count of bins but 7 is easier

    f->GetObject("luminocities",binLuminocity);

    int binCount = 0;

    //Create Hists
    TH1F *quarkPtTotal = new TH1F("hf_full","Heavy Quark from Hardest Process Cross Section;p_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 35, 0.0, 70.0);
    TH1F *quarkPtPart = new TH1F("quark_pt_part","", 35, 0.0, 70.0);
    TH1F *quarkPt1 = new TH1F("quark_pt_1","", 35, 0.0, 70.0);
    TH1F *quarkPt2 = new TH1F("quark_pt_2","", 35, 0.0, 70.0);
    TH1F *quarkPt3 = new TH1F("quark_pt_3","", 35, 0.0, 70.0);
    TH1F *quarkPt4 = new TH1F("quark_pt_4","", 35, 0.0, 70.0);
    TH1F *quarkPt5 = new TH1F("quark_pt_5","", 35, 0.0, 70.0);
    TH1F *quarkPt6 = new TH1F("quark_pt_6","", 35, 0.0, 70.0);
    

    for(std::vector<double>::iterator it = binLuminocity->begin(); it != binLuminocity->end(); ++it){
        
        //Fill Tuples
        quarkTuples[binCount] = (TNtuple*)f->Get(Form("quark%d", binCount));

        ////Fill Histograms

        //Quarks
        quarkPtPart->Reset();
        quarkTuples[binCount]->Draw("pt1>>quark_pt_part");
        quarkPtPart->Scale(1/(*it), "width");
        quarkPtTotal->Add(quarkPtPart);

        quarkPtPart->Reset();
        quarkTuples[binCount]->Draw("pt1>>quark_pt_part","decayMap>=1 && decayMap!=2");
        quarkPtPart->Scale(1/(*it), "width");
        quarkPt1->Add(quarkPtPart);

        quarkPtPart->Reset();
        quarkTuples[binCount]->Draw("pt1>>quark_pt_part","decayMap>=2");
        quarkPtPart->Scale(1/(*it), "width");
        quarkPt2->Add(quarkPtPart);

        quarkPtPart->Reset();
        quarkTuples[binCount]->Draw("pt1>>quark_pt_part","decayMap>=3");
        quarkPtPart->Scale(1/(*it), "width");
        quarkPt3->Add(quarkPtPart);

        quarkPtPart->Reset();
        quarkTuples[binCount]->Draw("pt1>>quark_pt_part","decayMap>=4");
        quarkPtPart->Scale(1/(*it), "width");
        quarkPt4->Add(quarkPtPart);

        quarkPtPart->Reset();
        quarkTuples[binCount]->Draw("pt1>>quark_pt_part","decayMap>=5");
        quarkPtPart->Scale(1/(*it), "width");
        quarkPt5->Add(quarkPtPart);

        quarkPtPart->Reset();
        quarkTuples[binCount]->Draw("pt1>>quark_pt_part","decayMap>=6");
        quarkPtPart->Scale(1/(*it), "width");
        quarkPt6->Add(quarkPtPart);
      
        binCount++;
    }


    ////Plotting
    // Quarks
    TFile *outf =  new TFile("sim_hists_quarks.root", "RECREATE");
    TCanvas *canvasQuark = new TCanvas("Quark_sigma","Quark_sigma");

    quarkPtTotal->SetLineColor(kBlack);
    quarkPtTotal->Draw();

    quarkPt1->SetLineColor(kRed);
    quarkPt1->SetStats(0);
    quarkPt1->Draw("SAME");

    quarkPt2->SetLineColor(kOrange);
    quarkPt2->SetStats(0);
    quarkPt2->Draw("SAME");

    quarkPt3->SetLineColor(kGreen);
    quarkPt3->SetStats(0);
    quarkPt3->Draw("SAME");

    quarkPt4->SetLineColor(kAzure+8);
    quarkPt4->SetStats(0);
    quarkPt4->Draw("SAME");

    quarkPt5->SetLineColor(kViolet);
    quarkPt5->SetStats(0);
    quarkPt5->Draw("SAME");

    quarkPt6->SetLineColor(kViolet+2);
    quarkPt6->SetStats(0);
    quarkPt6->Draw("SAME");


    auto legendQuark = new TLegend();
    legendQuark->AddEntry(quarkPtTotal,"Total events with heavy quark pairs","l");
    legendQuark->AddEntry(quarkPt1,"Events with heavy quarks decaying to #mu","l");
    legendQuark->AddEntry(quarkPt2,"Events with heavy quarks decaying to e","l");
    legendQuark->AddEntry(quarkPt3,"Events with heavy quarks decaying to e and #mu","l");
    legendQuark->AddEntry(quarkPt4,"Events with heavy quark pairs decaying to e and #mu","l");
    legendQuark->AddEntry(quarkPt5,"Events with heavy quark pairs decaying to opposite signed e and #mu","l");
    legendQuark->AddEntry(quarkPt6,"Events with heavy quark pairs decaying to opposite signed e and #mu in correct #eta range","l");
    
    legendQuark->Draw("SAME");

    canvasQuark->Write();

    delete outf;

}