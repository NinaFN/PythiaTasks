#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

void read_sim_quarks_ratios() {
    TFile *f = TFile::Open("sim_tuples_quarks.root","READ");
    
    
    vector<double> *binLuminocity;
    vector<double> *binEvCounts;
    vector<TNtuple*> quarkTuples(7); //should be count of bins but 7 is easier

    f->GetObject("luminocities",binLuminocity);
    f->GetObject("eventCounts",binEvCounts);

    int binCount = 0;

    //Create Hists
    TH1F *quarkPtTotal = new TH1F("hf_full","Events that produce c#bar{c} or b#bar{b} pairs in Pythia HardQCD processes;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 35, 0.0, 120.0);
    TH1F *quarkPtPart = new TH1F("quark_pt_part","", 35, 0.0, 120.0);
    TH1F *quarkPt1 = new TH1F("quark_pt_1",";#hat{p}_{T} (GeV/c);Ratio", 35, 0.0, 120.0);
    TH1F *quarkPt2 = new TH1F("quark_pt_2","", 35, 0.0, 120.0);
    TH1F *quarkPt3 = new TH1F("quark_pt_3","", 35, 0.0, 120.0);
    TH1F *quarkPt4 = new TH1F("quark_pt_4","", 35, 0.0, 120.0);
    TH1F *quarkPt5 = new TH1F("quark_pt_5","", 35, 0.0, 120.0);
    TH1F *quarkPt6 = new TH1F("quark_pt_6","", 35, 0.0, 120.0);
    

    for(std::vector<double>::iterator it = binLuminocity->begin(); it != binLuminocity->end(); ++it){
        
        //Fill Tuples
        quarkTuples[binCount] = (TNtuple*)f->Get(Form("quark%d", binCount));

        ////Fill Histograms

        //Quarks
        /*quarkPtPart->Reset();
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
        quarkPt6->Add(quarkPtPart);*/
      
        quarkPtPart->Reset();
        quarkTuples[binCount]->Draw("ptHat>>quark_pt_part");
        quarkPtPart->Scale(1/(*it), "width");
        quarkPtTotal->Add(quarkPtPart);

        quarkPtPart->Reset();
        quarkTuples[binCount]->Draw("ptHat>>quark_pt_part","decayMap>=1 && decayMap!=2");
        quarkPtPart->Scale(1/(*it), "width");
        quarkPt1->Add(quarkPtPart);

        quarkPtPart->Reset();
        quarkTuples[binCount]->Draw("ptHat>>quark_pt_part","decayMap>=2");
        quarkPtPart->Scale(1/(*it), "width");
        quarkPt2->Add(quarkPtPart);

        quarkPtPart->Reset();
        quarkTuples[binCount]->Draw("ptHat>>quark_pt_part","decayMap>=3");
        quarkPtPart->Scale(1/(*it), "width");
        quarkPt3->Add(quarkPtPart);

        quarkPtPart->Reset();
        quarkTuples[binCount]->Draw("ptHat>>quark_pt_part","decayMap>=4");
        quarkPtPart->Scale(1/(*it), "width");
        quarkPt4->Add(quarkPtPart);

        quarkPtPart->Reset();
        quarkTuples[binCount]->Draw("ptHat>>quark_pt_part","decayMap>=5");
        quarkPtPart->Scale(1/(*it), "width");
        quarkPt5->Add(quarkPtPart);

        quarkPtPart->Reset();
        quarkTuples[binCount]->Draw("ptHat>>quark_pt_part","decayMap>=6");
        quarkPtPart->Scale(1/(*it), "width");
        quarkPt6->Add(quarkPtPart);

        binCount++;
    }

    ////Plotting
    // Quarks
    TFile *outf =  new TFile("sim_hists_quarks_ratios.root", "RECREATE");
    TCanvas *canvasQuark = new TCanvas("Quark_sigma","Quark_sigma");
    
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetLogy();
    pad1->SetBottomMargin(0);
    quarkPtTotal->GetYaxis()->SetTitleOffset(1);
    pad1->Draw();
    pad1->cd(); 

    quarkPtTotal->SetLineColor(kBlack);
    quarkPtTotal->SetStats(0);
    quarkPtTotal->Draw();

    quarkPt1->SetLineColor(kRed);
    quarkPt1->SetStats(0);
    quarkPt1->DrawCopy("SAME");

    quarkPt2->SetLineColor(kOrange);
    quarkPt2->SetStats(0);
    quarkPt2->DrawCopy("SAME");

    quarkPt3->SetLineColor(kGreen);
    quarkPt3->SetStats(0);
    quarkPt3->DrawCopy("SAME");

    quarkPt4->SetLineColor(kAzure+8);
    quarkPt4->SetStats(0);
    quarkPt4->DrawCopy("SAME");

    quarkPt5->SetLineColor(kViolet+2);
    quarkPt5->SetStats(0);
    quarkPt5->DrawCopy("SAME");

    quarkPt6->SetLineColor(kViolet);
    quarkPt6->SetStats(0);
    quarkPt6->DrawCopy("SAME");


    auto legendQuark = new TLegend();
    legendQuark->AddEntry(quarkPtTotal,"All","l");
    legendQuark->AddEntry(quarkPt1,"Single quark decay to #mu","l");
    legendQuark->AddEntry(quarkPt2,"Single quark decay to e","l");
    legendQuark->AddEntry(quarkPt3,"Pair decay to e and #mu","l");
    legendQuark->AddEntry(quarkPt4,"Back to back decay to e and #mu","l");
    legendQuark->AddEntry(quarkPt5,"Back to back decay to opposite signed e and #mu","l");
    legendQuark->AddEntry(quarkPt6,"Back to back decay to opposite signed e and #mu in correct #eta range","l");
    
    legendQuark->Draw("SAME");

    canvasQuark->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);

    quarkPt1->GetXaxis()->SetTitleSize(13);
    quarkPt1->GetXaxis()->SetTitleOffset(1);
    quarkPt1->GetXaxis()->SetLabelSize(10);

    quarkPt1->GetYaxis()->SetTitleSize(13);
    quarkPt1->GetYaxis()->SetTitleOffset(1);
    quarkPt1->GetYaxis()->SetLabelSize(10);

    pad2->Draw("SAME");
    pad2->cd(); 

    quarkPt1->Divide(quarkPtTotal);
    quarkPt1->SetStats(0);
    quarkPt1->SetAxisRange(0,1,"Y");
    quarkPt1->DrawCopy();

    quarkPt2->Divide(quarkPtTotal);
    quarkPt2->SetStats(0);
    quarkPt2->DrawCopy("SAME");

    quarkPt3->Divide(quarkPtTotal);
    quarkPt3->SetStats(0);
    quarkPt3->DrawCopy("SAME");

    quarkPt4->Divide(quarkPtTotal);
    quarkPt4->SetStats(0);
    quarkPt4->DrawCopy("SAME");

    quarkPt5->Divide(quarkPtTotal);
    quarkPt5->SetStats(0);
    quarkPt5->DrawCopy("SAME");

    quarkPt6->Divide(quarkPtTotal);
    quarkPt6->SetStats(0);
    quarkPt6->DrawCopy("SAME");

    canvasQuark->Write();

    delete outf;

}