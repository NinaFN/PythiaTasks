#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

void read_sim() {
    TFile *f = TFile::Open("sim_tuples.root","READ");
    TFile *outf =  new TFile("sim_hists.root", "RECREATE");
    
    vector<double> *binLuminocity;
    vector<TNtuple*> muonTuples(7); //should be count of bins but 7 is easier
    vector<TNtuple*> elecTuples(7);

    f->GetObject("luminocities",binLuminocity);

    int binCount = 0;

    //Create Hists
    TH1F *muonPtTotal = new TH1F("muon_full","Produced Muon Cross-Section;p_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 35, 0.0, 70.0);
    TH1F *muonPtPart = new TH1F("muon_pt_part","", 35, 0.0, 70.0);
    TH1F *muonPtPrim = new TH1F("muon_pt_prim","", 35, 0.0, 70.0);
    

    for(std::vector<double>::iterator it = binLuminocity->begin(); it != binLuminocity->end(); ++it){
        
        //Fill Tuples
        muonTuples[binCount] = (TNtuple*)f->Get(Form("muon%d", binCount));
        elecTuples[binCount] = (TNtuple*)f->Get(Form("elec%d", binCount));

        ////Fill Histograms

        //Muons
        muonPtPart->Reset();
        muonTuples[binCount]->Draw("pt>>muon_pt_part");
        muonPtPart->Scale(1/(*it), "width");
        muonPtTotal->Add(muonPtPart);


      
        binCount++;
    }


    ////Plotting
    // Output Muons
    TCanvas *canvasMuon = new TCanvas("Muon_sigma","Muon_sigma");

    muonPtTotal->SetLineColor(kRed);

    muonPtTotal->Draw();

    muonPtPrim->SetLineColor(kOrange);
    muonPtPrim->SetStats(0);
    muonPtPrim->Draw("SAME");


    auto legendMuon = new TLegend();
    legendMuon->AddEntry(muonPtTotal,"Total #mu yield","l");
    legendMuon->AddEntry(muonPtPrim,"#mu yield from heavy quark primary decays","l");
    legendMuon->Draw("SAME");

    canvasMuon->Write();

    delete outf;

}