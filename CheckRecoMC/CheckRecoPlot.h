#include "TH1D.h"
#include "TH2D.h"

#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "Rtypes.h"
#include "TFile.h"

class CheckRecoPlot {
    public:
        TH1D *h1_n_ch[2]; // charged tracks multi.
        TH1D *h1_dtan_theta, *h1_dphi; // reco. - truth
        TH2D *h2_nseg_pz, *h2_dtan_theta_pz, *h2_dphi_pz;
        TCanvas *cvs;
        TLegend *lg;
        const TString str_cat[2] = {"reco", "truth"};
        const int str_color[2] = {kRed, kBlue};

        CheckRecoPlot();
        
        virtual void PlotHist(TString StoreAddress);
        virtual void StoreHist2ROOT(TString StoreAddress);

};


CheckRecoPlot::CheckRecoPlot(){

    for(int i=0; i<2; i++){
        h1_n_ch[i] = new TH1D("h1_n_ch_" + str_cat[i], "Charged Tracks Multi.", 20, 0, 20);
        h1_n_ch[i]->GetXaxis()->SetTitle("n_{ch}");
        h1_n_ch[i]->GetYaxis()->SetTitle("counts");
        h1_n_ch[i]->SetLineColor(str_color[i]);
    }

    h1_dtan_theta = new TH1D("h1_dtan_theta", "tan#theta_{reco} - tan#theta_{true}", 100, -0.05, 0.05);
    h1_dtan_theta->GetXaxis()->SetTitle("tan#theta_{reco} - tan#theta_{true}");
    h1_dtan_theta->GetYaxis()->SetTitle("counts");

    h1_dphi = new TH1D("h1_dphi", "#phi_{reco} - #phi_{true}", 400, -0.2, 0.2);
    h1_dphi->GetXaxis()->SetTitle("(#phi_{reco} - #phi_{true}) [rad]");
    h1_dphi->GetYaxis()->SetTitle("counts");

    h2_nseg_pz = new TH2D("h2_nseg_pz", "# seg versus P_{z}", 100, 0, 100, 80, 0, 80);
    h2_nseg_pz->GetXaxis()->SetTitle("P_{z} [GeV]");
    h2_nseg_pz->GetYaxis()->SetTitle("# seg.");

    h2_dphi_pz = new TH2D("h2_dphi_pz", "(#phi_{reco} - #phi_{true}) versus P_{z}", 100, 0, 100, 400, -0.2, 0.2);
    h2_dphi_pz->GetXaxis()->SetTitle("P_{z} [GeV]");
    h2_dphi_pz->GetYaxis()->SetTitle("(#phi_{reco} - #phi_{true}) [rad]");

    h2_dtan_theta_pz = new TH2D("h2_dtan_theta_pz", "", 100, 0, 100, 100, -0.05, 0.05);
    h2_dtan_theta_pz->GetXaxis()->SetTitle("P_{z} [GeV]");
    h2_dtan_theta_pz->GetYaxis()->SetTitle("tan#theta_{reco} - tan#theta_{true}");

};

void CheckRecoPlot::PlotHist(TString StoreAddress){

    cvs = new TCanvas();
    lg = new TLegend(.1, .7, .3, .9);

    cvs->cd();
    h1_n_ch[1]->Draw("HIST");
    h1_n_ch[0]->Draw("HIST SAME");
    lg->AddEntry(h1_n_ch[0], "reco.");
    lg->AddEntry(h1_n_ch[1], "truth");
    lg->SetFillStyle(0);
    lg->Draw("same");
    cvs->SetGrid();       
    cvs->SaveAs(StoreAddress+"/h1_n_ch.pdf");

    cvs->Clear();
    cvs->cd();
    h1_dtan_theta->Draw("HIST");
    cvs->SetGrid();
    cvs->SaveAs(StoreAddress+"/h1_dtan_theta.pdf");

    cvs->Clear();
    cvs->cd();
    h1_dphi->Draw("HIST");
    cvs->SetGrid();
    cvs->SaveAs(StoreAddress+"/h1_dphi.pdf");
/*
    cvs->Clear();
    cvs->cd();
    h2_nseg_pz->Draw("colz");
    cvs->SetGrid();
    cvs->SaveAs(StoreAddress+"/h2_nseg_pz.pdf");
*/
}

void CheckRecoPlot::StoreHist2ROOT(TString StoreAddress){

    TFile *f_output = new TFile(StoreAddress + "/CheckRecoPlot.root", "RECREATE");
    f_output->cd();

    h1_n_ch[0]->Write();
    h1_n_ch[1]->Write();
    h1_dtan_theta->Write();
    h1_dphi->Write();
    h2_nseg_pz->Write();
    h2_dphi_pz->Write();
    h2_dtan_theta_pz->Write();

    f_output->Save();
    f_output->Close();
    delete f_output;

}

