#include <iostream>
#include "TCanvas.h"
#include "TError.h"
#include "TPad.h"
#include "TString.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TF1.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TLatex.h"
#include "TString.h"  
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCut.h"
#include "TNtuple.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TLine.h"

void plot_fake(int nstep_cent=0,int nstep_accept=0,int nstep_pt=0,int nstep_rmin=0,double  bin_pt_min=0.4,double bin_pt_max=1,int cent_min=0,int cent_max=200,int nevents=1000,bool is_final=0,TString algo="akPu3PF",bool dogenjet=false){
TH1D::SetDefaultSumw2();
TH2D::SetDefaultSumw2(true);
 // int nstep_cent=1;
 // int nstep_accept=1;
 // int nstep_pt=1;

// TFile * f = new TFile(Form("../workDir/tracking_ntuples/new/track_ntuple_Dijet100_HydjetDrum_v27_mergedV1_%devts_cent_%d_accept_%d_pt_%d_ptmin%d_ptmax%d.root",nevents,nstep_cent,nstep_accept,nstep_pt,bin_pt_min,bin_pt_max));
TFile * f = new TFile(Form("../ntuples_%s/track_ntuple_HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_0_%devts_cent_%d_accept_%d_pt_%d_rmin_%d_ptmin%d_ptmax%d_cent%d_%d_%s_dogenjet%d.root",algo.Data(),nevents,nstep_cent,nstep_accept,nstep_pt,nstep_rmin,(int)bin_pt_min,(int)bin_pt_max,cent_min,cent_max,algo.Data(),dogenjet));
TString s=Form("../ntuples_%s/track_ntuple_HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_0_%devts_cent_%d_accept_%d_pt_%d_rmin_%d_ptmin%d_ptmax%d_cent%d_%d_%s_dogenjet%d.root",algo.Data(),nevents,nstep_cent,nstep_accept,nstep_pt,nstep_rmin,(int)bin_pt_min,(int)bin_pt_max,cent_min,cent_max,algo.Data(),dogenjet);
cout<<s.Data()<<endl;

TTree * nt = (TTree*)f->Get("nt_track");
//##########################################Initial fakeiciencies without Correction###########################################################################

/////////eta dependent fake////////////////
TProfile * p_eta= new TProfile("p_eta",";#eta;fake",20,-2.4,2.4);
TH1D * eta_gen = new TH1D("eta_gen",";#eta;N_{part}",20,-2.4,2.4);
TH1D * eta_reco = new TH1D("eta_reco",";#eta;N_{part}",20,-2.4,2.4);
nt->Draw("eta>>eta_gen","");
nt->Draw("eta>>eta_reco","trkfake");
nt->Draw("(trkfake):eta>>p_eta","trackselect && pt>1");
eta_gen->SetMarkerStyle(25);

TLegend *t5=new TLegend(0.4,0.85,0.9,0.95); 
t5->SetFillColor(0);
t5->SetBorderSize(0);
t5->SetFillStyle(0);
t5->SetTextFont(43);
t5->SetTextSize(22);
t5->AddEntry(eta_gen,"gen","p");
t5->AddEntry(eta_reco,"matched gen","p");

TCanvas * c1 = new TCanvas("c1","",600,600);
p_eta->Draw();
c1->SaveAs("../plots/fake_eta.png");

TCanvas *c1_2 = new TCanvas("c1_2","",600,600);
eta_gen->Draw();
eta_gen->SetMinimum(eta_reco->GetMinimum());
eta_gen->SetMaximum(1.2*eta_gen->GetMaximum());
eta_reco->Draw("same");
t5->Draw("same");
c1_2->SaveAs("../plots/eta_dists.png");


///////////acceptance dependent fake/////////
TProfile2D * p_eta_phi = new TProfile2D(Form("p_eta_phi"),";#phi;#eta;",20,-TMath::Pi(),TMath::Pi(),20,-2.4,2.4);
nt->Draw(Form("(trkfake):eta:phi>>p_eta_phi"),Form("(trackselect&& abs(eta)<2.4&& pt>%.1f && pt<%.1f)",bin_pt_min,bin_pt_max));

TCanvas * c2 = new TCanvas("c2","",600,600);
c2->SetRightMargin(0.2);
p_eta_phi->Draw("colz");
c2->SaveAs(Form("../plots/p_eta_phi.png"));
cout<<4<<endl;

//////////pt dependent fake////////////

const int ny=20;
double x[ny+1];
double inix=log(bin_pt_min)/log(10);
double delta=(log(bin_pt_max)-log(bin_pt_min))/(20*log(10));
int maxbin=ny;
for(int ix=0; ix<ny+1;ix++){
 x[ix]=pow(10,inix);
 if(x[ix]>100){
  x[ix]=bin_pt_max;
  maxbin=ix;
  break;
 }
 inix+=delta;
} 

TProfile * p_pt= new TProfile("p_pt",";p_{T}(GeV/c);fake",maxbin,x);
TH1D * pt_gen = new TH1D("pt_gen",";p_{T}(GeV/c);N_{part}",maxbin,x);
TH1D * pt_reco = new TH1D("pt_reco",";p_{T}(GeV/c);N_{part}",maxbin,x);
nt->Draw("pt>>pt_gen","");
nt->Draw("pt>>pt_reco","trkfake");
nt->Draw("(trkfake):pt>>p_pt","trackselect ");

TCanvas * c3 = new TCanvas("c3","",600,600);
c3->SetLogx();
p_pt->Draw();
c3->SaveAs("../plots/fake_pt.png");

TCanvas *c3_2 = new TCanvas("c3_2","",600,600);
c3_2->SetLogy();
pt_gen->SetMarkerStyle(25);
pt_gen->SetMinimum(pt_reco->GetMinimum());
pt_gen->SetMaximum(1.2*pt_gen->GetMaximum());
pt_gen->Draw();
pt_reco->Draw("same");
t5->Draw("same");
c3_2->SaveAs("../plots/pt_dists.png");


////////centrality dependent fake//////////

TProfile * p_cent= new TProfile("p_cent",";centrality bin;fake",(int)((cent_max-cent_min)/4),cent_min,cent_max);
TH1D * cent_gen = new TH1D("cent_gen",";centrality bin;N_{part}",(int)((cent_max-cent_min)/4),cent_min,cent_max);
TH1D * cent_reco = new TH1D("cent_reco",";centrality bin;N_{part}",(int)((cent_max-cent_min)/4),cent_min,cent_max);
nt->Draw("cent>>cent_gen","");
nt->Draw("cent>>cent_reco","trkfake");
nt->Draw("(trkfake):cent>>p_cent",Form("(trackselect && abs(eta)<2.4&& pt>%.1f && pt<%.1f)",bin_pt_min,bin_pt_max));

TH1D * fake_cent=(TH1D*)cent_reco->Clone("fake_cent");
fake_cent->Divide(cent_gen);

TCanvas * c4 = new TCanvas("c4","",600,600);
p_cent->Draw();
c4->SaveAs("../plots/fake_cent.png");


TCanvas *c4_2 = new TCanvas("c4_2","",600,600);
cent_gen->SetMarkerStyle(25);
cent_gen->SetMinimum(cent_gen->GetMinimum());
cent_gen->SetMaximum(1.2*cent_gen->GetMaximum());
cent_gen->Draw();
cent_reco->Draw("same");
t5->Draw("same");
c4_2->SaveAs("../plots/cent_dists.png");

//###############################################fakeiciencies after Correction###########################################################################

///////////////cent dependent///////////////////
TProfile * p_cent_fake= new TProfile("p_cent_fake",";centrality bin;fake",40,0,40);
// nt->Draw("(trkfake-fake):cent>>p_cent_fake",Form("(trackselect && abs(eta)<2.4&& pt>%.1f && pt<%.1f)",bin_pt_min,bin_pt_max));
nt->Draw("(trkfake-fake):cent>>p_cent_fake","(trackselect && abs(eta)<2.4&& trkfake<2)");


TCanvas * c5 = new TCanvas("c5","",600,600);
p_cent_fake->SetMaximum(1.1);
p_cent_fake->SetMinimum(-0.5);
p_cent_fake->Draw();
c5->SaveAs(Form("../plots/fake_cent_after_fake_nstep_cent%d_accept%d_pt%d.png",nstep_cent,nstep_accept,nstep_pt));


////////////pt dependent////////////////////////
TProfile * p_pt_fake= new TProfile("p_pt_fake",";p_{T}(GeV/c);fake",maxbin,x);
// nt->Draw("(trkfake-fake):pt>>p_pt_fake",Form("(trackselect && abs(eta)<2.4&& pt>%.1f && pt<%.1f)",bin_pt_min,bin_pt_max));
nt->Draw("(trkfake-fake):pt>>p_pt_fake","(trackselect && abs(eta)<2.4&& trkfake<2)");

TCanvas * c6 = new TCanvas("c6","",600,600);
c6->SetLogx();
p_pt_fake->SetMaximum(1.1);
p_pt_fake->SetMinimum(-0.5);
p_pt_fake->Draw();
c6->SaveAs(Form("../plots/fake_pt_after_fake_nstep_cent%d_accept%d_pt%d_rmin%d_pt%d_%d.png",nstep_cent,nstep_accept,nstep_pt,nstep_rmin,(int)bin_pt_min,(int)bin_pt_max));
c6->SaveAs(Form("../plots/fake_pt_after_fake_nstep_cent%d_accept%d_pt%d_rmin%d_pt%d_%d.pdf",nstep_cent,nstep_accept,nstep_pt,nstep_rmin,(int)bin_pt_min,(int)bin_pt_max));

//////////acceptance dependent///////////////////

TProfile2D * p_eta_phi_fake = new TProfile2D("p_eta_phi_fake",";#phi;#eta;",20,-TMath::Pi(),TMath::Pi(),20,-2.4,2.4);
// nt->Draw("(trkfake-fake):eta:phi>>p_eta_phi_fake",Form("(trackselect&& abs(eta)<2.4&& pt>%.1f && pt<%.1f)",bin_pt_min,bin_pt_max));
nt->Draw("(trkfake-fake):eta:phi>>p_eta_phi_fake","(trackselect&& abs(eta)<2.4 && trkfake<2)");

TCanvas * c7 = new TCanvas("c7","",600,600);
c7->SetRightMargin(0.2);
p_eta_phi_fake->Draw("colz");
c7->SaveAs(Form("../plots/fake_accept_after_fake_nstep_cent%d_accept%d_pt%d_rmin%d_pt%d_%d.png",nstep_cent,nstep_accept,nstep_pt,nstep_rmin,(int)bin_pt_min,(int)bin_pt_max));
c7->SaveAs(Form("../plots/fake_accept_after_fake_nstep_cent%d_accept%d_pt%d_rmin%d_pt%d_%d.pdf",nstep_cent,nstep_accept,nstep_pt,nstep_rmin,(int)bin_pt_min,(int)bin_pt_max));

///////////eta dependent/////////////////////////
TProfile * p_eta_fake = new TProfile("p_eta_fake",";#eta;fake",20,-2.4,2.4);

nt->Draw("(trkfake-fake):eta>>p_eta_fake","(trackselect&& abs(eta)<2.4 && trkfake<2)");

TCanvas * c8 = new TCanvas("c8","",600,600);
c8->SetRightMargin(0.2);
p_eta_fake->SetMaximum(1.1);
p_eta_fake->SetMinimum(-0.5);
p_eta_fake->Draw();
c8->SaveAs(Form("../plots/fake_eta_after_fake_nstep_cent%d_accept%d_pt%d_rmin%d_pt%d_%d.png",nstep_cent,nstep_accept,nstep_pt,nstep_rmin,(int)bin_pt_min,(int)bin_pt_max));
c8->SaveAs(Form("../plots/fake_eta_after_fake_nstep_cent%d_accept%d_pt%d_rmin%d_pt%d_%d.pdf",nstep_cent,nstep_accept,nstep_pt,nstep_rmin,(int)bin_pt_min,(int)bin_pt_max));


////////////phi dependent/////////////////////////
TProfile * p_phi_fake = new TProfile("p_phi_fake",";#phi;fake",20,-TMath::Pi(),TMath::Pi());
nt->Draw("(trkfake-fake):phi>>p_phi_fake","(trackselect&& abs(eta)<2.4 && trkfake<2)");

TCanvas * c9 = new TCanvas("c9","",600,600);
c9->SetRightMargin(0.2);

p_phi_fake->SetMaximum(1.1);
p_phi_fake->SetMinimum(-0.5);
p_phi_fake->Draw();
c9->SaveAs(Form("../plots/fake_phi_after_fake_nstep_cent%d_accept%d_pt%d_rmin%d_pt%d_%d.png",nstep_cent,nstep_accept,nstep_pt,nstep_rmin,(int)bin_pt_min,(int)bin_pt_max));
c9->SaveAs(Form("../plots/fake_phi_after_fake_nstep_cent%d_accept%d_pt%d_rmin%d_pt%d_%d.pdf",nstep_cent,nstep_accept,nstep_pt,nstep_rmin,(int)bin_pt_min,(int)bin_pt_max));

////////////rmin dependent/////////////////////////

const int n_rmin_bins=45;
double rmin_bins[n_rmin_bins+1];
for(int i=0; i<=n_rmin_bins; i++){
        if(i<=20) rmin_bins[i]=0.05*i;
        else if(i<=30) rmin_bins[i]=1+0.1*(i-20);
        else if(i<=45) rmin_bins[i]=2+0.2*(i-30);
}

TProfile * p_rmin_fake = new TProfile("p_rmin_fake",";r_{min};fake",n_rmin_bins,rmin_bins);
nt->Draw("(trkfake-fake):rmin>>p_rmin_fake","(trackselect&& abs(eta)<2.4 && trkfake<2)");

TCanvas * c10 = new TCanvas("c10","",600,600);
c10->SetRightMargin(0.2);

p_rmin_fake->SetMaximum(1.1);
p_rmin_fake->SetMinimum(-0.5);
p_rmin_fake->Draw();
c10->SaveAs(Form("../plots/fake_rmin_after_fake_nstep_cent%d_accept%d_pt%d_rmin%d_pt%d_%d.png",nstep_cent,nstep_accept,nstep_pt,nstep_rmin,(int)bin_pt_min,(int)bin_pt_max));
c10->SaveAs(Form("../plots/fake_rmin_after_fake_nstep_cent%d_accept%d_pt%d_rmin%d_pt%d_%d.pdf",nstep_cent,nstep_accept,nstep_pt,nstep_rmin,(int)bin_pt_min,(int)bin_pt_max));

/////output file to be used in the ntupler in the next step/////
TFile *outf = new TFile(Form("../stepwise_hists_%s/fake_nstep_cent%d_accept%d_pt%d_rmin%d_pt%d_%d_cent%d_%d.root",algo.Data(),nstep_cent,nstep_accept,nstep_pt,nstep_rmin,(int)bin_pt_min,(int)bin_pt_max,cent_min,cent_max),"recreate");
p_cent->Write();
p_pt->Write();
p_eta_phi->Write();
p_cent_fake->Write();
p_pt_fake->Write();
p_eta_phi_fake->Write();
p_rmin_fake->Write();

outf->Close();

////overall fake histograms///////////////////////////////
TProfile * p_fake_cent = new TProfile("p_fake_cent",";centrality bin;fake",(int)((cent_max-cent_min)/4),cent_min,cent_max);
nt->Draw("fake_cent:cent>>p_fake_cent","(trackselect&& abs(eta)<2.4&& trkfake<2)","prof");

TProfile * p_fake_pt = new TProfile("p_fake_pt",";p_{T} bin;fake",maxbin,x);
nt->Draw("fake_pt:pt>>p_fake_pt","(trackselect&& abs(eta)<2.4&& trkfake<2)","prof");

TProfile2D * p_fake_acceptance = new TProfile2D("p_fake_acceptance",";#phi;#eta;fake",20,-TMath::Pi(),TMath::Pi(),20,-2.4,2.4);
nt->Draw("fake_accept:eta:phi>>p_fake_acceptance","(trackselect&& abs(eta)<2.4&& trkfake<2)","prof");

TProfile * p_fake_rmin = new TProfile("p_fake_rmin",";#phi;#eta;fake",n_rmin_bins,rmin_bins);
nt->Draw("fake_rmin:rmin>>p_fake_rmin","(trackselect&& abs(eta)<2.4&& trkfake<2)","prof");


TFile *f_fake;
 if(is_final){
 f_fake= new TFile(Form("../final_hists_%s_20140312/fake_pt%d_%d_cent%d_%d_step_cent%daccept%dpt%drmin%d_%s_dogenjet%d.root",algo.Data(),(int)bin_pt_min,(int)bin_pt_max,cent_min,cent_max,nstep_cent,nstep_accept,nstep_pt,nstep_rmin,algo.Data(),dogenjet),"recreate");

 p_fake_cent->Write();
 p_fake_pt->Write();
 p_fake_acceptance->Write();
 p_fake_rmin->Write();

 f_fake->Close();
}

}
