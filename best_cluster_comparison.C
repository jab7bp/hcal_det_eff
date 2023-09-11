#include <iostream>
#include <fstream>
#include <numeric>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <string>
#include <chrono>
#include <TF1.h>
#include "TStopwatch.h"

// using namespace std::chrono;
#include "/w/halla-scshelf2102/sbs/jboyd/include/include_files.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/GEM_lookups.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/beam_variables.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/calc_functions.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/utility_functions.h"

TString tf_filename_base = "rootfiles/det_eff_histos_SBS4_LH2_mag0_bgSub_W2BinFactor200_dxBinFactor200_dySigMult2_5_w2SigMult4_0_scoreTest";

TFile *tf_e = new TFile(Form("%s_1.root", tf_filename_base.Data()), "READ");
TFile *tf_dx = new TFile(Form("%s_2.root", tf_filename_base.Data()), "READ");
TFile *tf_ADC = new TFile(Form("%s_3.root", tf_filename_base.Data()), "READ");

TH1D *h_dx_wcut_e, *h_dx_wcut_dx, *h_dx_wcut_ADC;

void best_cluster_comparison(){

	h_dx_wcut_e = static_cast<TH1D*>(tf_e->Get("h_dx_wcut"));
	h_dx_wcut_dx = static_cast<TH1D*>(tf_dx->Get("h_dx_wcut"));
	h_dx_wcut_ADC = static_cast<TH1D*>(tf_ADC->Get("h_dx_wcut"));

	TCanvas *c = new TCanvas("c", "c", 600, 500);

	h_dx_wcut_dx->SetLineColor(6);
	h_dx_wcut_dx->GetXaxis()->SetRangeUser(-0.2, 0.4);
	h_dx_wcut_dx->Draw();

	h_dx_wcut_e->SetLineColor(kSpring);
	h_dx_wcut_e->SetFillColorAlpha(kSpring, 0.25);
	h_dx_wcut_e->SetFillStyle(3005);
	h_dx_wcut_e->Draw("same");

	h_dx_wcut_ADC->SetLineColor(kBlue);
	h_dx_wcut_ADC->Draw("same");

	TLegend *tl_dx_best_cluster = new TLegend(0.60, 0.6, 0.75, 0.75, "Default best cluster");
	tl_dx_best_cluster->AddEntry(h_dx_wcut_dx, "Min dx");
	tl_dx_best_cluster->AddEntry(h_dx_wcut_e, "Max hcal_e");
	tl_dx_best_cluster->AddEntry(h_dx_wcut_ADC, "Min ADC time");
	tl_dx_best_cluster->Draw("same");

	c->Modified();

}