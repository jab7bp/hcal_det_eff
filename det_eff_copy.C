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

using namespace std::chrono;
#include "/w/halla-scshelf2102/sbs/jboyd/include/include_files.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/GEM_lookups.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/beam_variables.h"

int polN = 6;
bool use_parsed = false;

bool use_bbcal_cuts = false;
bool fiducial_cut = false;
bool correct_beam_energy = false;
bool apply_fcut = false;

bool fit_only = true;
bool calc_W_only = false;
bool fit_with_reject = true;

TF1 *tf1_W2_interpolate;

TGraph *gr_gaus1, *gr_polyBG;

Double_t fit_expo(Double_t *x, Double_t *par){
	double expo = 0.0;
	expo = par[0] - (par[0] - par[1])*exp(-par[2]*x[0]);
	return expo;
}

Double_t fitPol4(Double_t *x, Double_t *par){
	return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] + par[4]*x[0]*x[0]*x[0]*x[0];
}

Double_t fitPol5(Double_t *x, Double_t *par){
	return par[0] + par[1]*x[0] + par[2]*pow(x[0], 2) + par[3]*pow(x[0], 3) + par[4]*pow(x[0], 4) + par[5]*pow(x[0], 5);
}

Double_t fitPol6(Double_t *x, Double_t *par){
	return par[0] + par[1]*x[0] + par[2]*pow(x[0], 2) + par[3]*pow(x[0], 3) + par[4]*pow(x[0], 4) + par[5]*pow(x[0], 5) + par[6]*pow(x[0], 6);
}

Double_t fitPol7(Double_t *x, Double_t *par){
	return par[0] + par[1]*x[0] + par[2]*pow(x[0], 2) + par[3]*pow(x[0], 3) + par[4]*pow(x[0], 4) + par[5]*pow(x[0], 5) + par[6]*pow(x[0], 6) + par[7]*pow(x[0], 7);
}

Double_t fitPol8(Double_t *x, Double_t *par){
	return par[0] + par[1]*x[0] + par[2]*pow(x[0], 2) + par[3]*pow(x[0], 3) + par[4]*pow(x[0], 4) + par[5]*pow(x[0], 5) + par[6]*pow(x[0], 6) + par[7]*pow(x[0], 7) + par[8]*pow(x[0], 8);
}

Double_t fitPol10(Double_t *x, Double_t *par){
	return par[0] + par[1]*x[0] + par[2]*pow(x[0], 2) + par[3]*pow(x[0], 3) + par[4]*pow(x[0], 4) + par[5]*pow(x[0], 5) + par[6]*pow(x[0], 6) + par[7]*pow(x[0], 7) + par[8]*pow(x[0], 8) + par[9]*pow(x[0], 9) + par[10]*pow(x[0], 10);
}


Double_t fitPol4_with_reject(Double_t *x, Double_t *par){
	
	double reject_min = 0.50217680;
	double reject_max = 1.09;

	if( x[0] > (reject_min) && x[0] < (reject_max) ){
		TF1::RejectPoint();
		return 0;
	}
	return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] + par[4]*x[0]*x[0]*x[0]*x[0];
}

Double_t fitPol6_with_reject(Double_t *x, Double_t *par){

	double reject_min = 0.50217680;
	double reject_max = 1.09;

	if( x[0] > (reject_min) && x[0] < (reject_max) ){
		TF1::RejectPoint();
		return 0;
	}
	return par[0] + par[1]*x[0] + par[2]*pow(x[0], 2) + par[3]*pow(x[0], 3) + par[4]*pow(x[0], 4) + par[5]*pow(x[0], 5) + par[6]*pow(x[0], 6);
}

Double_t fitPol8_with_reject(Double_t *x, Double_t *par){

	double reject_min = 0.50217680;
	double reject_max = 1.09;

	if( x[0] > (reject_min) && x[0] < (reject_max) ){
		TF1::RejectPoint();
		return 0;
	}
	return par[0] + par[1]*x[0] + par[2]*pow(x[0], 2) + par[3]*pow(x[0], 3) + par[4]*pow(x[0], 4) + par[5]*pow(x[0], 5) + par[6]*pow(x[0], 6) + par[7]*pow(x[0], 7) + par[8]*pow(x[0], 8);
}

Double_t W2_BG_fit_with_reject( Double_t *x, Double_t *par ){
	if( x[0] > 0.55 && x[0] < 1.2 ){
		TF1::RejectPoint();
		return 0;
	}
	return exp(par[0] +par[1]*x[0]);
}

Double_t W2_BG_fit_NO_reject( Double_t *x, Double_t *par ){
	return exp(par[0] +par[1]*x[0]);
}

Double_t fit_gaus(Double_t * x, Double_t *par){

	Double_t g = 0.0;

	g = par[0]*exp((-0.5)*pow(((x[0] -  par[1])/par[2]),2));
	
	return g;
}

Double_t fullFitFunction_reject( Double_t *x, Double_t *par ){
	// return fit_gaus(x, par) + W2_BG_fit_NO_reject(x, &par[3]);
	return fit_gaus(x, par) + fitPol6_with_reject(x, &par[3]);
	// return fit_gaus(x, par) + fitPol8_with_reject(x, &par[3]);
}

Double_t fullFitFunction( Double_t *x, Double_t *par ){
	// return fit_gaus(x, par) + W2_BG_fit_NO_reject(x, &par[3]);
	// return fit_gaus(x, par) + fitPol4(x, &par[3]);
	// return fit_gaus(x, par) + fitPol5(x, &par[3]);
	// return fit_gaus(x, par) + fitPol6(x, &par[3]);
	// return fit_gaus(x, par) + fitPol7(x, &par[3]);
	return fit_gaus(x, par) + fitPol8(x, &par[3]);
	// return fit_gaus(x, par) + fitPol10(x, &par[3]);
}

Double_t plot_gaus(Double_t *par){
	Double_t fit_x[2400], gaus1_y[2400];
	Int_t fit_n = 800;
	int fit_cnt = 0;

	for(double i = -2.5; i < 2.5; i+=.01){
		fit_x[fit_cnt] = i;
		gaus1_y[fit_cnt] = par[0]*exp((-0.5)*pow(((i -  par[1])/par[2]),2));
		fit_cnt++;
	}

	TGraph *gr_gaus1 = new TGraph(fit_n, fit_x, gaus1_y);
	gr_gaus1->SetLineColor(1);
	gr_gaus1->SetLineStyle(7);

	gr_gaus1->Draw("same");
	
	return 1;
}

Double_t plot_sig_and_BG(Double_t *par){
	Double_t fit_x[180], gaus1_y[180], polyBG_y[180];
	Int_t fit_n = 180;
	int fit_cnt = 0;

	for(double i = 0.4; i < 1.4; i+=.01){
		fit_x[fit_cnt] = i;
		gaus1_y[fit_cnt] = par[0]*exp((-0.5)*pow(((i -  par[1])/par[2]),2));
		polyBG_y[fit_cnt] = par[11]*pow(i, 8) + par[10]*pow(i, 7) + par[9]*pow(i, 6) + par[8]*pow(i, 5) + par[7]*pow(i, 4) + par[6]*pow(i, 3) + par[5]*pow(i, 2) + par[4]*i + par[3];
		fit_cnt++;
	}

	gr_gaus1 = new TGraph(fit_n, fit_x, gaus1_y);
	// cout << "Max element of sig: " << gr_gaus1->GetMaximum();
	gr_gaus1->SetLineColor(1);
	gr_gaus1->SetLineStyle(7);
	gr_polyBG = new TGraph(fit_n, fit_x, polyBG_y);
	gr_polyBG->SetLineColor(7);

	gr_gaus1->Draw("same");
	gr_polyBG->Draw("same");
	
	return 1;

}

Double_t plot_polN(Double_t *par, Int_t polOrd){
	const int fit_n = 180;
	Double_t fit_x[fit_n], polN_y[fit_n];
	int fit_cnt = 0;

	for(double i = 0.4; i < 1.4; i += 0.01 ){
		fit_x[fit_cnt] = i;
		for( int ord = 0; ord <= polOrd; ord ++){
			polN_y[fit_cnt] += par[ord]*pow(i, ord);
		}
		fit_cnt++;
	}

	TGraph *gr_polN = new TGraph(fit_n, fit_x, polN_y);
	gr_polN->SetLineColor(3);
	gr_polN->Draw("same");
	return 1;

}

template<typename T>
double VectorMean(std::vector<T> const& v){
	if(v.empty()){
		return 0;
	}
	return std::accumulate(v.begin(), v.end(), 0.0)/v.size();
}

void list_files(TString directory, vector<TString> &filenames, TString ext){
	const char *dirname = directory.Data();

	TSystemDirectory dir(dirname, dirname);
	TList *files = dir.GetListOfFiles();

	if(files){
		TSystemFile *file;
		TString fname;
		TIter next(files);
		while( (file = (TSystemFile*)next() ) ){
			fname = file->GetName();
			if( !file->IsDirectory() && fname.BeginsWith( ext.Data() )) {
				filenames.push_back( Form("%s/%s", directory.Data(), fname.Data()) );
				// cout << fname.Data() << endl;
			}
		}
	}
}


vector<int> runnum_vec;

//Run info and lookups
TString run_target = "LH2";
int kine = 4;
int sbsfieldscale = 0;

TString experiment = "gmn";

int pass;
double W2_yield, W2_yield_reject;

//Experimental Lookup Parameters
double E_beam = lookup_beam_energy_from_kine(kine); //Electron beam energy (electron energy) in GeV
double SBS_field = sbsfieldscale; //Strength (in percentage) of SBS magnet

double BB_dist, BB_theta, W_mean, W_sigma;
double dx_p, dx_p_sigma, dy_p, dy_p_sigma, dx_n, dx_n_sigma, dy_n, dy_n_sigma, dx_pn_max;
double sig_integral, tf_sig_integral, tf_dx_sig_integral, tf_sig_integral_rej, sig_integral_rej;

// TString rootfile_dir = "/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/11449";
TString rootfile_dir, infile_name, outfile_name;
// TString input_rootfile;

TFile *outfile, *infile, *temp_file, *histo_infile;
TChain *TC = new TChain("T");
vector<TString> master_cut_vec, input_filenames;
TString master_cut_string;

TString elastic_yield_str = "";
TCut master_cut = "";

bool theta_pq_cut = false;
double theta_pq_p_thresh = 0.04;
double theta_pq_n_thresh = 0.04;

//Experimental Constants, Thresholds, cuts, etc DEFINITIONS
const double pi = TMath::Pi();
const double Mp = 0.938272; //Mass of proton [GeV]
const double Mn = 0.939565; //Mass of neutron [GeV]
const double Me = 0.00051; //Mass of electron [GeV]

//SBS Magnet
const Double_t Dgap = 48.0*2.54/100.0; //about 1.22 m
const Double_t maxSBSfield = 1.26; //Tesla
const Double_t SBSdist = 2.25; //m
const Double_t dipGap = 1.22; //m
const Double_t sbsmaxfield = 3.1 * atan( 0.85/(11.0 - 2.25 - 1.22/2.0 ))/0.3/1.22/0.7;

double W2_mean; //Invariant Mass-squared (mean val) {With perfect optics W2 = Mp. Can be calculated run-by-run}
double W2_sigma; //Invariant Mass-squared sigma {Reasonable default/guess. Can be calculated run-by-run from W plot}

//HCal constants and stuff
double tdiff = 510;		//Time difference between BBCal and HCal signals
double tdiff_max = 10;	//Maximum time difference from coincidences through tdctrig cut
double HCal_dist; 	//Distace from HCal face to target chamber
double HCal_theta;		//Theta angle for HCal from downstream beamline
double scint_intersect, x_expected_HCal, y_expected_HCal;
double hcal_y_fmin, hcal_y_fmax, hcal_x_fmin, hcal_x_fmax, ADC_time_min, ADC_time_max;

//Scattered kinematics
double e_prime_theta; //Scattered electron theta angle
double e_prime_phi; //Scattered electron phi angle
double p_el, nu, pp, nucleon_theta, nucleon_phi, E_ep, p_ep, Q2, W, W2, E_pp, E_nucleon, KE_p, dx, dy;

//Static Detector Parameters
const int maxTracks = 1000; // Reasonable limit on tracks to be stored per event
const int maxTdcChan = 10; // Set to accomodate original 5 TDCTrig channels with buffer
// const double hcalheight = -0.2897; // Height of HCal above beamline
double hcalheight;

const Double_t sampfrac = 0.077; 	//Estimate of the sampling fraction from MC
const Int_t kNcell = 288; // Total number of HCal modules
const Int_t kNrows = 24; // Total number of HCal rows
const Int_t kNcols = 12; // Total number of HCal columns
const Int_t kNtrack = 100; // Reasonable max number of tracks per event
const Int_t kNtdc = 1000; // Reasonable max number of tdc signals per event
const Double_t Xi = -2.20; // Distance from beam center to top of HCal in m
const Double_t Xf = 1.47; // Distance from beam center to bottom of HCal in m
const Double_t Yi = -0.853; // Distance from beam center to opposite-beam side of HCal in m
const Double_t Yf = 0.853; // Distance from beam center to beam side of HCal in m

//Static Target Parameters
const double l_tgt = 0.15; // Length of the target (m)
const double rho_tgt = 0.0723; // Density of target (g/cc)
const double rho_Al = 2.7; // Density of aluminum windows (g/cc)
const double cell_diameter = 1.6*2.54; //cm, right now this is a guess
const double Ztgt = 1.0;
const double Atgt = 1.0;
const double Mmol_tgt = 1.008; //g/mol

//For energy-loss correction to beam energy:
const double dEdx_tgt=0.00574; //According to NIST ESTAR, the collisional stopping power of hydrogen is about 5.74 MeV*cm2/g at 2 GeV energy
const double dEdx_Al = 0.0021; //According to NIST ESTAR, the collisional stopping power of Aluminum is about 2.1 MeV*cm2/g between 1-4 GeV
const double uwallthick_LH2 = 0.0145; //cm
const double dwallthick_LH2 = 0.015; //cm
const double cellthick_LH2 = 0.02; //cm, this is a guess;
const double Alshieldthick = 2.54/8.0; //= 1/8 inch * 2.54 cm/inch

double p_recon, nu_recon, E_loss, E_corr, theta_pq_n, theta_pq_p;

int useAlshield = 0;

//Declare vars
Double_t atime[kNcell], row[kNcell], col[kNcell], tdctime[kNcell], cblkid[kNcell], cblke[kNcell];
Double_t nblk, nclus, SH_nclus, PS_nclus, hcal_x, hcal_y, hcal_e;

Double_t par[14], parRej[14];

double bb_tr_p[maxTracks], bb_tr_px[maxTracks], bb_tr_py[maxTracks], bb_tr_pz[maxTracks];
double bb_tr_vx[maxTracks], bb_tr_vy[maxTracks], bb_tr_vz[maxTracks], bb_tr_chi2[maxTracks];
double bb_fp_x[maxTracks], bb_fp_y[maxTracks], bb_fp_th[maxTracks], bb_fp_ph[maxTracks];
double bb_tgt_x[maxTracks], bb_tgt_y[maxTracks], bb_tgt_th[maxTracks], bb_tgt_ph[maxTracks];
double hcal_clusblk_ADC_time[15]; //Maximum number of blocks in a cluster is 15 as per S. Seeds
double bb_tr_n, bb_ps_x, bb_ps_y, bb_ps_e, bb_sh_x, bb_sh_y, bb_sh_e;
double bb_ntracks[maxTracks];

Double_t TDCT_id[kNtdc], TDCT_tdc[kNtdc], hodo_tmean[kNtdc]; 
Int_t TDCTndata;

Long64_t nTCevents, Nevents;

//INITIALIZE ALL HISTOGRAMS:
TH1D *hin_bb_gem_Ntracks, *h_W2_full_rand;
TH1D *h_atime, *h_W, *h_W2, *h_W2copy, *h_W2recon, *h_KE_p, *h_KE_low, *h_Diff, *h_X, *h_Y, *h_E_eloss, *h_hcal_clusblk_ADC_time;
TH1D *h_Wfull, *h_W2full, *h_W2_residual;
TH1D *h_W_cut, *h_W_fcut, *h_vz_cut;
TH1D *h_Q2, *h_E_ep, *h_E_pp;
TH1D *h_dy, *h_dy_cut, *h_dy_wcut, *h_dx, *h_dx_cut, *h_dx_wcut, *h_dx_fcut, *h_dx_wcut_fcut, *h_dy_wcut_fcut;
TH1D *h_dx_w2cut, *h_dx_w2cut_fcut;
TH1D *h_Nevents;
 
TH2D *h_E_ecorr_vs_vert;
TH2D *h_dxdy, *h_dxdy_cut, *h_dxdy_wcut, *h_dxdy_ncut, *h_dxdy_pcut, *h_dxdy_fcut, *h_dxdy_wcut_fcut;
TH2D *h_xy, *h_xy_cut, *h_xy_fcut, *h_xy_cut_p, *h_xy_cut_n, *h_PAngleCorr_theta, *h_PAngleCorr_phi;

double n_integral, p_integral, n_center, n_sigma, p_center, p_sigma, dx_bg_max;
double p_Beam, E_loss_outgoing;
double Eloss, E_beam_final, p_corr;
int n_counts, p_counts, elastic_yield;
int n_hcal_clusblk_atime_cut = 0;
int N_gem_tracks_found = 0;
int N_bb_tracks_found = 0;

TString pq_cut_String = "";

TF1 *tf1_W2_bg_pol8, *tf1_W2_sig, *tf1_W2_fullFit, *tf1_W2_bg_rej, *tf1_W2_bg_rej_full, *tf1_W2_bg_rej_init, *tf1_W2_bg_rej_init, *tf1_W2_bg_pol8_rej_full, *tf1_W2_bg_pol8_rej, *tf1_W2_sig_rej, *tf1_W2_fullFit_rej;
TF1 *tf1_dx_fullFit, *tf1_dx_bg_pol8, *tf1_dx_sig;
TH1 *h_tf_W2_bg, *h_tf_W2_sig, *h_tf_W2_bg_rej, *h_tf_W2_sig_rej, *h_tf_W2_full;
TH1 *h_tf_dx_bg, *h_tf_dx_sig, *h_tf_dx_full;

void det_eff(){

	auto total_time_start = high_resolution_clock::now();
	TStopwatch *StopWatch = new TStopwatch();

	gStyle->SetPalette(55);
	cout << "--------------------------------------" << endl;
	cout << "Analysis started. " << endl;
	cout << "--------------------------------------" << endl;

	cout << "Run parameters: " << endl;
	cout << "Kinematic: SBS" << kine << endl;
	cout << "Target: " << run_target.Data() << endl;
	cout << "Beam Energy: " << E_beam << endl;
	cout << "SBS Field: " << sbsfieldscale << "%" << endl;
	cout << "-----------------------------------" << endl;
	cout << "BB angle [deg]: " << (180/pi)*BB_theta << endl;
	cout << "SBS angle: " << lookup_SBS_angle_by_kine(kine, "deg") << endl;
	cout << "HCal angle [deg]: " << (180/pi)*HCal_theta << endl;
	cout << "HCal distance: " << HCal_dist << endl;
	cout << "-----------------------------------" << endl << endl;

	if( kine == 4 ){
		pass = 0;
		hcalheight = -0.312479; // Height of HCal above beamline
		//offset: -1.37521e-01
	}
	if( kine == 8 ){
		pass = 1;
		hcalheight = -0.450; // Height of HCal above beamline
	}
	
	outfile_name = Form("rootfiles/det_eff_histos_SBS%i_%s_mag%i.root", kine, run_target.Data(), sbsfieldscale);

	if( !fit_only ){

		BB_dist = lookup_BB_dist_by_kine(kine);
		BB_theta = lookup_BB_angle_by_kine(kine, "rad");
		HCal_dist = lookup_HCal_dist_by_kine(kine);
		HCal_theta = lookup_HCal_angle_by_kine(kine, "rad");
		// W_mean = lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "W_mean");
		// W_sigma = lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "W_sigma");

		if( kine == 4 ){
			W_mean = 0.918;
			W_sigma = 0.053;
			W2_mean = 0.835;
			W2_sigma = 0.1;
		}

		hin_bb_gem_Ntracks = new TH1D("hin_bb_gem_Ntracks", "Number of tracks found", 11, -0.5, 10.5);
		// ADC_time_min = lookup_ADC_time_cut(run_target, kine, sbsfieldscale, "ADC_time_min");
		// ADC_time_max = lookup_ADC_time_cut(run_target, kine, sbsfieldscale, "ADC_time_max");

		if( use_parsed ){
			cout << endl << "-------- Parsed run mode --------" << endl;
			rootfile_dir = "/lustre19/expphy/volatile/halla/sbs/jboyd/analysis_rootfiles/jboyd_parsed";
			infile_name = Form("gmn_parsed_%s_SBS%i_mag%i.root", run_target.Data(), kine, sbsfieldscale);
			cout << endl << "Adding files to TChain..." << endl;
			TC->Add(Form("%s/%s", rootfile_dir.Data(), infile_name.Data()));
			cout << "Finished adding to TChain. " << endl;
		}

		else{
			cout << endl << "-------- Raw Un-parsed run mode --------" << endl;
			rootfile_dir = Form("/w/halla-scshelf2102/sbs/sbs-gmn/pass%i/SBS%i/%s/rootfiles", pass, kine, run_target.Data());
			int num_runs_config = lookup_parsed_runs_cnt( run_target, kine, sbsfieldscale); 
			cout << "Adding runnum files to TChain..." << endl;
			for( int run = 0; run < num_runs_config; run++ ){
				runnum_vec.push_back(lookup_parsed_runnums( run_target, kine, sbsfieldscale, run ));
				infile_name = Form("e1209019_fullreplay_%i_stream0_seg*", runnum_vec[run]);
				cout << runnum_vec[run] << "  ";
				TC->Add(Form("%s/%s", rootfile_dir.Data(), infile_name.Data()));
				list_files( rootfile_dir, input_filenames, infile_name.ReplaceAll("*", "") );
			}
			cout << "Finished adding " << num_runs_config << " files to the TChain. " << endl << endl;

			cout << "Pulling histograms from files... " << endl;

			for( int file = 0; file < input_filenames.size(); file++ ){
				temp_file = new TFile(input_filenames[file].Data(), "READ");

				hin_bb_gem_Ntracks->Add( static_cast<TH1D*>(temp_file->Get("hbb_gem_Ntracks")) );	

			}
		}
		
		nTCevents = TC->GetEntries();
		cout << "--------------------------------------" << endl;
		cout << "--------------------------------------" << endl << endl;
		cout << "Number of events in TChain: " << nTCevents << endl  << endl;
		cout << "--------------------------------------" << endl;
		cout << "--------------------------------------" << endl << endl;

		TC->SetBranchStatus( "*", 0 );

		// HCal
		TC->SetBranchStatus( "sbs.hcal.x", 1 );
		TC->SetBranchStatus( "sbs.hcal.y", 1 );
		TC->SetBranchStatus( "sbs.hcal.e", 1 );
		TC->SetBranchStatus( "sbs.hcal.nclus", 1);

		// BB track
		TC->SetBranchStatus( "bb.tr.chi2", 1 );
		TC->SetBranchStatus( "bb.tr.n", 1 );
		TC->SetBranchStatus( "bb.tr.px", 1 );
		TC->SetBranchStatus( "bb.tr.py", 1 );
		TC->SetBranchStatus( "bb.tr.pz", 1 );    
		TC->SetBranchStatus( "bb.tr.p", 1 );
		TC->SetBranchStatus( "bb.tr.vx", 1 );
		TC->SetBranchStatus( "bb.tr.vy", 1 );
		TC->SetBranchStatus( "bb.tr.vz", 1 );
		TC->SetBranchStatus( "bb.tr.r_x", 1 );
		TC->SetBranchStatus( "bb.tr.r_y", 1 );
		TC->SetBranchStatus( "bb.tr.r_th", 1 );
		TC->SetBranchStatus( "bb.tr.r_ph", 1 );
		TC->SetBranchStatus( "bb.tr.tg_x", 1 );
		TC->SetBranchStatus( "bb.tr.tg_y", 1 );
		TC->SetBranchStatus( "bb.tr.tg_th", 1 );
		TC->SetBranchStatus( "bb.tr.tg_ph", 1 );
		TC->SetBranchStatus( "bb.gem.track.ntrack", 1);
		TC->SetBranchStatus( "bb.gem.track.nhits", 1);
		TC->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1);

		// BBCal shower preshower
		TC->SetBranchStatus( "bb.ps.e", 1 );
		TC->SetBranchStatus( "bb.ps.x", 1 );
		TC->SetBranchStatus( "bb.ps.y", 1 );
		TC->SetBranchStatus( "bb.sh.e", 1 );
		TC->SetBranchStatus( "bb.sh.x", 1 );
		TC->SetBranchStatus( "bb.sh.y", 1 );
		TC->SetBranchStatus( "bb.sh.nclus", 1 );
		TC->SetBranchStatus( "bb.ps.nclus", 1 );

		// Trigger TDC
		TC->SetBranchStatus( "bb.tdctrig.tdc", 1 );
		TC->SetBranchStatus( "bb.tdctrig.tdcelemID", 1 );
		TC->SetBranchStatus( "Ndata.bb.tdctrig.tdcelemID", 1 );

	// Set BRANCH ADDRESSES
		// HCal
		TC->SetBranchAddress( "sbs.hcal.x", &hcal_x );
		TC->SetBranchAddress( "sbs.hcal.y", &hcal_y );
		TC->SetBranchAddress( "sbs.hcal.e", &hcal_e );
		TC->SetBranchAddress( "sbs.hcal.nclus", &nclus );
		TC->SetBranchAddress( "sbs.hcal.clus_blk.atime", &hcal_clusblk_ADC_time );

		// BB track
		TC->SetBranchAddress( "bb.tr.chi2", bb_tr_chi2 );
		TC->SetBranchAddress( "bb.tr.n", &bb_tr_n );
		TC->SetBranchAddress( "bb.tr.px", bb_tr_px );
		TC->SetBranchAddress( "bb.tr.py", bb_tr_py );
		TC->SetBranchAddress( "bb.tr.pz", bb_tr_pz );
		TC->SetBranchAddress( "bb.tr.p", bb_tr_p );
		TC->SetBranchAddress( "bb.tr.vx", bb_tr_vx );
		TC->SetBranchAddress( "bb.tr.vy", bb_tr_vy );
		TC->SetBranchAddress( "bb.tr.vz", bb_tr_vz );
		TC->SetBranchAddress( "bb.tr.r_x", bb_fp_x );
		TC->SetBranchAddress( "bb.tr.r_y", bb_fp_y );
		TC->SetBranchAddress( "bb.tr.r_th", bb_fp_th );
		TC->SetBranchAddress( "bb.tr.r_ph", bb_fp_ph );
		TC->SetBranchAddress( "bb.tr.tg_x", bb_tgt_x );
		TC->SetBranchAddress( "bb.tr.tg_y", bb_tgt_y );
		TC->SetBranchAddress( "bb.tr.tg_th", bb_tgt_th );
		TC->SetBranchAddress( "bb.tr.tg_ph", bb_tgt_ph );
		TC->SetBranchAddress( "bb.gem.track.ntrack", bb_ntracks );

		// BBCal shower preshower
		TC->SetBranchAddress( "bb.ps.e", &bb_ps_e );
		TC->SetBranchAddress( "bb.ps.x", &bb_ps_x );
		TC->SetBranchAddress( "bb.ps.y", &bb_ps_y );
		TC->SetBranchAddress( "bb.sh.e", &bb_sh_e );
		TC->SetBranchAddress( "bb.sh.x", &bb_sh_x );
		TC->SetBranchAddress( "bb.sh.y", &bb_sh_y );
		TC->SetBranchAddress( "bb.sh.nclus", &SH_nclus );
		TC->SetBranchAddress( "bb.ps.nclus", &PS_nclus );

		// Trigger TDC
		TC->SetBranchAddress( "bb.tdctrig.tdcelemID", TDCT_id );
		TC->SetBranchAddress( "bb.tdctrig.tdc", TDCT_tdc );
		TC->SetBranchAddress( "Ndata.bb.tdctrig.tdcelemID", &TDCTndata );
		cout << " done. " << endl;
		cout << "--------------------------------------" << endl;


		if( use_bbcal_cuts ){
			cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << endl;
			cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
			cout << "		Using BBCal Cuts" << endl;
			cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << endl;
			cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
			master_cut_vec = {
				"sbs.hcal.nclus>0",
				"bb.ps.nclus>0",
				"bb.sh.nclus>0",
				"abs(bb.tr.vz[0])<0.076",
				"bb.gem.track.nhits[0]>4",
				"bb.tr.n==1",	
				"bb.ps.e>0.15"	
			};

			for(size_t cut = 0; cut < master_cut_vec.size(); cut++){
				if(cut == master_cut_vec.size() - 1){
					master_cut_string.Append(Form("%s", master_cut_vec[cut].Data()));
				}
				else{
					master_cut_string.Append(Form("%s%s", master_cut_vec[cut].Data(), "&&"));
				}
			}

			master_cut = Form("%s", master_cut_string.Data());
			cout << "--------------------------------------" << endl;
			cout << "--------------------------------------" << endl << endl;
			cout << "Number of RAW entries: " << TC->GetEntries() << endl << endl;
			cout << "--------------------------------------" << endl;
			cout << "--------------------------------------" << endl;
			cout << "Applying Master Cut: " << endl;
			cout << master_cut << endl;

		}
		else{
			cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << endl;
			cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
			cout << "		No BBCal (master) Cuts applied" << endl;
			cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << endl;
			cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
		}

		TEventList *ev_list = new TEventList("ev_list", "Elastic Events List");
		TC->Draw(">>ev_list", master_cut);
		Nevents = ev_list->GetN();

		outfile = new TFile(outfile_name.Data(), "RECREATE");

		h_E_eloss = new TH1D("E_eloss", Form("Scattered Electron Energy Loss in Target - SBS%i %i, %s", kine, sbsfieldscale, run_target.Data()), 500, 0.0, (0.1)*E_beam);
		h_E_ecorr_vs_vert = new TH2D("h_E_ecorr_vs_vert", Form("Corrected Beam Energy vs Vertex - SBS%i %i, %s; E_{e} (GeV); Z_{vertex} (m)", kine, sbsfieldscale, run_target.Data()), 250, -0.125, 0.125, 500, 0, 0.001);
		h_Q2 = new TH1D("h_Q2", Form("Momentum Transfer Q^2 - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 250, 0.5, 3.0);
		h_E_ep = new TH1D("h_E_ep", Form("Scattered Electron Energy - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 500, 0.0, 1.5*E_beam);
		h_E_pp = new TH1D("h_E_pp", Form("Scattered Proton Energy - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 500, 0.0, 1.5*E_beam);
		h_W = new TH1D("h_W", Form("Invariant Mass W - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 280, 0.0, 1.4);
		h_W_cut = new TH1D("h_W_cut", Form("Invariant Mass W (Coin & Vert Cuts) - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 600, 0.0, 3.0);
		h_W_fcut = new TH1D("h_W_fcut", Form("Invariant Mass W (Fiduc. Cuts) - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 600, 0.0, 3.0);
		h_W2 = new TH1D("h_W2", Form("Invariant Mass W^2 - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 280, 0.0, 1.4);

		h_Wfull = new TH1D("h_Wfull", Form("Invariant Mass W Full Range- SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 600, 0.0, 3.0);
		h_W2full = new TH1D("h_W2full", Form("Invariant Mass W^2 Full Range - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 600, 0.0, 3.0);


		h_KE_p = new TH1D("h_KE_p", Form("Scattered Proton Kinetic Energy - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 500, 0.0, 1.5*E_beam);

		h_dx = new TH1D("h_dx",Form("dx (NO CUTS) - SBS%i %i, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 1000, -2.5, 2.5);

		h_dx_cut = new TH1D("h_dx_cut",Form("dx (Basic CUTS) - SBS%i %i, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5,2.5);
		h_dx_wcut = new TH1D("h_dx_wcut",Form("dx (W cut) - SBS%i %i, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5,2.5);
		h_dx_w2cut = new TH1D("h_dx_w2cut",Form("dx (W^2 cut) - SBS%i %i, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5,2.5);
		h_dx_fcut = new TH1D("h_dx_fcut",Form("dx (f cut) - SBS%i %i, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5,2.5);
		h_dx_wcut_fcut = new TH1D("h_dx_wcut_fcut",Form("dx (W & Fiduc. Cuts) - SBS%i %i, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5,2.5);
		h_dx_w2cut_fcut = new TH1D("h_dx_w2cut_fcut",Form("dx (W^2 & Fiduc. Cuts) - SBS%i %i, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5,2.5);
		h_dy = new TH1D("h_dy",Form("dy (NO CUTS) - SBS%i %i, %s; y_{HCal} - y_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -1.25,1.25);
		h_dy_cut = new TH1D("h_dy_cut",Form("dy (Basic Cuts) - SBS%i %i, %s; y_{HCal} - y_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -1.25,1.25);  
		h_dy_wcut = new TH1D("h_dy_wcut",Form("dy (W Cuts) - SBS%i %i, %s; y_{HCal} - y_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -1.25,1.25);
		h_dy_wcut_fcut = new TH1D("h_dy_wcut_fcut",Form("dy (W & Fiduc. Cuts) - SBS%i %i, %s; y_{HCal} - y_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -1.25,1.25);    

		h_dxdy = new TH2D("h_dxdy", Form("Hadron Spot(s) on HCal (NO CUTS) - SBS%i %i, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), 300, -1.5, 1.5, 500, -2.5, 2.5 );
		h_dxdy_wcut = new TH2D("h_dxdy_wcut", Form("Hadron Spot(s) on HCal (W cut) - SBS%i %i, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), 300, -1.5, 1.5, 500, -2.5, 2.5 );
		h_dxdy_wcut_fcut = new TH2D("h_dxdy_wcut_fcut", Form("Hadron Spot(s) on HCal (W & Fiduc. Cuts) - SBS%i %i, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), 300, -1.5, 1.5, 500, -2.5, 2.5 );
		h_dxdy_cut = new TH2D("h_dxdy_cut", Form("Hadron Spot(s) on HCal (Basic cuts) - SBS%i %i, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), 300, -1.5, 1.5, 500, -2.5, 2.5 );
		h_dxdy_ncut = new TH2D("h_dxdy_ncut", Form("Hadron Spot(s) on HCal (n cut) - SBS%i %i, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), 250, -1.25, 1.25, 250, -2.5, 2.5 );
		h_dxdy_pcut = new TH2D("h_dxdy_pcut", Form("Hadron Spot(s) on HCal (p cut) - SBS%i %i, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), 250, -1.25, 1.25, 250, -2.5, 2.5 );
		h_dxdy_fcut = new TH2D("h_dxdy_fcut", Form("Hadron Spot(s) on HCal (f cut) - SBS%i %i, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), 250, -1.25, 1.25, 250, -2.5, 2.5 );
		h_xy = new TH2D("h_xy",Form("HCal Hadron Spots (x, y) (NO CUTS) - SBS%i %i, %s;y_{HCal} (m); x_{HCal} (m)", kine, sbsfieldscale, run_target.Data()),12,-0.9,0.9,24,-2.165,1.435);
		h_xy_cut = new TH2D("h_xy_cut", Form("HCal Hadron Spots (x, y) (BASIC CUTS) - SBS%i %i, %s;y_{HCal} (m); x_{HCal} (m)", kine, sbsfieldscale, run_target.Data()),12,-0.9,0.9,24,-2.165,1.435);
		h_xy_fcut = new TH2D("h_xy_fcut", Form("HCal Hadron Spots (x, y) (Fiduc. CUTS) - SBS%i %i, %s;y_{HCal} (m); x_{HCal} (m)", kine, sbsfieldscale, run_target.Data()),12,-0.9,0.9,24,-2.165,1.435);
		h_xy_cut_p = new TH2D("h_xy_cut_p", Form("HCal Hadron Spots (x, y) (p CUT) - SBS%i %i, %s;y_{HCal} (m); x_{HCal} (m)", kine, sbsfieldscale, run_target.Data()),12,-0.9,0.9,24,-2.165,1.435);
		h_xy_cut_n = new TH2D("h_xy_cut_n", Form("HCal Hadron Spots (x, y) (n CUT) - SBS%i %i, %s;y_{HCal} (m); x_{HCal} (m)", kine, sbsfieldscale, run_target.Data()),12,-0.9,0.9,24,-2.165,1.435);

		h_hcal_clusblk_ADC_time = new TH1D("h_hcal_clusblk_ADC_time", Form("ADC time of the highest energy block in the largest cluster - SBS%i %i, %s; ADC Time (ns)", kine, sbsfieldscale, run_target.Data()), 300, -100, 200);
		h_PAngleCorr_theta = new TH2D( "h_PAngCorr_theta",Form("BB theta vs HCal theta - SBS%i %i, %s", kine, sbsfieldscale, run_target.Data()), 200, 0.55, 0.75, 300, 0.35, 0.65 );
		h_PAngleCorr_phi = new TH2D( "h_PAngCorr_phi",Form("BB phi vs HCal phi - SBS%i %i, %s", kine, sbsfieldscale, run_target.Data()), 500, -0.4, 0.1, 500, 2.7, 3.2 );
		h_vz_cut = new TH1D("h_vz_cut",Form("BB phi vs HCal phi - SBS%i %i, %s; vertex z (m);", kine, sbsfieldscale, run_target.Data()), 250,-0.125,0.125);


		cout << "--------------------------------------" << endl;
		cout << "Number of events to analyze: " << Nevents << endl;
		cout << "--------------------------------------" << endl;
		cout << "--------------------------------------" << endl;
		cout << "Starting analysis loop on events..... " << endl;

	//Basic Energy calcs
		p_Beam = E_beam/(1.0 + E_beam/Mp*(1.0 - cos(BB_theta)));

		E_loss_outgoing = cell_diameter/2.0/sin(BB_theta)*rho_tgt*dEdx_tgt; //Should be about 1 MeV
		if( useAlshield !=0 ) E_loss_outgoing += Alshieldthick*rho_Al*dEdx_Al;

	//FIDUCIAL CUT:
		hcal_y_fmin = -0.75;
		hcal_y_fmax = 0.75;
		hcal_x_fmin = -2.015;
		hcal_x_fmax = 1.285;

		elastic_yield = 0;

		int watch_cnt = 0;	
		int five_percent = int(0.05*Nevents);
		vector<double> time_for_five;
		double average_time = 0.0, time_remaining;
		StopWatch->Start();
		
		for(Long64_t nevent = 0; nevent < Nevents; nevent++){
			TC->GetEntry( ev_list->GetEntry( nevent ));

			if( nevent%five_percent == 0){		
				StopWatch->Stop();

				if( watch_cnt == 0){
					cout << "Evt: " << nevent <<"/" << Nevents << Form("(%.0f/100%%)", 100.0*double(1.0*nevent/Nevents)) << ". Elastic yield = " << elastic_yield << ". "<< endl;
				}

				if( watch_cnt > 0 ){
					time_for_five.push_back(StopWatch->RealTime());	
					average_time = VectorMean(time_for_five);
					// cout << "average time for 5 = " << average_time << endl;
					time_remaining = average_time*( 1.0 - double(nevent)/double(Nevents));
					cout << "Evt: " << nevent <<"/" << Nevents << Form("(%.0f/100%%)", 100.0*double(1.0*nevent/Nevents)) << ". Elastic yield = " << elastic_yield << ". Time left: " << time_remaining << endl;
				}
				watch_cnt++;
				StopWatch->Reset();
				StopWatch->Continue();
			}

		// PULL NUMBER OF TRACKS FROM BB AND FROM GEMS
			N_gem_tracks_found += bb_ntracks[0];
			N_bb_tracks_found += bb_tr_n;

		//Timing stuff
			// h_hcal_clusblk_ADC_time->Fill(hcal_clusblk_ADC_time[0]);

			// bool b_hcal_clusblk_ADC_cut = false;

			// if( hcal_clusblk_ADC_time[0] < ADC_time_min || hcal_clusblk_ADC_time[0] > ADC_time_max ){
			// 	b_hcal_clusblk_ADC_cut = true;
			// 	n_hcal_clusblk_atime_cut++;
			// 	continue;
			// }

		//Sanity check
			if( use_bbcal_cuts ){
		      	if( (int)bb_tr_n!=1 ){
		      		cout << "**************************************************************" << endl;
		      		cout << "--------------------------------------------------------------" << endl;
		      		cout << endl << endl << "WARNING: Total tracks not as expected from global cut. Check globalcut for errors." << endl << endl;
		      		cout << "--------------------------------------------------------------" << endl;
		      		cout << "**************************************************************" << endl;
		      	} 			
			}

		    if( correct_beam_energy ){
	      		Eloss = (bb_tr_vz[0]+l_tgt/2.0) * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //approximately 3 MeV
	      		h_E_eloss->Fill( Eloss );

	      		E_beam_final = E_beam - Eloss;
	      		h_E_ecorr_vs_vert->Fill( bb_tr_vz[0], E_beam_final);
	      	}
	      	if( !correct_beam_energy){
	      		Eloss = (bb_tr_vz[0]+l_tgt/2.0) * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //approximately 3 MeV
	      		h_E_eloss->Fill( Eloss );

	      		E_beam_final = E_beam;
	      		h_E_ecorr_vs_vert->Fill( bb_tr_vz[0], E_beam_final);
	      	}

	    ////Corrections
		    //Correct the beam energy with energy loss in target using vertex position
		    Double_t E_loss = (bb_tr_vz[0]+l_tgt/2.0) * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //approximately 3 MeV
		    Double_t E_corr = E_beam - Eloss;	

	      	p_corr = bb_tr_p[0] - E_loss_outgoing; //Neglecting mass of e'


	    //Proceed only if at least one track exists in BB arm - lowest chi2 track always first element
	      	if( bb_tr_n > 1){
	      		continue;
	      	}

	//-------------
			p_Beam = E_beam/(1.0 + E_beam/Mp*(1.0 - cos(BB_theta)));
	      	
	      	e_prime_theta = acos( bb_tr_pz[0]/bb_tr_p[0] ); //Ucorrected track momenutm to reconstruct e' theta
	      	e_prime_phi = atan2( bb_tr_py[0], bb_tr_px[0]);

	      	TVector3 vertex( 0, 0, bb_tr_vz[0] ); // z location of vertex in hall coordinates
			TLorentzVector P_beam( 0, 0, E_beam_final, E_beam_final ); //Mass of e negligable
			TLorentzVector k_prime( bb_tr_px[0], bb_tr_py[0], bb_tr_pz[0], bb_tr_p[0] );
			TLorentzVector P_targ( 0, 0, 0, Mp );

			p_el = E_beam_final/( 1.0+E_beam_final/Mp*( 1.0-cos(e_prime_theta) ) );
			nu = E_beam_final - bb_tr_p[0];
			pp = sqrt( pow(nu,2)+2.*Mp*nu );
			nucleon_phi = e_prime_phi + pi; //assume coplanarity
			//double thetanucleon = acos( (E_corr - BBtr_p[0]*cos(etheta))/pp ); //use elastic constraint on nucleon kinematics
			nucleon_theta = acos( (E_beam_final - bb_tr_pz[0])/pp ); //use elastic constraint on nucleon kinematics

			TVector3 pNhat( sin(nucleon_theta)*cos(nucleon_phi), sin(nucleon_theta)*sin(nucleon_phi), cos(nucleon_theta) );

			//Define HCal coordinate system
			TVector3 HCAL_zaxis( sin(-HCal_theta ), 0, cos(-HCal_theta) );
			TVector3 HCAL_xaxis( 0, -1, 0 );
			TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();

			TVector3 HCAL_origin = HCal_dist*HCAL_zaxis + hcalheight*HCAL_xaxis;

			TVector3 HCAL_pos = HCAL_origin + (hcal_x*HCAL_xaxis) + (hcal_y*HCAL_yaxis);

			//Define intersection points for hadron vector
			scint_intersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / (pNhat.Dot( HCAL_zaxis ) );
			TVector3 HCAL_intersect = vertex + scint_intersect * pNhat;

			//Define the expected position of hadron on HCal from BB track
			x_expected_HCal = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis );
			y_expected_HCal = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );

	//--------------------------------------------------------------
	//Calculate theta pq variables
	      	//Reconstructed momentum, corrected for mean loss exiting the target
			p_recon = bb_tr_p[0] + E_loss_outgoing; 

			TLorentzVector k_prime_recon(p_recon*bb_tr_px[0]/bb_tr_p[0], p_recon*bb_tr_py[0]/bb_tr_p[0], p_recon*bb_tr_pz[0]/bb_tr_p[0], p_recon);
			TLorentzVector q_recon = P_beam - k_prime_recon;
			TVector3 qvec_recon = q_recon.Vect();

		//Calculate q vector as beam momentum - scattered k
			TLorentzVector q = P_beam - k_prime;

		//Expected neutron direction
			TVector3 Neutron_Direction = (HCAL_pos - vertex).Unit();

		//Expected proton direction
			//Need to incorporate deflection due to SBS magnet
			double Bdl = sbsfieldscale*maxSBSfield*Dgap;
			double Proton_Deflection = tan( 0.3*Bdl/qvec_recon.Mag() )*(HCal_dist - (SBSdist + Dgap/2.0) ); 

			TVector3 Proton_Direction = (HCAL_pos + Proton_Deflection*HCAL_xaxis - vertex).Unit();

			theta_pq_n = acos( Neutron_Direction.Dot( qvec_recon.Unit() ) );
			theta_pq_p = acos( Proton_Direction.Dot( qvec_recon.Unit() ) );

			if( theta_pq_cut ){
				if( run_target == "LD2" ){
					if( theta_pq_p < theta_pq_p_thresh && theta_pq_n < theta_pq_n_thresh){
						continue;
					}
					// if( theta_pq_n > theta_pq_n_thresh ){
					// 	continue;
					// }
				}
			}

	//--------------------------------------------------------------

			TLorentzVector P_gammaN = P_targ + q; //(-px, -py, ebeam - pz, Mp + ebeam - p)

			E_ep = sqrt( pow(Me,2) + pow(bb_tr_p[0],2) ); // Obtain the scattered electron energy
			h_E_ep->Fill( E_ep );

			p_ep = bb_tr_p[0];

			Q2 = 2*E_beam_final*E_ep*( 1-(bb_tr_pz[0]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta
			h_Q2->Fill( Q2 );

			//Get invariant mass transfer W from the four-momentum of the scattered nucleon
			W = P_gammaN.M();
			W2 = pow(W, 2);

		//fill work histograms
			h_W->Fill( W );
			h_Wfull->Fill( W );

			if( W2 < ( (h_W2->GetXaxis()->GetNbins())*(h_W2->GetXaxis()->GetBinWidth(1)) ) ){
				h_W2->Fill( W2 );	
			}
			h_W2full->Fill( W2 );

			if( !calc_W_only ){

				//Use the electron kinematics to predict the proton momedntum assuming elastic scattering on free proton at rest (will need to correct for fermi motion):
				E_pp = nu + Mp; // Get energy of the proton
				E_nucleon = sqrt(pow(pp,2)+pow(Mp,2)); // Check on E_pp, same
				h_E_pp->Fill( E_pp ); // Fill histogram

				KE_p = nu; // For elastics
				h_KE_p->Fill( KE_p );

				dx = hcal_x - x_expected_HCal;
				dy = hcal_y - y_expected_HCal;

				//Resolve the hadron spots without cuts
				h_dx->Fill( dx );
				h_dy->Fill( dy );
				h_dxdy->Fill( dy, dx );
				h_xy->Fill( hcal_y, hcal_x );

			// Coincidence timing cut and vertex cut to resolve W well
				// if( fabs(diff - tdiff)<tdiff_max && fabs(bb_tr_vz[0])<=0.075 ){
				// 	h_W_cut->Fill( W );
				// } 

				// Preliminary HCal projections with single cut on W
				if( fabs(W - W_mean) < 4*W_sigma ){
				// if( ( W > (W_mean - 0.75*W_sigma) ) && ( W < (W_mean + 0.75*W_sigma) ) ){
				// if( W2 < 1.05 && W2 > 0.85){
					h_dx_wcut->Fill( dx );
					h_dy_wcut->Fill ( dy );
					h_dxdy_wcut->Fill( dy, dx );
					h_W_cut->Fill( W );
				}
				if( fabs(W2 - W2_mean) < 4*W2_sigma ){
					h_dx_w2cut->Fill( dx );
				}
				//Populate position histograms with cuts
				h_dxdy_cut->Fill( dy, dx );
				h_dx_cut->Fill( dx );
				h_dy_cut->Fill( dy );

				//Populate BB/HCal correlation histograms from elastics
				h_PAngleCorr_phi->Fill( e_prime_phi, nucleon_phi );
				h_PAngleCorr_theta->Fill( e_prime_theta, nucleon_theta );

				//Fill vertex position histogram for cut on tracks
		    	h_vz_cut->Fill( bb_tr_vz[0] );



				//Check "elastic" events on center HCal for id with spot checks
				bool HCal_on = false;
				bool is_p = false;
				bool is_n = false;

			//FIDUCIAL Cut

				if( fiducial_cut ){
					if( hcal_y>hcal_y_fmin && hcal_y<hcal_y_fmax && hcal_x>hcal_x_fmin && hcal_x<hcal_x_fmax ){
						HCal_on = true;
					}

					apply_fcut = ((y_expected_HCal - dy_p_sigma) > hcal_y_fmin) && ((y_expected_HCal + dy_p_sigma) < hcal_y_fmax) && ((x_expected_HCal - dx_pn_max - dx_p_sigma) > hcal_x_fmin) && ((x_expected_HCal + dx_p_sigma) < hcal_x_fmax);

					if( ( pow( (hcal_x - x_expected_HCal - dx_p)/dx_p_sigma,2) + pow( (hcal_y - y_expected_HCal - dy_p)/dy_p_sigma,2) ) <= pow(1.5,2) ){
						is_p = true;
					}

					if( ( pow( (hcal_x - x_expected_HCal - dx_n)/dx_n_sigma,2) + pow( (hcal_y - y_expected_HCal - dy_n)/dy_n_sigma,2) ) <= pow(1.5,2) ){
						is_n = true;
					}

			//Fill respective histograms for these checks.
					if( HCal_on && is_n && apply_fcut ) h_dxdy_ncut->Fill( dy, dx );
					if( HCal_on && is_p && apply_fcut ) h_dxdy_pcut->Fill( dy, dx );

			//----------neutron
					if( HCal_on && is_n && apply_fcut ){
						if( (hcal_x - dx_pn_max )>hcal_x_fmin ){
							h_dxdy_fcut->Fill( dy, dx );
							h_dx_fcut->Fill( dx );
							h_W_fcut->Fill( W );
							h_xy_fcut->Fill( hcal_y, hcal_x );
							h_xy_cut_n->Fill( hcal_y, hcal_x );

							//including work cut
							if( fabs(W - W_mean) < 4*W_sigma ){
								h_dx_wcut_fcut->Fill( dx );
								h_dy_wcut_fcut->Fill ( dy );
								h_dxdy_wcut_fcut->Fill( dy, dx );
							}
							if( fabs(W2 - W2_mean) < 4*W2_sigma ){
								h_dx_w2cut_fcut->Fill(dx);
							}

							elastic_yield++;
						}
					}
			//----------proton
					else if( HCal_on && is_p && apply_fcut ){
						if( (hcal_x + dx_pn_max)<hcal_x_fmax ){
							h_dxdy_fcut->Fill( dy, dx );
							h_dx_fcut->Fill( dx );
							h_W_fcut->Fill( W );
							h_xy_fcut->Fill( hcal_y, hcal_x );
							h_xy_cut_p->Fill( hcal_y, hcal_x );

							//including work cut
							if( fabs(W - W_mean) < W_sigma ){
								h_dx_wcut_fcut->Fill( dx );
								h_dy_wcut_fcut->Fill ( dy );
								h_dxdy_wcut_fcut->Fill( dy, dx );
							}

							elastic_yield++;
						}
					}
				}
			//END OF FIDUCIAL Cut
				//Still should count elastic yields if we got this far.....
				if( !fiducial_cut ){
					elastic_yield++;
				}
			}

			if( calc_W_only ){
				elastic_yield++;
			}
	//END OF ENTRY LOOP
		}	

		cout << "---------------------------------------" << endl;
		cout << "-----Finished going through events-----" << endl;
		cout << "---------------------------------------" << endl;
		outfile->Write();	
	}

//--------------------------------------------------
//------------------- FITS -----------------------
//-----------------------------------------------

	histo_infile = new TFile(outfile_name.Data(), "READ");

	if( fit_only ){
		h_W2 = static_cast<TH1D*>(histo_infile->Get("h_W2"));
		h_W2copy = static_cast<TH1D*>(histo_infile->Get("h_W2copy"));
		h_dx = static_cast<TH1D*>(histo_infile->Get("h_dx"));

		h_W2_full_rand = new TH1D("h_W2_full_rand", Form("Invariant Mass W^2 Full Histo Filled Randomly- SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 280, 0.0, 1.4);
		h_W2_residual = new TH1D("h_W2_residual", Form("Invariant Mass W^2 Fit Residual- SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 280, 0.0, 1.4);

		TCanvas *c_W2_fullFit = new TCanvas("c_W2_fullFit", "c_W2_fullFit", 600, 500);
		h_W2->Draw();
		tf1_W2_fullFit = new TF1("tf1_W2_fullFit", fullFitFunction, 0, 1.4, 12);
		tf1_W2_fullFit->SetNpx(280);
		tf1_W2_fullFit->SetParName(0, "W2 Gaus Norm (No Rej)" );
		tf1_W2_fullFit->SetParName(1, "W2 Gaus Mean (No Rej)" );
		tf1_W2_fullFit->SetParName(2, "W2 Gaus Sigma (No Rej)");
		tf1_W2_fullFit->SetParName(3, "W2 BG Pol8 p0 (No Rej)");
		tf1_W2_fullFit->SetParName(4, "W2 BG Pol8 p1 (No Rej)");
		tf1_W2_fullFit->SetParName(5, "W2 BG Pol8 p2 (No Rej)");
		tf1_W2_fullFit->SetParName(6, "W2 BG Pol8 p3 (No Rej)");
		tf1_W2_fullFit->SetParName(7, "W2 BG Pol8 p4 (No Rej)");
		tf1_W2_fullFit->SetParName(8, "W2 BG Pol8 p5 (No Rej)");
		tf1_W2_fullFit->SetParName(9, "W2 BG Pol8 p6 (No Rej)");
		tf1_W2_fullFit->SetParName(10, "W2 BG Pol8 p7 (No Rej)");
		tf1_W2_fullFit->SetParName(11, "W2 BG Pol8 p8 (No Rej)");

		tf1_W2_fullFit->SetParLimits(0, 4000, 6000);
		tf1_W2_fullFit->SetParLimits(1, 0.8, 0.9);
		tf1_W2_fullFit->SetParLimits(2, 0, 0.1);

		h_W2->Fit("tf1_W2_fullFit", "R0+");
		tf1_W2_fullFit->GetParameters(&par[0]);

//---------------------------------------------
//----------------NO REJECT----------------------
//---------------------------------------------

		tf1_W2_sig = new TF1("tf1_W2_sig", fit_gaus, 0, 1.4, 3);
		tf1_W2_sig->SetNpx(280);
		tf1_W2_sig->SetParameter(0, par[0]);
		tf1_W2_sig->SetParameter(1, par[1]);
		tf1_W2_sig->SetParameter(2, par[2]);


		W2_yield = (tf1_W2_sig->GetParameter(0))*(tf1_W2_sig->GetParameter(2))*sqrt(2*pi)/h_W2->GetBinWidth(1);

		h_tf_W2_sig = (TH1*)tf1_W2_sig->GetHistogram();
		h_tf_W2_sig->SetLineColor(6);
		h_tf_W2_sig->Scale(1./h_tf_W2_sig->Integral(), "width");
		h_tf_W2_sig->Scale(1./h_tf_W2_sig->Integral(), "height");
		h_tf_W2_sig->Scale(par[0]/(h_tf_W2_sig->GetMaximum()));
		// h_tf_W2_sig->Draw("hist+same");

	//------------------------------------------------------------
		tf1_W2_bg_pol8 = new TF1("tf1_W2_bg_pol8", fitPol8, 0, 1.4, 9);
		tf1_W2_bg_pol8->SetNpx(280);
		tf1_W2_bg_pol8->FixParameter(0, par[3]);
		tf1_W2_bg_pol8->FixParameter(1, par[4]);
		tf1_W2_bg_pol8->FixParameter(2, par[5]);
		tf1_W2_bg_pol8->FixParameter(3, par[6]);
		tf1_W2_bg_pol8->FixParameter(4, par[7]);
		tf1_W2_bg_pol8->FixParameter(5, par[8]);
		tf1_W2_bg_pol8->FixParameter(6, par[9]);
		tf1_W2_bg_pol8->FixParameter(7, par[10]);
		tf1_W2_bg_pol8->FixParameter(8, par[11]);

		tf1_W2_bg_pol8->SetLineColor(9);
		h_W2->Fit("tf1_W2_bg_pol8", "R0+");

		tf_sig_integral = h_tf_W2_sig->Integral();
		int W2_minus_bg_cnt = h_W2->GetEntries() - int(tf_sig_integral);
		TH1D *h_bg_rand = new TH1D("h_bg_rand", "Randomly filled BG from tf for BG", 280, 0, 1.4);
		cout << "Filling bg_rand with: " << W2_minus_bg_cnt << " events." << endl;
		h_bg_rand->FillRandom("tf1_W2_bg_pol8", W2_minus_bg_cnt);

		h_tf_W2_bg = (TH1*)tf1_W2_bg_pol8->GetHistogram();
		h_tf_W2_bg->SetLineColor(9);
		h_tf_W2_bg->Scale(1.0/h_tf_W2_bg->Integral(), "width");
		h_tf_W2_bg->Scale(1.0/h_tf_W2_bg->Integral(), "height");
		h_tf_W2_bg->Scale( (h_bg_rand->GetMaximum())/(h_tf_W2_bg->GetMaximum()) );

		TH1D *h_W2_sig = (TH1D*)h_W2->Clone("h_W2_sig");
		h_W2_sig->SetLineColor(6);
		for( int bin = 0; bin <= 280; bin++ ){
			double x = bin*0.005;
			double min_x = par[1] - 4*par[2];
			double max_x = par[1] + 4*par[2];

			if( x >= min_x && x <= max_x ){
				h_W2_sig->SetBinContent( bin, (h_W2->GetBinContent(bin)) - (h_tf_W2_bg->GetBinContent(bin)) );	
			}
			else{
				h_W2_sig->SetBinContent( bin, 0 );
			}			
		}

		h_W2_sig->Draw("hist+same");
		sig_integral = h_W2_sig->Integral();

		TLegend *tl_no_reject = new TLegend(0.15, 0.70, 0.45, 0.85,"Full Range BG Fit");
		tl_no_reject->AddEntry(tf1_W2_fullFit, "Total W2 fit");
		tl_no_reject->AddEntry(h_W2_sig, "W2 peak after BG sub");
		tl_no_reject->AddEntry(h_tf_W2_bg, "BG from full range (no rej.)");
		// tl_no_reject->AddEntry((TObject*)0, Form("Fit Chi-Sq = %.2f", tf1_W2_fullFit->GetChisquare()), "");
		tl_no_reject->Draw("same");

		TPaveText *tp_W2 = new TPaveText(0.15, 0.64, 0.45, 0.69, "NDCbr");
		tp_W2->AddText(Form("Estimated real W2 count: %i", int(sig_integral)));
		tp_W2->Draw("same");

		cout << "-------------------------------------" << endl;
		cout << "tf_sig_integral: " << tf_sig_integral << ", sig_integral: " << sig_integral << endl;
		cout << "-------------------------------------" << endl;

		//RESIDUAL

		h_tf_W2_full = (TH1*)tf1_W2_fullFit->GetHistogram();
		h_tf_W2_full->SetLineColor(7);
		h_tf_W2_full->Scale(1.0/h_tf_W2_full->Integral(), "width");
		h_tf_W2_full->Scale(1.0/h_tf_W2_full->Integral(), "height");
		h_tf_W2_full->Scale( (h_W2->GetMaximum())/(h_tf_W2_full->GetMaximum()) );

		h_W2_full_rand->FillRandom("tf1_W2_fullFit", h_W2->GetEntries());
		h_W2_full_rand->SetLineColor(2);

		for(int bin = 0; bin <= 280; bin++){
			h_W2_residual->SetBinContent( bin, ( h_W2->GetBinContent(bin) - h_W2_full_rand->GetBinContent(bin) ) );
		}


//---------------------------------------------
//-----------------REJECTING--------------------
//---------------------------------------------
		if( fit_with_reject ){
			h_W2copy = (TH1D*)h_W2->Clone("h_W2copy");
			if( true ){
				TCanvas *c_W2_fullFit_reject = new TCanvas("c_W2_fullFit_reject", "c_W2_fullFit_reject", 600, 500);
				h_W2copy->Draw();

		//----- For the reject BG fit we first need to find the BG fit with the reject in place
		//----- Then we take the parameters from that and fix the Full Fit BG parametes to match

				if( polN == 4 ){
					tf1_W2_bg_rej_init = new TF1("tf1_W2_bg_rej_init", fitPol4_with_reject, 0, 1.4, polN+1);
				}

				if( polN == 6 ){
					tf1_W2_bg_rej_init = new TF1("tf1_W2_bg_rej_init", fitPol6_with_reject, 0, 1.4, polN+1);
				}

				if( polN == 8 ){
					tf1_W2_bg_rej_init = new TF1("tf1_W2_bg_rej_init", fitPol8_with_reject, 0, 1.4, polN+1);
				}
				tf1_W2_bg_rej_init->SetNpx(280);
				h_W2copy->Fit("tf1_W2_bg_rej_init", "R0+");

				tf1_W2_bg_rej_init->GetParameters(&parRej[0]);

				if( polN == 4 ){
					tf1_W2_bg_rej = new TF1("tf1_W2_bg_rej", fitPol4, 0, 1.4, polN+1);					
				}

				if( polN == 6 ){
					tf1_W2_bg_rej = new TF1("tf1_W2_bg_rej", fitPol6, 0, 1.4, polN+1);					
				}

				if( polN == 8 ){
					tf1_W2_bg_rej = new TF1("tf1_W2_bg_rej", fitPol8, 0, 1.4, polN+1);					
				}

				tf1_W2_bg_rej->SetNpx(280);

				for( int param = 0; param <= polN; param++ ){
					tf1_W2_bg_rej->FixParameter(param, parRej[param]);					
				}
				tf1_W2_bg_rej->SetLineColor(3);
				// tf1_W2_bg_pol8_rej->Draw("hist+same");

				//Need a new pol8 without the rejection to draw across reject zone:

				if( polN == 4 ){
					tf1_W2_bg_rej_full = new TF1("tf1_W2_bg_rej_full", fitPol4, 0, 1.4, polN+1);					
				}

				if( polN == 6 ){
					tf1_W2_bg_rej_full = new TF1("tf1_W2_bg_rej_full", fitPol6, 0, 1.4, polN+1);					
				}

				if( polN == 8 ){
					tf1_W2_bg_rej_full = new TF1("tf1_W2_bg_rej_full", fitPol8, 0, 1.4, polN+1);					
				}

				tf1_W2_bg_rej_full->SetNpx(280);

				for( int param = 0; param <= polN; param++ ){
					tf1_W2_bg_rej_full->FixParameter(param, parRej[param]);					
				}	

				tf1_W2_bg_rej_full->SetLineColor(9);	
				h_W2copy->Fit("tf1_W2_bg_rej_full", "R0+");
		//----------
			// Using regular fullFitFucntion because we got the BG with Reject just above
				tf1_W2_fullFit_rej = new TF1("tf1_W2_fullFit_rej", fullFitFunction, 0, 1.4, polN+4);
				tf1_W2_fullFit_rej->SetNpx(280);
				tf1_W2_fullFit_rej->SetLineColor(2);

				tf1_W2_fullFit_rej->SetParName(0, "W2 Gaus Norm (Rej)" );
				tf1_W2_fullFit_rej->SetParName(1, "W2 Gaus Mean (Rej)" );
				tf1_W2_fullFit_rej->SetParName(2, "W2 Gaus Sigma (Rej)");

				for( int param = 3; param < polN+4; param++ ){
					tf1_W2_fullFit_rej->SetParName(param, Form("W2 BG Pol%i p%i (Rej)", polN, param-3));
				}

			//Fix parameters from the BG_reject fit
				for(int param = 3; param < polN+4; param++ ){
					tf1_W2_fullFit_rej->FixParameter(param, parRej[param-3]);
				}

				tf1_W2_fullFit_rej->SetParLimits(0, 4000, 6000);
				tf1_W2_fullFit_rej->SetParLimits(1, 0.8, 0.9);
				tf1_W2_fullFit_rej->SetParLimits(2, 0, 0.1);

				h_W2copy->Fit("tf1_W2_fullFit_rej", "R0+");
				tf1_W2_fullFit_rej->GetParameters(&par[0]);
			}

			tf1_W2_sig_rej = new TF1("tf1_W2_sig_rej", fit_gaus, 0, 1.4, 3);
			tf1_W2_sig_rej->SetNpx(280);
			tf1_W2_sig_rej->SetParameter(0, par[0]);
			tf1_W2_sig_rej->SetParameter(1, par[1]);
			tf1_W2_sig_rej->SetParameter(2, par[2]);

			W2_yield_reject = (tf1_W2_sig_rej->GetParameter(0))*(tf1_W2_sig_rej->GetParameter(2))*sqrt(2*pi)/h_W2copy->GetBinWidth(1);

			h_tf_W2_sig_rej = (TH1*)tf1_W2_sig_rej->GetHistogram();
			h_tf_W2_sig_rej->SetLineColor(6);
			h_tf_W2_sig_rej->Scale(1./h_tf_W2_sig_rej->Integral(), "width");
			h_tf_W2_sig_rej->Scale(1./h_tf_W2_sig_rej->Integral(), "height");
			h_tf_W2_sig_rej->Scale(par[0]/(h_tf_W2_sig_rej->GetMaximum()));

		//------------------------------------------------------------


			tf_sig_integral_rej = h_tf_W2_sig_rej->Integral();
			int W2_minus_bg_cnt_rej = h_W2copy->GetEntries() - int(tf_sig_integral_rej);
			

			TH1D *h_bg_rand_rej = new TH1D("h_bg_rand_rej", "Randomly filled BG from tf for BG with Reject", 280, 0, 1.4);
			cout << "Filling Reject bg_rand with: " << W2_minus_bg_cnt_rej << " events." << endl;
			h_bg_rand_rej->FillRandom("tf1_W2_bg_pol8_rej_full", W2_minus_bg_cnt_rej);

			h_tf_W2_bg_rej = (TH1*)tf1_W2_bg_rej_full->GetHistogram();
			h_tf_W2_bg_rej->SetLineColor(9);
			h_tf_W2_bg_rej->Scale(1.0/h_tf_W2_bg_rej->Integral(), "width");
			h_tf_W2_bg_rej->Scale(1.0/h_tf_W2_bg_rej->Integral(), "height");
			h_tf_W2_bg_rej->Scale( (h_bg_rand_rej->GetMaximum())/(h_tf_W2_bg_rej->GetMaximum()) );

			TH1D *h_W2_sig_rej = (TH1D*)h_W2copy->Clone("h_W2_sig_rej");
			h_W2_sig_rej->SetLineColor(6);
			for( int bin = 0; bin <= 280; bin++ ){

				TAxis *W2_xaxis = h_W2->GetXaxis();
				// double x = bin*0.005;
				double min_x = par[1] - 4*par[2];
				double max_x = par[1] + 4*par[2];

				Int_t min_bin_x = W2_xaxis->FindBin(min_x);
				Int_t max_bin_x = W2_xaxis->FindBin(max_x);

				if( bin >= min_bin_x && bin <= max_bin_x ){
					h_W2_sig_rej->SetBinContent( bin, (h_W2copy->GetBinContent(bin)) - (h_tf_W2_bg_rej->GetBinContent(bin)) );	
				}
				else{
					h_W2_sig_rej->SetBinContent( bin, 0 );
				}			
			}

			h_W2_sig_rej->Draw("hist+same");
			sig_integral_rej = h_W2_sig_rej->Integral();

			TLegend *tl_reject = new TLegend(0.15, 0.70, 0.45, 0.85, "BG with Reject Region");
			tl_reject->AddEntry(tf1_W2_fullFit_rej, "Total W2 fit");
			tl_reject->AddEntry(h_W2_sig_rej, "W2 peak after BG sub");
			tl_reject->AddEntry(h_tf_W2_bg_rej, "BG from full range (with rej.)");
			// tl_reject->AddEntry((TObject*)0, Form("Fit Chi-Sq = %.2f", tf1_W2_fullFit_rej->GetChisquare()), "");
			tl_reject->Draw("same");		

			TPaveText *tp_W2_rej = new TPaveText(0.15, 0.64, 0.45, 0.69, "NDCbr");
			tp_W2_rej->AddText(Form("Estimated real W2 count: %i", int(sig_integral_rej)));
			tp_W2_rej->Draw("same");

			cout << "-------------------------------------" << endl;
			cout << " WITH REJECT POINTS " << endl;	
			cout << "-------------------------------------" << endl;
			cout << "tf_sig_integral_rej: " << tf_sig_integral_rej << ", sig_integral_rej: " << sig_integral_rej << endl;
			cout << "-------------------------------------" << endl;
		}

		cout << " NO REJECT " << endl;
		cout << "-------------------------------------" << endl;
		cout << "tf_sig_integral: " << tf_sig_integral << ", sig_integral: " << sig_integral << endl;

//------------------------------------------------------
//---------------------- DX STUFF ----------------------
//------------------------------------------------------

		cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
		cout << "        Working on dx plots   " << endl;
		TCanvas *c_dx = new TCanvas("c_dx", "c_dx", 600, 500);
		h_dx->Draw();

		tf1_dx_fullFit = new TF1("tf1_dx_fullFit", fullFitFunction, -2.5, 2.5, 12);
		tf1_dx_fullFit->SetNpx(1000);
		tf1_dx_fullFit->SetParName(0, "dx Gaus Norm (No Rej)" );
		tf1_dx_fullFit->SetParName(1, "dx Gaus Mean (No Rej)" );
		tf1_dx_fullFit->SetParName(2, "dx Gaus Sigma (No Rej)");
		tf1_dx_fullFit->SetParName(3, "dx BG Pol8 p0 (No Rej)");
		tf1_dx_fullFit->SetParName(4, "dx BG Pol8 p1 (No Rej)");
		tf1_dx_fullFit->SetParName(5, "dx BG Pol8 p2 (No Rej)");
		tf1_dx_fullFit->SetParName(6, "dx BG Pol8 p3 (No Rej)");
		tf1_dx_fullFit->SetParName(7, "dx BG Pol8 p4 (No Rej)");
		tf1_dx_fullFit->SetParName(8, "dx BG Pol8 p5 (No Rej)");
		tf1_dx_fullFit->SetParName(9, "dx BG Pol8 p6 (No Rej)");
		tf1_dx_fullFit->SetParName(10, "dx BG Pol8 p7 (No Rej)");
		tf1_dx_fullFit->SetParName(11, "dx BG Pol8 p8 (No Rej)");

		tf1_dx_fullFit->SetParLimits(0, 0.5*h_dx->GetMaximum(), h_dx->GetMaximum());
		tf1_dx_fullFit->SetParLimits(1, -0.05, 0.05);
		tf1_dx_fullFit->SetParLimits(2, 0, 0.07);

		h_dx->Fit("tf1_dx_fullFit", "R0+");
		tf1_dx_fullFit->GetParameters(&par[0]);		

		tf1_dx_sig = new TF1("tf1_dx_sig", fit_gaus, -2.5, 2.5, 3);
		tf1_dx_sig->SetNpx(1000);
		tf1_dx_sig->SetParameter(0, par[0]);
		tf1_dx_sig->SetParameter(1, par[1]);
		tf1_dx_sig->SetParameter(2, par[2]);

		h_tf_dx_sig = (TH1*)tf1_dx_sig->GetHistogram();
		h_tf_dx_sig->SetLineColor(6);
		h_tf_dx_sig->Scale(1./h_tf_dx_sig->Integral(), "width");
		h_tf_dx_sig->Scale(1./h_tf_dx_sig->Integral(), "height");
		h_tf_dx_sig->Scale(par[0]/(h_tf_dx_sig->GetMaximum()));

	//------------------------------------------------------------
		tf1_dx_bg_pol8 = new TF1("tf1_dx_bg_pol8", fitPol8, -2.5, 2.5, 9);
		tf1_dx_bg_pol8->SetNpx(1000);
		tf1_dx_bg_pol8->FixParameter(0, par[3]);
		tf1_dx_bg_pol8->FixParameter(1, par[4]);
		tf1_dx_bg_pol8->FixParameter(2, par[5]);
		tf1_dx_bg_pol8->FixParameter(3, par[6]);
		tf1_dx_bg_pol8->FixParameter(4, par[7]);
		tf1_dx_bg_pol8->FixParameter(5, par[8]);
		tf1_dx_bg_pol8->FixParameter(6, par[9]);
		tf1_dx_bg_pol8->FixParameter(7, par[10]);
		tf1_dx_bg_pol8->FixParameter(8, par[11]);

		tf1_dx_bg_pol8->SetLineColor(9);
		h_dx->Fit("tf1_dx_bg_pol8", "R0+");

		tf_dx_sig_integral = h_tf_dx_sig->Integral();
		int dx_minus_bg_cnt = h_dx->GetEntries() - int(tf_dx_sig_integral);
		TH1D *h_dx_bg_rand = new TH1D("h_dx_bg_rand", "dx: Randomly filled BG from tf for BG", 1000, -2.5, 2.5);
		cout << "Filling dx bg_rand with: " << dx_minus_bg_cnt << " events." << endl;
		h_dx_bg_rand->FillRandom("tf1_dx_bg_pol8", dx_minus_bg_cnt);

		h_dx->GetXaxis()->SetRangeUser(-2.5, -0.3);
		// dx_bg_max = h_dx->GetMaximum();
		dx_bg_max = tf1_dx_bg_pol8->GetMaximum();
		h_dx->GetXaxis()->SetRangeUser(-2.5, 2.5);

		h_tf_dx_bg = (TH1*)tf1_dx_bg_pol8->GetHistogram();
		h_tf_dx_bg->SetLineColor(7);
		h_tf_dx_bg->Scale(1.0/h_tf_dx_bg->Integral(), "width");
		h_tf_dx_bg->Scale(1.0/h_tf_dx_bg->Integral(), "height");
		h_tf_dx_bg->Scale( (dx_bg_max)/(h_tf_dx_bg->GetMaximum()) );

		TH1D *h_dx_sig = (TH1D*)h_dx->Clone("h_dx_sig");
		h_dx_sig->SetLineColor(6);
		for( int bin = 0; bin <= 1000; bin++ ){
			
			TAxis *dx_xaxis = h_dx->GetXaxis();
			// double x = bin*0.005;
			double min_x = par[1] - 4*par[2];
			double max_x = par[1] + 4*par[2];

			Int_t min_bin_x = dx_xaxis->FindBin(min_x);
			Int_t max_bin_x = dx_xaxis->FindBin(max_x);

			if( bin >= min_bin_x && bin <= max_bin_x ){
				h_dx_sig->SetBinContent( bin, (h_dx->GetBinContent(bin)) - (h_tf_dx_bg->GetBinContent(bin)) );	
			}
			else{
				h_dx_sig->SetBinContent( bin, 0 );
			}			
		}

		h_dx_sig->Draw("hist+same");
		sig_integral = h_dx_sig->Integral();

		TLegend *tl_dx = new TLegend(0.15, 0.70, 0.45, 0.85,"Full Range BG Fit");
		tl_dx->AddEntry(tf1_dx_fullFit, "Total dx fit");
		tl_dx->AddEntry(h_dx_sig, "dx peak after BG sub");
		tl_dx->AddEntry(h_tf_dx_bg, "BG from full range (no rej.)");
		// tl_no_reject->AddEntry((TObject*)0, Form("Fit Chi-Sq = %.2f", tf1_dx_fullFit->GetChisquare()), "");
		tl_dx->Draw("same");

		TPaveText *tp_dx = new TPaveText(0.15, 0.64, 0.45, 0.69, "NDCbr");
		tp_dx->AddText(Form("Estimated real dx count: %i", int(sig_integral)));
		tp_dx->Draw("same");

	}
		auto total_time_end = high_resolution_clock::now();
		auto total_time_duration = duration_cast<minutes>(total_time_end - total_time_start);
		cout << "Total time for analysis: " << total_time_duration.count() << " minutes. " << endl;
		cout << endl << "Outfile: " << outfile_name.Data() << endl;
}