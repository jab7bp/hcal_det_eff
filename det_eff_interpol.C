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

TFitResultPtr r;

//Calc Method for kinematic-based variables:
//1) Use four-momentum functions
//2) Use tree kinematic variables (e.kine.**)
//3) Use reconstructed variable functions

int calc_method = 0; //1 isn't good

bool fit_with_expo = false;

int polN = 4;
int interpolN = 99;
int dypolN = 8;
bool use_parsed = false;

bool use_bbcal_cuts = false;
bool fiducial_cut = true;
bool correct_beam_energy = false;
bool apply_fcut = false;

bool fit_only = true;
bool calc_W_only = false;
bool fit_with_reject = true;

bool acceptance_match = true;

bool plot_dxdy = true;
bool plot_dxdy_anticut = true;

TF1 *tf1_W2_interpolate, *tf1_W2_elastics, *tf1_W2_elastics_init, *tf1_W2_interpol_bg, *tf1_W2_interpol_sig;
TF1 *tf1_W2_interpolate_full;
TF1 *tf1_W2_bgCut_pol8, *tf1_W2_bgCut_sig, *tf1_W2_bgCut_full, *tf1_W2_bgCut_gaus;
TF1 *tf1_W2_anticut_init, *tf1_W2_anticut;
TH1 *h_tf_W2_bgCut_bg;
TH1D *h_W2_bgCut, *h_W2_bgCut_sig, *h_W2_bgCut_bg, *h_W2_bgCut_bg_sub;
double dx_sig_integral;
double bgCut_sig_cnt, bgCut_sig_gaus_cnt;
double bgCut_full_integral, bgCut_pol8_integral, sigma_bgCut_pol8_integral;
double bgCut_sig_integral, sigma_bgCut_full_integral, sigma_bgCut_sig_integral, sigma_bgCut_TOTAL;
double dx_full_integral, sigma_dx_full_integral, sigma_dx_sig_integral, dx_bg_integral, sigma_dx_bg_integral, sigma_dx_TOTAL;
double det_eff_FINAL, sigma_det_eff_FINAL;
double bgCut_gaus_min, bgCut_gaus_max, dx_sig_gaus_min, dx_sig_gaus_max;

TH1D *h_dx_subrange_sig;

TFitResultPtr tfrp_W2_bgCut_gaus;
TF1 *tf1_W2_bgCut_gaus_err;
TMatrixDSym cov_W2_bgCut_gaus;

Double_t W2_bgCut_chi2, W2_bgCut_gaus_par0, W2_bgCut_gaus_par1, W2_bgCut_gaus_par2, W2_bgCut_gaus_parErr0, W2_bgCut_gaus_parErr1, W2_bgCut_gaus_parErr2;

TH1D *h_W2elastics, *h_W2elastics_thr2, *h_W2elastics_fcut, *h_W2elastics_init, *h_W2elastics_gaus, *h_W2elastics_corr;
TH1D *h_dy_sig, *h_W2_sig_rej, *h_W2_sig;
TH1 *h_tf_W2_interpol_bg;
TH1D *h_W2_interpol_bg, *h_W2_interpol_sig;
double parInterpol[14];
double dy_sig_mean, dy_sig_sigma, dy_sig_min, dy_sig_max;
double dx_sig_mean, dx_sig_sigma;
double reject_min, reject_max;
vector<vector<double>> iter_xmin_xmax_chi2(1000, vector<double>(3));

#include "/w/halla-scshelf2102/sbs/jboyd/include/fit_functions.h"

Double_t W2FullInterpol(Double_t *x, Double_t *par){ //get poly fit to bg with scaled fit to "pure elastics"
  Double_t W_2 = x[0];
  Double_t sig_scale = par[0];
  Double_t signal = sig_scale * h_W2elastics->Interpolate(W_2);

  if( interpolN == 99 ){
  	return signal + fit_expo(x, &par[1]);  	
  }
  if( interpolN == 3 ){
  	return signal + fitPol3(x, &par[1]);  	
  }
  if( interpolN == 4 ){
  	return signal + fitPol4(x, &par[1]);  	
  }
   if( interpolN == 5 ){
  	return signal + fitPol5(x, &par[1]);  	
  }
  if( interpolN == 6 ){
  	return signal + fitPol6(x, &par[1]);  	
  }
  if( interpolN == 8 ){
  	return signal + fitPol8(x, &par[1]);  	
  }

}

double gaus_integral(double par0, double par2, double binWidth){
	return (par0)*(par2)*sqrt(2*TMath::Pi())/(binWidth);
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
double dx_scale = 3.0;
double dy_scale = 2.0;
double W2_scale = 2.0;
double dx_min, dx_max, dy_min, dy_max, W2_min, W2_max;
bool global_cut_fail = false;
bool hit_on_HCal = false;

double x_min_W2 = 0.0;
double x_max_W2 = 1.4;
int nBins_W2 = int((x_max_W2 - x_min_W2)*200.0);

double W2_fit_min = 0.0;
double W2_fit_max = 1.4;
int nBins_fit_W2 = int((W2_fit_max - W2_fit_min)*200.0);

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

bool theta_pq_cut = true;
double theta_pq_p_thresh;
double theta_pq_n_thresh;

//Experimental Constants, Thresholds, cuts, etc DEFINITIONS
const double pi = TMath::Pi();
const double Mp = 0.938272; //Mass of proton [GeV]
const double Mn = 0.939565; //Mass of neutron [GeV]
double Nucleon_mass;
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
double p_el, nu, p_Nucleon, nucleon_theta, nucleon_phi, E_ep, p_ep, Q2, W, W2, E_pp, E_nucleon, KE_p, dx, dy;

TLorentzVector pN;
TVector3 pNhat;

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

double p_recon, nu_recon, E_loss, E_corr, theta_pq_n, theta_pq_p, theta_pq_n_pNhat, theta_pq_p_pNhat, thetapq_n, thetapq_p;

int useAlshield = 0;

//Declare vars
Double_t atime[kNcell], row[kNcell], col[kNcell], tdctime[kNcell], cblkid[kNcell], cblke[kNcell];
Double_t nblk, nclus, SH_nclus, PS_nclus, hcal_x, hcal_y, hcal_e;

Double_t par[14], parRej[14], parRejSave[14], parDY[14], parAntiCut[14];

double bb_tr_p[maxTracks], bb_tr_px[maxTracks], bb_tr_py[maxTracks], bb_tr_pz[maxTracks];
double bb_tr_vx[maxTracks], bb_tr_vy[maxTracks], bb_tr_vz[maxTracks], bb_tr_chi2[maxTracks];
double bb_fp_x[maxTracks], bb_fp_y[maxTracks], bb_fp_th[maxTracks], bb_fp_ph[maxTracks];
double bb_tgt_x[maxTracks], bb_tgt_y[maxTracks], bb_tgt_th[maxTracks], bb_tgt_ph[maxTracks];
Double_t e_kine_W2, e_kine_Q2, e_kine_theta_eprime, e_kine_nu;
Double_t e_kine_omega, e_kine_ph_q, e_kine_theta_Q;
double hcal_clusblk_ADC_time[15]; //Maximum number of blocks in a cluster is 15 as per S. Seeds
double bb_tr_n, bb_ps_x, bb_ps_y, bb_ps_e, bb_sh_x, bb_sh_y, bb_sh_e;
double bb_ntracks[maxTracks], bb_nhits[maxTracks];

Double_t TDCT_id[kNtdc], TDCT_tdc[kNtdc], hodo_tmean[kNtdc]; 
Int_t TDCTndata;

Long64_t nTCevents, Nevents;

//INITIALIZE ALL HISTOGRAMS:
TH1D *hin_bb_gem_Ntracks, *h_W2_full_rand;
TH1D *h_atime, *h_W, *h_W2, *h_W2copy, *h_W2interpol, *h_W2recon, *h_KE_p, *h_KE_low, *h_Diff, *h_X, *h_Y, *h_E_eloss, *h_hcal_clusblk_ADC_time;
TH1D *h_Wfull, *h_W2full, *h_W2_residual;
TH1D *h_W_cut, *h_W_fcut, *h_vz_cut, *h_theta_pq_p, *h_theta_pq_n, *h_theta_pq_p_pNhat, *h_theta_pq_n_pNhat;
TH1D *h_Q2, *h_E_ep, *h_E_pp;
TH1D *h_dy, *h_dy_cut, *h_dy_wcut, *h_dx, *h_dx_cut, *h_dx_wcut, *h_dx_fcut, *h_dx_wcut_fcut, *h_dy_wcut_fcut;
TH1D *h_dx_w2cut, *h_dx_w2cut_fcut;
TH1D *h_dy_subrange, *h_dx_subrange;
TH1D *h_Nevents;

TH1D *h_rej_bg, *h_rej_sig_gaus;

//Interpolation "tight" histograms
TH1D *h_W2antiCut, *h_W2fullCut, *h_W2elasticsFull, *h_W2antiCutFull, *h_dx_allCuts;
TH1D *h_W2antiCut_fcut, *h_W2fullCut_fcut, *h_W2elasticsFull_fcut, *h_W2antiCutFull_fcut, *h_dx_allCuts_fcut;

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

TF1 *tf1_W2_bg_pol8, *tf1_W2_sig, *tf1_W2_fullFit, *tf1_W2_bg_rej, *tf1_W2_bg_rej_full, *tf1_W2_bg_rej_init, *tf1_W2_bg_pol8_rej_full, *tf1_W2_bg_pol8_rej, *tf1_W2_sig_rej, *tf1_W2_fullFit_rej;
TF1 *tf1_dy_bg_rej, *tf1_dy_bg_rej_full, *tf1_dy_bg_rej_init, *tf1_dy_sig_rej;
TF1 *tf1_dx_fullFit, *tf1_dx_bg_pol8, *tf1_dx_sig, *tf1_dx_bg, *tf1_W2_bg;
TF1 *tf1_dy_fullFit, *tf1_dy_bg_full, *tf1_dy_bg, *tf1_dy_sig;
TH1 *h_tf_W2_bg, *h_tf_W2_sig, *h_tf_W2_bg_rej, *h_tf_W2_sig_rej, *h_tf_W2_full;
TH1 *h_tf_dx_bg, *h_tf_dx_sig, *h_tf_dx_full;
TH1 *h_tf_dy_bg, *h_tf_dy_sig, *h_tf_dy_bg_rej, *h_tf_dy_sig_rej, *h_tf_dy_full;
TH1* h_tf_W2elastics_init, *h_tf_W2antiCut_init;
TH1D *h_dx_sig, *h_dx_bg;

int polNfit, interpolNfit;
int min_max_iter_cnt;

void det_eff_interpol(){

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
	cout << "HCal angle [deg]: " << (180/pi)*lookup_HCal_angle_by_kine(kine, "deg") << endl;
	cout << "HCal distance: " << HCal_dist << endl;
	cout << "-----------------------------------" << endl << endl;

	if( run_target == "LH2" ){
		Nucleon_mass = Mp;
	}
	if( run_target == "LD2" ){
		Nucleon_mass = 0.5*(Mp + Mn);
	}

	if( polN == 99 ){ //99 is exponential
		polN = 3;
	}

	if( kine == 4 ){
		pass = 0;
		hcalheight = -0.312479; // Height of HCal above beamline

		W_mean = 0.918;
		W_sigma = 0.053;
		W2_mean = 0.835;
		W2_sigma = 0.1;

		//offset: -1.37521e-01
		dx_min = -0.0061781621 - (dx_scale)*0.057494823;
		dx_max = -0.0061781621 + (dx_scale)*0.057494823;
		dx_p = -0.0061781621;
		dx_p_sigma = 0.057494823;
		dx_n = dx_p;
		dx_n_sigma = dx_n_sigma;

		dy_min = -0.0374881 - (dy_scale)*0.056;
		dy_min = -0.0374881 + (dy_scale)*0.056;

		dy_p = -0.0374881;
		dy_p_sigma = 0.056;
		dy_n = dy_p;
		dy_n_sigma = dy_p_sigma;

		dx_pn_max = 0;

		W2_min = W2_mean - (W2_scale)*W2_sigma;
		W2_max = W2_mean + (W2_scale)*W2_sigma;
	}
	if( kine == 8 ){
		pass = 1;
		hcalheight = -0.450; // Height of HCal above beamline
	}
	
	outfile_name = Form("rootfiles/det_eff_histos_SBS%i_%s_mag%i.root", kine, run_target.Data(), sbsfieldscale);

	cout << "Outfile name: " << outfile_name.Data() << endl;

	if( !fit_only ){

		BB_dist = lookup_BB_dist_by_kine(kine);
		BB_theta = lookup_BB_angle_by_kine(kine, "rad");
		HCal_dist = lookup_HCal_dist_by_kine(kine);
		HCal_theta = lookup_HCal_angle_by_kine(kine, "rad");
		// W_mean = lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "W_mean");
		// W_sigma = lookup_parsed_cut(run_target.Data(), kine, sbsfieldscale, "W_sigma");

		hin_bb_gem_Ntracks = new TH1D("hin_bb_gem_Ntracks", "Number of tracks found", 11, -0.5, 10.5);
		ADC_time_min = lookup_ADC_time_cut(run_target, kine, sbsfieldscale, "ADC_time_min");
		ADC_time_max = lookup_ADC_time_cut(run_target, kine, sbsfieldscale, "ADC_time_max");

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

		//Kinematic variables
		TC->SetBranchStatus( "e.kine.W2", 1 );
		TC->SetBranchStatus( "e.kine.Q2", 1 );
		TC->SetBranchStatus( "e.kine.angle", 1 );
		TC->SetBranchStatus( "e.kine.nu", 1 );
		TC->SetBranchStatus( "e.kine.omega", 1 );
		TC->SetBranchStatus( "e.kine.ph_q", 1 );
		TC->SetBranchStatus( "e.kine.th_q", 1 );

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

		//Kinematic variables
		TC->SetBranchAddress( "e.kine.W2", &e_kine_W2 );
		TC->SetBranchAddress( "e.kine.Q2", &e_kine_Q2 );
		TC->SetBranchAddress( "e.kine.angle", &e_kine_theta_eprime );
		TC->SetBranchAddress( "e.kine.nu", &e_kine_nu );
		TC->SetBranchAddress( "e.kine.omega", &e_kine_omega );
		TC->SetBranchAddress( "e.kine.ph_q", &e_kine_ph_q );
		TC->SetBranchAddress( "e.kine.th_q", &e_kine_theta_Q );

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
		TC->SetBranchAddress( "bb.gem.track.nhits", bb_nhits );

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
		h_W = new TH1D("h_W", Form("Invariant Mass W - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), nBins_W2, x_min_W2, x_max_W2);
		h_W_cut = new TH1D("h_W_cut", Form("Invariant Mass W (Coin & Vert Cuts) - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 600, 0.0, 3.0);
		h_W_fcut = new TH1D("h_W_fcut", Form("Invariant Mass W (Fiduc. Cuts) - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 600, 0.0, 3.0);
		h_W2 = new TH1D("h_W2", Form("Invariant Mass W^2 - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), nBins_W2, x_min_W2, x_max_W2);

		h_Wfull = new TH1D("h_Wfull", Form("Invariant Mass W Full Range- SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 600, 0.0, 3.0);
		h_W2full = new TH1D("h_W2full", Form("Invariant Mass W^2 Full Range - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 600, 0.0, 3.0);

	//---histograms for tight cuts
		h_W2elastics = new TH1D("h_W2elastics", Form("Tight Cut for Real Elastics of Invariant Mass W^2 - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), nBins_W2, x_min_W2, x_max_W2);
		h_W2elastics_thr2 = new TH1D("h_W2elastics_thr2", Form("Thresh 2 Tight Cut for Real Elastics of Invariant Mass W^2 - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), nBins_W2, x_min_W2, x_max_W2);


		h_W2elasticsFull = new TH1D("h_W2elasticsFull", Form("(Full-Range) Tight Cut for Real Elastics of Invariant Mass W^2 - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 600, 0.0, 3.0);

		h_W2antiCut = new TH1D("h_W2antiCut", Form("Anti-Cut for Real Elastics of Invariant Mass W^2 - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), nBins_W2, x_min_W2, x_max_W2);
		h_W2antiCutFull = new TH1D("h_W2antiCutFull", Form("(Full-Range) Anti-Cut for Real Elastics of Invariant Mass W^2 - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 600, 0.0, 3.0);
		h_W2fullCut = new TH1D("h_W2fullCut", Form("(Full-Range) Tight Cut for Real Elastics of Invariant Mass W^2 - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 600, 0.0, 3.0);
		
		h_dx_allCuts = new TH1D("h_dx_allCuts", Form("dx (All Tight Elastic Cuts) - SBS%i %i, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 1000, -2.5, 2.5);

		h_W2elastics_fcut = new TH1D("h_W2elastics_fcut", Form("Tight Cut for Real Elastics of Invariant Mass W^2 (with FCut) - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), nBins_W2, x_min_W2, x_max_W2);
		h_W2elasticsFull_fcut = new TH1D("h_W2elasticsFull_fcut", Form("(Full-Range) Tight Cut for Real Elastics of Invariant Mass W^2 (with FCut) - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 600, 0.0, 3.0);

		h_W2antiCut_fcut = new TH1D("h_W2antiCut_fcut", Form("Anti-Cut for Real Elastics of Invariant Mass W^2 (with FCut) - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), nBins_W2, x_min_W2, x_max_W2);
		h_W2antiCutFull_fcut = new TH1D("h_W2antiCutFull_fcut", Form("(Full-Range) Anti-Cut for Real Elastics of Invariant Mass W^2 (with FCut) - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 600, 0.0, 3.0);
		h_W2fullCut_fcut = new TH1D("h_W2fullCut_fcut", Form("(Full-Range) Tight Cut for Real Elastics of Invariant Mass W^2 (with FCut) - SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), 600, 0.0, 3.0);
		
		h_dx_allCuts_fcut = new TH1D("h_dx_allCuts_fcut", Form("dx (All Tight Elastic Cuts) (with FCut) - SBS%i %i, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 1000, -2.5, 2.5);

	//---

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

		h_dxdy = new TH2D("h_dxdy", Form("Hadron Spot(s) on HCal (NO CUTS) - SBS%i %i, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), 600, -1.5, 1.5, 1000, -2.5, 2.5 );
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

		h_theta_pq_n = new TH1D("h_theta_pq_n", Form("Theta pq for neutron - SBS%i %i, %s; theta_pq_n (rad);", kine, sbsfieldscale, run_target.Data()), 120, 0.0, 0.6);
		h_theta_pq_p = new TH1D("h_theta_pq_p", Form("Theta pq for proton - SBS%i %i, %s; theta_pq_p (rad);", kine, sbsfieldscale, run_target.Data()), 120, 0.0, 0.6);

		h_theta_pq_n_pNhat = new TH1D("h_theta_pq_n_pNhat", Form("Theta pq for neutron (pNhat) - SBS%i %i, %s; theta_pq_n (rad);", kine, sbsfieldscale, run_target.Data()), 120, 0.0, 0.6);
		h_theta_pq_p_pNhat = new TH1D("h_theta_pq_p_pNhat", Form("Theta pq for proton (pNhat) - SBS%i %i, %s; theta_pq_p (rad);", kine, sbsfieldscale, run_target.Data()), 120, 0.0, 0.6);


		cout << "--------------------------------------" << endl;
		cout << "Number of events to analyze: " << Nevents << endl;
		cout << "--------------------------------------" << endl;
		cout << "--------------------------------------" << endl;
		cout << "Starting analysis loop on events..... " << endl;

	//Basic Energy calcs

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

		//Correct the beam energy with energy loss in target using vertex position
      		Double_t Eloss = (bb_tr_vz[0]+l_tgt/2.0) * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //approximately 3 MeV
      		h_E_eloss->Fill( Eloss );

		    if( correct_beam_energy ){
	      		E_beam_final = E_beam - Eloss;
	      	}

	      	if( !correct_beam_energy){
	      		E_beam_final = E_beam;
	      	}

	      	h_E_ecorr_vs_vert->Fill( bb_tr_vz[0], E_beam_final);
	    ////Corrections
		    Double_t E_corr = E_beam - Eloss;	

	      	p_corr = bb_tr_p[0] - E_loss_outgoing; //Neglecting mass of e'


	    //Proceed only if at least one track exists in BB arm - lowest chi2 track always first element
	      	if( bb_tr_n > 1){
	      		continue;
	      	}

	//-------------
 	
	      	
	      	// e_prime_theta = acos( bb_tr_pz[0]/bb_tr_p[0] ); //Uncorrected track momenutm to reconstruct e' theta
	      	// e_prime_phi = atan2( bb_tr_py[0], bb_tr_px[0]);

	      	TVector3 vertex( 0, 0, bb_tr_vz[0] ); // z location of vertex in hall coordinates
			TLorentzVector P_beam( 0, 0, E_beam_final, E_beam_final ); //Mass of e negligable
			TLorentzVector k_prime( bb_tr_px[0], bb_tr_py[0], bb_tr_pz[0], bb_tr_p[0] );
			TLorentzVector P_targ( 0, 0, 0, Nucleon_mass );
		//Calculate q vector as beam momentum - scattered k
			TLorentzVector q = P_beam - k_prime; 

			e_prime_theta = acos( k_prime.Pz() / k_prime.E() );
			e_prime_phi = atan2( k_prime.Py(), k_prime.Px() );
			nucleon_phi = e_prime_phi + pi; //assume coplanarity	

			p_center = E_beam_final/(1.0 + (E_beam_final/Nucleon_mass)*( 1.0 - cos(e_prime_theta)));

			if( calc_method == 0 ){
				theta_pq_p_thresh = 0.0175;
				theta_pq_n_thresh = 0.0175;

				TLorentzVector P_gammaN = P_targ + q; //(-px, -py, ebeam - pz, Mp + ebeam - p)

				p_el = E_beam_final/( 1.0+E_beam_final/Mp*( 1.0-cos(e_prime_theta) ) );
				//double thetanucleon = acos( (E_corr - BBtr_p[0]*cos(etheta))/p_Nucleon ); //use elastic constraint on nucleon kinematics
				nu = E_beam_final - bb_tr_p[0];
				p_Nucleon = sqrt( pow(nu,2.0)+2.0*Nucleon_mass*nu );

				nucleon_theta = acos( (E_beam_final - bb_tr_pz[0])/p_Nucleon ); //use elastic constraint on nucleon kinematics

				pNhat.SetX(sin(nucleon_theta)*cos(nucleon_phi));
				pNhat.SetY(sin(nucleon_theta)*sin(nucleon_phi));
				pNhat.SetZ(cos(nucleon_theta));

				E_ep = sqrt( pow(Me,2) + pow(bb_tr_p[0],2) ); // Obtain the scattered electron energy
				h_E_ep->Fill( E_ep );

				p_ep = bb_tr_p[0];
				W = P_gammaN.M();
				W2 = pow(W, 2);	
				Q2 = 2*E_beam_final*E_ep*( 1-(bb_tr_pz[0]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta	
			}

			if( calc_method == 1){
				theta_pq_p_thresh = 0.04;
				theta_pq_n_thresh = 0.04;

				p_Beam = E_beam_final/(1.0 + E_beam_final/Nucleon_mass*(1.0 - cos(BB_theta)));
				pN = q + P_targ;
				pNhat = pN.Vect().Unit();
				Q2 = -q.M2();
				W2 = pN.M2();
				nu = q.E();
			}
			if( calc_method == 2){
				theta_pq_p_thresh = 0.065;
				theta_pq_n_thresh = 0.065;

				Q2 = e_kine_Q2;
				W2 = e_kine_W2;
				nu = e_kine_nu;
				p_Nucleon = sqrt( pow(nu, 2.0) + 2.0*(Nucleon_mass*nu) );
				nucleon_theta = acos( (P_beam.E() - p_center*cos(e_prime_theta))/p_Nucleon );
				pNhat.SetX( sin(nucleon_theta)*cos(nucleon_phi) );
				pNhat.SetY( sin(nucleon_theta)*sin(nucleon_phi) );
				pNhat.SetZ( cos(nucleon_theta) );
				pN.SetPxPyPzE( p_Nucleon*pNhat.X(), p_Nucleon*pNhat.Y(), p_Nucleon*pNhat.Z(), nu + P_targ.E() );
			}
			if( calc_method == 3){
				theta_pq_p_thresh = 0.02;
				theta_pq_n_thresh = 0.02;

				nu = P_beam.E() - p_center;
				p_Nucleon = sqrt( pow(nu, 2.0) + 2.0*(Nucleon_mass*nu) );
				nucleon_theta = acos( (P_beam.E() - p_center*cos(e_prime_theta)) / p_Nucleon);
				pNhat.SetX( sin(nucleon_theta)*cos(nucleon_phi) );
				pNhat.SetY( sin(nucleon_theta)*sin(nucleon_phi) );
				pNhat.SetZ( cos(nucleon_theta) );
				pN.SetPxPyPzE( p_Nucleon*pNhat.X(), p_Nucleon*pNhat.Y(), p_Nucleon*pNhat.Z(), nu + P_targ.E() );
				Q2 = 2.0*P_beam.E()*(k_prime.E())*(1.0 - cos(e_prime_theta));
				W2 = pow( Nucleon_mass, 2.0 ) + (2.0*(Nucleon_mass)*( P_beam.E() - k_prime.E() )) - Q2;
			}
			if( calc_method == 4){
				theta_pq_p_thresh = 0.065;
				theta_pq_n_thresh = 0.065;

				Q2 = e_kine_Q2;
				W2 = e_kine_W2;
				nu = e_kine_nu;
				p_Nucleon = sqrt( pow(nu, 2.0) + 2.0*(Nucleon_mass*nu) );
				nucleon_theta = acos( (P_beam.E() - p_center*cos(e_prime_theta)) / p_Nucleon);
				pNhat.SetX( sin(nucleon_theta)*cos(nucleon_phi) );
				pNhat.SetY( sin(nucleon_theta)*sin(nucleon_phi) );
				pNhat.SetZ( cos(nucleon_theta) );
				pN.SetPxPyPzE( p_Nucleon*pNhat.X(), p_Nucleon*pNhat.Y(), p_Nucleon*pNhat.Z(), nu + P_targ.E() );
				Q2 = 2.0*P_beam.E()*(k_prime.E())*(1.0 - cos(e_prime_theta));
				W2 = pow( Nucleon_mass, 2.0 ) + (2.0*(Nucleon_mass)*( P_beam.E() - k_prime.E() )) - Q2;
			}

			h_Q2->Fill( Q2 );
			// cout << "P_beam.E(): " <<  P_beam.E() << endl;
			// cout << "p_center: " <<  p_center << endl;
			// cout << "e_prime_theta: " << e_prime_theta << endl;
			// cout << "p_nucleon: " << p_Nucleon << endl;

			// cout << "nucleon theta: " << nucleon_theta << endl;
			// cout << "nucleon_phi: " << nucleon_phi << endl;
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

		//Expected neutron direction
			TVector3 Neutron_Direction = (HCAL_pos - vertex).Unit();

		//Expected proton direction
			//Need to incorporate deflection due to SBS magnet
			double Bdl = sbsfieldscale*maxSBSfield*Dgap/100.0;
			double Proton_Deflection = tan( 0.3*Bdl/q.Mag() )*(HCal_dist - (SBSdist + Dgap/2.0) ); 

			TVector3 Proton_Direction = (HCAL_pos + Proton_Deflection*HCAL_xaxis - vertex).Unit();

			// theta_pq_n = acos( Neutron_Direction.Dot( qvec_recon.Unit() ) );
			// theta_pq_p = acos( Proton_Direction.Dot( qvec_recon.Unit() ) );

		//SECONDARY METHOD FOR THETA_PQ
			theta_pq_p = acos( Proton_Direction.Dot( pNhat ));
			theta_pq_n = acos( Neutron_Direction.Dot( pNhat ));

			h_theta_pq_n->Fill(theta_pq_n);
			h_theta_pq_p->Fill(theta_pq_p);



			// cout << "Proton_Direction: " << Proton_Direction.X() << ", "<< Proton_Direction.Y() << ", "<< Proton_Direction.Z() << ", " << endl;
			// cout << "Neutron_Direction: " << Neutron_Direction.X() << ", "<< Neutron_Direction.Y() << ", "<< Neutron_Direction.Z() << ", " << endl;
			// cout << "pNhat: " << pNhat.X() << ", " << pNhat.Y() << ", " << pNhat.Z() << endl;

			// h_theta_pq_n_pNhat->Fill(thetapq_n);
			// h_theta_pq_p_pNhat->Fill(thetapq_p);

		//---------------------------
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



			//Get invariant mass transfer W from the four-momentum of the scattered nucleon
		
			W = pow(W2, 0.5);
			// W = P_gammaN.M();
			// W2 = pow(W, 2);				

		//IS NUCLEON ON HCAL FACE?
			if( acceptance_match ){
				hit_on_HCal = false;
				if( hcal_y>hcal_y_fmin && hcal_y<hcal_y_fmax && hcal_x>hcal_x_fmin && hcal_x<hcal_x_fmax ){
					hit_on_HCal = true;
				}
			}
			if( !acceptance_match ){
				hit_on_HCal = true;
			}


		//fill work histograms
			h_W->Fill( W );
			h_Wfull->Fill( W );

			if( (hit_on_HCal) && ( W2 < ( (h_W2->GetXaxis()->GetNbins())*(h_W2->GetXaxis()->GetBinWidth(1)) ) ) ){
				h_W2->Fill( W2 );	
			}
			h_W2full->Fill( W2 );
			

			if( !calc_W_only ){

				//Use the electron kinematics to predict the proton momedntum assuming elastic scattering on free proton at rest (will need to correct for fermi motion):
				E_pp = nu + Nucleon_mass; // Get energy of the proton
				E_nucleon = sqrt(pow(p_Nucleon,2)+pow(Nucleon_mass,2)); // Check on E_pp, same
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

///// PREVIOUSLY HAD ACCEPTANCE MATCH CUT

				//Check if the global cut WOULD HAVE failed. 
				//This would select elastics typically so we can check here to help make a filter on "real elastics"
				if( !( (nclus>0) && (PS_nclus>0) && (SH_nclus>0) && (abs(bb_tr_vz[0])<=0.075) && (bb_nhits[0]>4) && (bb_tr_n==1) && (bb_ps_e>0.15) && (hcal_e>0.01)) && ((bb_ps_e+bb_sh_e)/(bb_tr_p[0]) > 0.75) && ((bb_ps_e+bb_sh_e)/(bb_tr_p[0]) < 1.2) ){
					global_cut_fail = true;
				}
				else{
					global_cut_fail = false;
				}

				//Looking to cut out as much of the "real elastics" here to try and pull the BG shape out
				if( (theta_pq_p > theta_pq_p_thresh) && (dx<dx_min || dx>dx_max) && (dy<dy_min || dy>dy_max) && (global_cut_fail) && hit_on_HCal ){
					h_W2antiCut->Fill( W2 );
					h_W2antiCutFull->Fill( W2 );
				}

				//Looking to cut out all the background and only focus in on the "real elastics" signal here
				if( (theta_pq_p < theta_pq_p_thresh) && (dx>dx_min || dx<dx_max) && (dy>dy_min || dy<dy_max) && (!global_cut_fail) && hit_on_HCal ){

					h_W2elastics->Fill( W2 );
					h_W2elasticsFull->Fill( W2 );
				}
				if( (theta_pq_p < (0.30*theta_pq_p_thresh)) && (dx>dx_min || dx<dx_max) && (dy>dy_min || dy<dy_max) && (!global_cut_fail) && hit_on_HCal ){

					h_W2elastics_thr2->Fill( W2 );

				}
				if( (W2>W2_min) && (W2<W2_max) && (dy>dy_min) && (dy<dy_max) && (!global_cut_fail) && hit_on_HCal ){
					cout << "filling h_dx_allCuts... " << endl;
					h_dx_allCuts->Fill(dx);
				}

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

							//Looking to cut out as much of the "real elastics" here to try and pull the BG shape out
							if( (theta_pq_p > theta_pq_p_thresh) && (dx<dx_min || dx>dx_max) && (dy<dy_min || dy>dy_max) && (global_cut_fail) ){
								h_W2antiCut_fcut->Fill( W2 );
								h_W2antiCutFull_fcut->Fill( W2 );
							}

							//Looking to cut out all the background and only focus in on the "real elastics" signal here
							if( (thetapq_p < theta_pq_p_thresh) && (dx>dx_min || dx<dx_max) && (dy>dy_min || dy<dy_max) && (!global_cut_fail) ){

								h_W2elastics_fcut->Fill( W2 );
								h_W2elasticsFull_fcut->Fill( W2 );
							}
							if( (W2>W2_min) && (W2<W2_max) && (dy>dy_min) && (dy<dy_max) && (!global_cut_fail) ){
								cout << "filling h_dx_allCuts... " << endl;
								h_dx_allCuts_fcut->Fill(dx);
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

							//Looking to cut out as much of the "real elastics" here to try and pull the BG shape out
							if( (theta_pq_p > theta_pq_p_thresh) && (dx<dx_min || dx>dx_max) && (dy<dy_min || dy>dy_max) && (global_cut_fail) ){
								h_W2antiCut_fcut->Fill( W2 );
								h_W2antiCutFull_fcut->Fill( W2 );
							}

							//Looking to cut out all the background and only focus in on the "real elastics" signal here
							if( (theta_pq_p < theta_pq_p_thresh) && (dx>dx_min || dx<dx_max) && (dy>dy_min || dy<dy_max) && (!global_cut_fail) ){

								h_W2elastics_fcut->Fill( W2 );
								h_W2elasticsFull_fcut->Fill( W2 );
							}
							if( (W2>W2_min) && (W2<W2_max) && (dy>dy_min) && (dy<dy_max) && (!global_cut_fail) ){
								cout << "filling h_dx_allCuts... " << endl;
								h_dx_allCuts_fcut->Fill(dx);
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

		h_theta_pq_p->GetXaxis()->SetRangeUser(0, 0.3);
		cout << "theta_pq_p minimum: " << h_theta_pq_p->GetXaxis()->GetBinLowEdge(h_theta_pq_p->GetMinimumBin()) << endl;
		h_theta_pq_p->GetXaxis()->SetRangeUser(0, 0.6);

		cout << "---------------------------------------" << endl;
		cout << "-----Finished going through events-----" << endl;
		cout << "---------------------------------------" << endl;
		outfile->Write();	
	}

//--------------------------------------------------
//------------------- FITS -----------------------
//-----------------------------------------------
	cout << "---------------------------------------" << endl;
	cout << "Starting FITS" << endl;
	cout << "---------------------------------------" << endl;
	histo_infile = new TFile(outfile_name.Data(), "READ");

	if( fit_only ){

		polNfit = polN+1;

		if( interpolN != 99 ){
			interpolNfit = interpolN + 2;			
		}
		if( interpolN == 99 ){
			interpolNfit = 4;
		}


		cout << "polNfit: " << polNfit << endl;
		cout << "Fetching histograms. " << endl;

		h_W2 = static_cast<TH1D*>(histo_infile->Get("h_W2"));
		h_W2copy = static_cast<TH1D*>(histo_infile->Get("h_W2copy"));
		h_W2interpol = (TH1D*)h_W2->Clone("h_W2interpol");
		h_dx = static_cast<TH1D*>(histo_infile->Get("h_dx"));
		h_dy = static_cast<TH1D*>(histo_infile->Get("h_dy"));
		h_W2elastics = static_cast<TH1D*>(histo_infile->Get("h_W2elastics_thr2"));
		h_W2elastics_thr2 = static_cast<TH1D*>(histo_infile->Get("h_W2elastics_thr2"));
		h_W2elastics->Scale(2.0);
		h_W2elastics_fcut = static_cast<TH1D*>(histo_infile->Get("h_W2elastics_fcut"));
		// h_W2elastics_init = static_cast<TH1D*>(histo_infile->Get("h_W2elastics"));
		h_W2antiCut = static_cast<TH1D*>(histo_infile->Get("h_W2antiCut"));
		h_dxdy = static_cast<TH2D*>(histo_infile->Get("h_dxdy"));

		cout << "making histograms" << endl;
		h_W2_full_rand = new TH1D("h_W2_full_rand", Form("Invariant Mass W^2 Full Histo Filled Randomly- SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), nBins_W2, x_min_W2, x_max_W2);
		h_W2_residual = new TH1D("h_W2_residual", Form("Invariant Mass W^2 Fit Residual- SBS%i %i, %s; GeV", kine, sbsfieldscale, run_target.Data()), nBins_W2, x_min_W2, x_max_W2);

		TCanvas *c_W2_fullFit = new TCanvas("c_W2_fullFit", "c_W2_fullFit", 600, 500);
		h_W2->Draw();
		tf1_W2_fullFit = new TF1("tf1_W2_fullFit", fullFitFunction, W2_fit_min, W2_fit_max, 12);
		tf1_W2_fullFit->SetNpx(nBins_fit_W2);
		tf1_W2_fullFit->SetParName(0, "W2 Gaus Norm (No Rej)" );
		tf1_W2_fullFit->SetParName(1, "W2 Gaus Mean (No Rej)" );
		tf1_W2_fullFit->SetParName(2, "W2 Gaus Sigma (No Rej)");

		for(int param = 3; param <= polN+4; param++){
			tf1_W2_fullFit->SetParName(param, Form("W2 BG Pol%i p%i (No Rej)", polN, param-3));
		}

		tf1_W2_fullFit->SetParLimits(0, 4000, 6000);
		tf1_W2_fullFit->SetParLimits(1, 0.8, 0.9);
		tf1_W2_fullFit->SetParLimits(2, 0, 0.1);

		h_W2->Fit("tf1_W2_fullFit", "RM0+");
		tf1_W2_fullFit->GetParameters(&par[0]);


//---------------------------------------------
//----------------NO REJECT----------------------
//---------------------------------------------
		cout << "---------------------------------------" << endl;
		cout << " Plotting without reject" << endl;
		cout << "---------------------------------------" << endl;

		tf1_W2_sig = new TF1("tf1_W2_sig", fit_gaus, W2_fit_min, W2_fit_max, 3);
		tf1_W2_sig->SetNpx(nBins_fit_W2);
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

		if( polN == 99 ){
			tf1_W2_bg = new TF1("tf1_W2_bg", fit_expo, W2_fit_min, W2_fit_max, polN+1);
		}

		if( polN == 3 ){
			tf1_W2_bg = new TF1("tf1_W2_bg", fitPol3, W2_fit_min, W2_fit_max, polN+1);
		}

		if( polN == 4 ){
			tf1_W2_bg = new TF1("tf1_W2_bg", fitPol4, W2_fit_min, W2_fit_max, polN+1);
		}

		if( polN == 5 ){
			tf1_W2_bg = new TF1("tf1_W2_bg", fitPol5, W2_fit_min, W2_fit_max, polN+1);
		}

		if( polN == 6 ){
			tf1_W2_bg = new TF1("tf1_W2_bg", fitPol6, W2_fit_min, W2_fit_max, polN+1);
		}

		if( polN == 8 ){
			tf1_W2_bg = new TF1("tf1_W2_bg", fitPol8, W2_fit_min, W2_fit_max, polN+1);
		}

		tf1_W2_bg->SetNpx(nBins_fit_W2);
		for( int param = 0; param <= polN; param++ ){
			tf1_W2_bg->FixParameter(param, par[param+3]);
		}

		tf1_W2_bg->SetLineColor(9);
		h_W2->Fit("tf1_W2_bg", "RM0+");

		tf_sig_integral = h_tf_W2_sig->Integral();
		int W2_minus_bg_cnt = h_W2->GetEntries() - int(tf_sig_integral);
		TH1D *h_bg_rand = new TH1D("h_bg_rand", "Randomly filled BG from tf for BG", nBins_fit_W2, W2_fit_min, W2_fit_max);
		cout << "Filling bg_rand with: " << W2_minus_bg_cnt << " events." << endl;
		h_bg_rand->FillRandom("tf1_W2_bg", W2_minus_bg_cnt);

		h_tf_W2_bg = (TH1*)tf1_W2_bg->GetHistogram();
		h_tf_W2_bg->SetLineColor(9);
		h_tf_W2_bg->Scale(1.0/h_tf_W2_bg->Integral(), "width");
		h_tf_W2_bg->Scale(1.0/h_tf_W2_bg->Integral(), "height");
		h_tf_W2_bg->Scale( (h_bg_rand->GetMaximum())/(h_tf_W2_bg->GetMaximum()) );

		h_W2_sig = (TH1D*)h_W2->Clone("h_W2_sig");
		h_W2_sig->SetLineColor(6);
		for( int bin = 0; bin <= nBins_fit_W2; bin++ ){
			double x = bin*0.005;
			double min_x = par[1] - 4.5*par[2];
			double max_x = par[1] + 4.5*par[2];

			double sub_val = (h_W2->GetBinContent(bin)) - (h_tf_W2_bg->GetBinContent(bin));
			// h_W2_sig->SetBinContent(bin, sub_val);
			// if( sub_val >= 0){
			// 	h_W2_sig->SetBinContent( bin, (h_W2->GetBinContent(bin)) - (h_tf_W2_bg->GetBinContent(bin)) );
			// }
			// if( sub_val < 0 ){
			// 	h_W2_sig->SetBinContent( bin, 0 );
			// }

			if( x >= min_x && x <= max_x ){
				h_W2_sig->SetBinContent( bin, sub_val );	
			}
			else{
				h_W2_sig->SetBinContent( bin, 0 );
			}			
		}

		h_W2_sig->Draw("hist+same");

		double W2_sig_min_bin = h_W2_sig->GetXaxis()->FindBin(0.4);
		double W2_sig_max_bin = h_W2_sig->GetXaxis()->FindBin(1.2);

		sig_integral = h_W2_sig->Integral(W2_sig_min_bin, W2_sig_max_bin);

		TLegend *tl_no_reject = new TLegend(0.15, 0.70, 0.45, 0.85,"Full Range BG Fit");
		tl_no_reject->AddEntry(tf1_W2_fullFit, "Total W2 fit");
		tl_no_reject->AddEntry(h_W2_sig, "W2 peak after BG sub");
		tl_no_reject->AddEntry(h_tf_W2_bg, "BG from full range (no rej.)");
		// tl_no_reject->AddEntry((TObject*)0, Form("Fit Chi-Sq = %.2f", tf1_W2_fullFit->GetChisquare()), "");
		tl_no_reject->Draw("same");

		TPaveText *tp_W2 = new TPaveText(0.15, 0.58, 0.45, 0.69, "NDCbr");
		if( polN == 99 ){
			tp_W2->AddText("Fit with Expo");
		}
		if( polN != 99 ){
			tp_W2->AddText(Form("Fit with pol%i", polN));
		}
		tp_W2->AddText(Form("Estimated real W2 count: %i", int(sig_integral)));
		tp_W2->AddText(Form("Total fit chi-sq: %.2f", tf1_W2_fullFit->GetChisquare()));
		tp_W2->Draw("same");

	//------MANUALLY PLOTTING FITS/HISTOGRAMS/ETC
			c_W2_fullFit->cd();
			tf1_W2_fullFit->Draw("hist+same");
			h_tf_W2_bg->Draw("hist+same");

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

		for(int bin = 0; bin <= nBins_fit_W2; bin++){
			h_W2_residual->SetBinContent( bin, ( h_W2->GetBinContent(bin) - h_W2_full_rand->GetBinContent(bin) ) );
		}


//---------------------------------------------
//-----------------REJECTING--------------------
//---------------------------------------------
		cout << "---------------------------------------" << endl;
		cout << " Plotting WITH reject" << endl;
		cout << "---------------------------------------" << endl;

		if( fit_with_reject ){
			h_W2copy = (TH1D*)h_W2->Clone("h_W2copy");
		
			TCanvas *c_W2_fullFit_reject = new TCanvas("c_W2_fullFit_reject", "c_W2_fullFit_reject", 600, 500);
			h_W2copy->Draw();

	//----- For the reject BG fit we first need to find the BG fit with the reject in place
	//----- Then we take the parameters from that and fix the Full Fit BG parametes to match
			if( false ) {
				min_max_iter_cnt = 0;
				for( double min = 0.4; min < 0.55; min += 0.005 ){
					for(double max = 1.05; max < 1.20; max += 0.005 ){

						reject_min = min;
						reject_max = max;

						if( polN == 99 ){
							// reject_min = 0.50217680;
							// reject_max = 1.15;

							tf1_W2_bg_rej_init = new TF1("tf1_W2_bg_rej_init", fit_expo_with_reject, W2_fit_min, W2_fit_max, polNfit);
						}

						if( polN != 99 ){
							// reject_min = 0.50217680;
							// reject_max = 1.075;
						}

						if( polN == 3 ){
							tf1_W2_bg_rej_init = new TF1("tf1_W2_bg_rej_init", fitPol3_with_reject, W2_fit_min, W2_fit_max, polNfit);
						}

						if( polN == 4 ){
							tf1_W2_bg_rej_init = new TF1("tf1_W2_bg_rej_init", fitPol4_with_reject, W2_fit_min, W2_fit_max, polNfit);
						}

						if( polN == 5 ){
							tf1_W2_bg_rej_init = new TF1("tf1_W2_bg_rej_init", fitPol5_with_reject, W2_fit_min, W2_fit_max, polNfit);
						}

						if( polN == 6 ){
							tf1_W2_bg_rej_init = new TF1("tf1_W2_bg_rej_init", fitPol6_with_reject, W2_fit_min, W2_fit_max, polNfit);
						}

						if( polN == 8 ){
							tf1_W2_bg_rej_init = new TF1("tf1_W2_bg_rej_init", fitPol8_with_reject, W2_fit_min, W2_fit_max, polNfit);
						}
						tf1_W2_bg_rej_init->SetNpx(nBins_fit_W2);
						h_W2copy->Fit("tf1_W2_bg_rej_init", "RM0+");

						tf1_W2_bg_rej_init->GetParameters(&parRejSave[0]);
						tf1_W2_bg_rej_init->GetParameters(&parRej[0]);

						if( polN == 99 ){
							tf1_W2_bg_rej = new TF1("tf1_W2_bg_rej", fit_expo, W2_fit_min, W2_fit_max, polNfit);					
						}

						if( polN == 3 ){
							tf1_W2_bg_rej = new TF1("tf1_W2_bg_rej", fitPol3, W2_fit_min, W2_fit_max, polNfit);					
						}

						if( polN == 4 ){
							tf1_W2_bg_rej = new TF1("tf1_W2_bg_rej", fitPol4, W2_fit_min, W2_fit_max, polNfit);					
						}

						if( polN == 5 ){
							tf1_W2_bg_rej = new TF1("tf1_W2_bg_rej", fitPol5, W2_fit_min, W2_fit_max, polNfit);					
						}

						if( polN == 6 ){
							tf1_W2_bg_rej = new TF1("tf1_W2_bg_rej", fitPol6, W2_fit_min, W2_fit_max, polNfit);					
						}

						if( polN == 8 ){
							tf1_W2_bg_rej = new TF1("tf1_W2_bg_rej", fitPol8, W2_fit_min, W2_fit_max, polNfit);					
						}

						tf1_W2_bg_rej->SetNpx(nBins_fit_W2);

						for( int param = 0; param <= polN; param++ ){
							tf1_W2_bg_rej->FixParameter(param, parRej[param]);					
						}
						tf1_W2_bg_rej->SetLineColor(3);
						// tf1_W2_bg_pol8_rej->Draw("hist+same");

						//Need a new pol8 without the rejection to draw across reject zone:
						if( polN == 99 ){
							tf1_W2_bg_rej_full = new TF1("tf1_W2_bg_rej_full", fit_expo, W2_fit_min, W2_fit_max, polNfit);					
						}

						if( polN == 3 ){
							tf1_W2_bg_rej_full = new TF1("tf1_W2_bg_rej_full", fitPol3, W2_fit_min, W2_fit_max, polNfit);					
						}

						if( polN == 4 ){
							tf1_W2_bg_rej_full = new TF1("tf1_W2_bg_rej_full", fitPol4, W2_fit_min, W2_fit_max, polNfit);					
						}

						if( polN == 5 ){
							tf1_W2_bg_rej_full = new TF1("tf1_W2_bg_rej_full", fitPol5, W2_fit_min, W2_fit_max, polNfit);					
						}

						if( polN == 6 ){
							tf1_W2_bg_rej_full = new TF1("tf1_W2_bg_rej_full", fitPol6, W2_fit_min, W2_fit_max, polNfit);					
						}

						if( polN == 8 ){
							tf1_W2_bg_rej_full = new TF1("tf1_W2_bg_rej_full", fitPol8, W2_fit_min, W2_fit_max, polNfit);					
						}

						tf1_W2_bg_rej_full->SetNpx(nBins_fit_W2);

						for( int param = 0; param <= polN; param++ ){
							tf1_W2_bg_rej_full->FixParameter(param, parRej[param]);					
						}	

						tf1_W2_bg_rej_full->SetLineColor(9);	
						h_W2copy->Fit("tf1_W2_bg_rej_full", "RM0+");
				//----------
					// Using regular fullFitFucntion because we got the BG with Reject just above
						tf1_W2_fullFit_rej = new TF1("tf1_W2_fullFit_rej", fullFitFunction, W2_fit_min, W2_fit_max, polN+4);
						tf1_W2_fullFit_rej->SetNpx(nBins_fit_W2);
						tf1_W2_fullFit_rej->SetLineColor(2);

						tf1_W2_fullFit_rej->SetParName(0, "W2 Gaus Norm (Rej)" );
						tf1_W2_fullFit_rej->SetParName(1, "W2 Gaus Mean (Rej)" );
						tf1_W2_fullFit_rej->SetParName(2, "W2 Gaus Sigma (Rej)");

						for( int param = 3; param <= polN+4; param++ ){
							tf1_W2_fullFit_rej->SetParName(param, Form("W2 BG Pol%i p%i (Rej)", polN, param-3));
						}

					//Fix parameters from the BG_reject fit
						for(int param = 3; param < polN+4; param++ ){
							tf1_W2_fullFit_rej->FixParameter(param, parRej[param-3]);
						}

						tf1_W2_fullFit_rej->SetParLimits(0, 4000, 6000);
						tf1_W2_fullFit_rej->SetParLimits(1, 0.8, 0.9);
						tf1_W2_fullFit_rej->SetParLimits(2, 0, 0.1);

						h_W2copy->Fit("tf1_W2_fullFit_rej", "RM0");
						tf1_W2_fullFit_rej->GetParameters(&par[0]);

						tf1_W2_sig_rej = new TF1("tf1_W2_sig_rej", fit_gaus, W2_fit_min, W2_fit_max, 3);
						tf1_W2_sig_rej->SetNpx(nBins_fit_W2);
						tf1_W2_sig_rej->FixParameter(0, par[0]);
						tf1_W2_sig_rej->FixParameter(1, par[1]);
						tf1_W2_sig_rej->FixParameter(2, par[2]);

						h_W2copy->Fit("tf1_W2_fullFit_rej", "RM");

						iter_xmin_xmax_chi2[min_max_iter_cnt] = {reject_min, reject_max, tf1_W2_fullFit_rej->GetChisquare()};
						min_max_iter_cnt++;
					}

				}				
			}
			
			if( true ){
				if( polN == 99 ){
					// reject_min = 0.50217680;
					// reject_max = 1.15;

					reject_min = 0.545;
					reject_max = 1.06;

					tf1_W2_bg_rej_init = new TF1("tf1_W2_bg_rej_init", fit_expo_with_reject, W2_fit_min, W2_fit_max, polNfit);
				}

				if( polN != 99 ){
					reject_min = 0.54500000;
					reject_max = 1.07;
				}

				if( polN == 3 ){
					tf1_W2_bg_rej_init = new TF1("tf1_W2_bg_rej_init", fitPol3_with_reject, W2_fit_min, W2_fit_max, polNfit);
				}

				if( polN == 4 ){
					tf1_W2_bg_rej_init = new TF1("tf1_W2_bg_rej_init", fitPol4_with_reject, W2_fit_min, W2_fit_max, polNfit);
				}

				if( polN == 5 ){
					tf1_W2_bg_rej_init = new TF1("tf1_W2_bg_rej_init", fitPol5_with_reject, W2_fit_min, W2_fit_max, polNfit);
				}

				if( polN == 6 ){
					tf1_W2_bg_rej_init = new TF1("tf1_W2_bg_rej_init", fitPol6_with_reject, W2_fit_min, W2_fit_max, polNfit);
				}

				if( polN == 8 ){
					tf1_W2_bg_rej_init = new TF1("tf1_W2_bg_rej_init", fitPol8_with_reject, W2_fit_min, W2_fit_max, polNfit);
				}
				tf1_W2_bg_rej_init->SetNpx(nBins_fit_W2);
				h_W2copy->Fit("tf1_W2_bg_rej_init", "RM0+");

				tf1_W2_bg_rej_init->GetParameters(&parRejSave[0]);
				tf1_W2_bg_rej_init->GetParameters(&parRej[0]);

				if( polN == 99 ){
					tf1_W2_bg_rej = new TF1("tf1_W2_bg_rej", fit_expo, W2_fit_min, W2_fit_max, polNfit);					
				}

				if( polN == 3 ){
					tf1_W2_bg_rej = new TF1("tf1_W2_bg_rej", fitPol3, W2_fit_min, W2_fit_max, polNfit);					
				}

				if( polN == 4 ){
					tf1_W2_bg_rej = new TF1("tf1_W2_bg_rej", fitPol4, W2_fit_min, W2_fit_max, polNfit);					
				}

				if( polN == 5 ){
					tf1_W2_bg_rej = new TF1("tf1_W2_bg_rej", fitPol5, W2_fit_min, W2_fit_max, polNfit);					
				}

				if( polN == 6 ){
					tf1_W2_bg_rej = new TF1("tf1_W2_bg_rej", fitPol6, W2_fit_min, W2_fit_max, polNfit);					
				}

				if( polN == 8 ){
					tf1_W2_bg_rej = new TF1("tf1_W2_bg_rej", fitPol8, W2_fit_min, W2_fit_max, polNfit);					
				}

				tf1_W2_bg_rej->SetNpx(nBins_fit_W2);

				for( int param = 0; param <= polN; param++ ){
					tf1_W2_bg_rej->FixParameter(param, parRej[param]);					
				}
				tf1_W2_bg_rej->SetLineColor(3);
				// tf1_W2_bg_pol8_rej->Draw("hist+same");

				//Need a new pol8 without the rejection to draw across reject zone:
				if( polN == 99 ){
					tf1_W2_bg_rej_full = new TF1("tf1_W2_bg_rej_full", fit_expo, W2_fit_min, W2_fit_max, polNfit);					
				}

				if( polN == 3 ){
					tf1_W2_bg_rej_full = new TF1("tf1_W2_bg_rej_full", fitPol3, W2_fit_min, W2_fit_max, polNfit);					
				}

				if( polN == 4 ){
					tf1_W2_bg_rej_full = new TF1("tf1_W2_bg_rej_full", fitPol4, W2_fit_min, W2_fit_max, polNfit);					
				}

				if( polN == 5 ){
					tf1_W2_bg_rej_full = new TF1("tf1_W2_bg_rej_full", fitPol5, W2_fit_min, W2_fit_max, polNfit);					
				}

				if( polN == 6 ){
					tf1_W2_bg_rej_full = new TF1("tf1_W2_bg_rej_full", fitPol6, W2_fit_min, W2_fit_max, polNfit);					
				}

				if( polN == 8 ){
					tf1_W2_bg_rej_full = new TF1("tf1_W2_bg_rej_full", fitPol8, W2_fit_min, W2_fit_max, polNfit);					
				}

				tf1_W2_bg_rej_full->SetNpx(nBins_fit_W2);

				for( int param = 0; param <= polN; param++ ){
					tf1_W2_bg_rej_full->FixParameter(param, parRej[param]);					
				}	

				tf1_W2_bg_rej_full->SetLineColor(9);	
				h_W2copy->Fit("tf1_W2_bg_rej_full", "RM0+");
		//----------
			// Using regular fullFitFucntion because we got the BG with Reject just above
				tf1_W2_fullFit_rej = new TF1("tf1_W2_fullFit_rej", fullFitFunction, W2_fit_min, W2_fit_max, polN+4);
				tf1_W2_fullFit_rej->SetNpx(nBins_fit_W2);
				tf1_W2_fullFit_rej->SetLineColor(2);

				tf1_W2_fullFit_rej->SetParName(0, "W2 Gaus Norm (Rej)" );
				tf1_W2_fullFit_rej->SetParName(1, "W2 Gaus Mean (Rej)" );
				tf1_W2_fullFit_rej->SetParName(2, "W2 Gaus Sigma (Rej)");

				for( int param = 3; param <= polN+4; param++ ){
					tf1_W2_fullFit_rej->SetParName(param, Form("W2 BG Pol%i p%i (Rej)", polN, param-3));
				}

			//Fix parameters from the BG_reject fit
				for(int param = 3; param < polN+4; param++ ){
					tf1_W2_fullFit_rej->FixParameter(param, parRej[param-3]);
				}

				tf1_W2_fullFit_rej->SetParLimits(0, 4000, 6000);
				tf1_W2_fullFit_rej->SetParLimits(1, 0.8, 0.9);
				tf1_W2_fullFit_rej->SetParLimits(2, 0, 0.1);

				h_W2copy->Fit("tf1_W2_fullFit_rej", "RM0+");
				tf1_W2_fullFit_rej->GetParameters(&par[0]);
				iter_xmin_xmax_chi2[min_max_iter_cnt] = {reject_min, reject_max, tf1_W2_fullFit_rej->GetChisquare()};

				tf1_W2_sig_rej = new TF1("tf1_W2_sig_rej", fit_gaus, W2_fit_min, W2_fit_max, 3);
				tf1_W2_sig_rej->SetNpx(nBins_fit_W2);
				tf1_W2_sig_rej->SetParameter(0, par[0]);
				tf1_W2_sig_rej->SetParameter(1, par[1]);
				tf1_W2_sig_rej->SetParameter(2, par[2]);

				h_W2copy->Fit("tf1_W2_fullFit_rej", "RM+");
			}

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
			h_bg_rand_rej->FillRandom("tf1_W2_bg_rej_full", W2_minus_bg_cnt_rej);

			h_tf_W2_bg_rej = (TH1*)tf1_W2_bg_rej_full->GetHistogram();
			h_tf_W2_bg_rej->SetLineColor(9);
			h_tf_W2_bg_rej->Scale(1.0/h_tf_W2_bg_rej->Integral(), "width");
			h_tf_W2_bg_rej->Scale(1.0/h_tf_W2_bg_rej->Integral(), "height");
			h_tf_W2_bg_rej->Scale( (tf1_W2_bg_rej_full->GetMaximum())/(h_tf_W2_bg_rej->GetMaximum()) );

			// TH1D *h_W2_sig_rej = (TH1D*)h_W2copy->Clone("h_W2_sig_rej");
			h_W2_sig_rej = new TH1D("h_W2_sig_rej", "h_W2_sig_rej", h_W2copy->GetXaxis()->GetNbins(), 0, 1.4);
			h_W2_sig_rej->SetLineColor(6);
			for( int bin = 0; bin <= nBins_fit_W2; bin++ ){

				TAxis *W2_xaxis = h_W2->GetXaxis();
				// double x = bin*0.005;
				double min_x = par[1] - 4.5*par[2];
				double max_x = par[1] + 4.5*par[2];

				Int_t min_bin_x = W2_xaxis->FindBin(min_x);
				Int_t max_bin_x = W2_xaxis->FindBin(max_x);

				double sub_val = (h_W2copy->GetBinContent(bin)) - (h_tf_W2_bg_rej->GetBinContent(bin));

				// if( sub_val >= 0 ){
				// 	// h_W2_sig_rej->SetBinContent( bin, (h_W2copy->GetBinContent(bin)) - (h_tf_W2_bg_rej->GetBinContent(bin)) );	
				// 	h_W2_sig_rej->SetBinContent(bin, sub_val);
				// }
				// else{
				// 	h_W2_sig_rej->SetBinContent( bin, 0 );
				// }	

				if( bin >= min_bin_x && bin <= max_bin_x ){
					h_W2_sig_rej->SetBinContent( bin, sub_val );	
				}
				else{
					h_W2_sig_rej->SetBinContent( bin, 0 );
				}			
			}

			h_W2_sig_rej->Draw("hist+same");
			double W2_sig_rej_min_bin = h_W2_sig_rej->GetXaxis()->FindBin(0.4);
			double W2_sig_rej_max_bin = h_W2_sig_rej->GetXaxis()->FindBin(1.2);
			sig_integral_rej = h_W2_sig_rej->Integral(W2_sig_rej_min_bin, W2_sig_rej_max_bin);

			TLegend *tl_reject = new TLegend(0.15, 0.70, 0.45, 0.85, "BG with Reject Region");
			tl_reject->AddEntry(tf1_W2_fullFit_rej, "Total W2 fit");
			tl_reject->AddEntry(h_W2_sig_rej, "W2 peak after BG sub");
			tl_reject->AddEntry(h_tf_W2_bg_rej, "BG from full range (with rej.)");
			// tl_reject->AddEntry((TObject*)0, Form("Fit Chi-Sq = %.2f", tf1_W2_fullFit_rej->GetChisquare()), "");
			tl_reject->Draw("same");		

			TPaveText *tp_W2_rej = new TPaveText(0.15, 0.58, 0.45, 0.69, "NDCbr");
			if( polN == 99 ){
				tp_W2_rej->AddText("Fit with Expo & Reject()");
			}
			if( polN != 99 ){
				tp_W2_rej->AddText(Form("Fit with pol%i & Reject()", polN));
			}
			tp_W2_rej->AddText(Form("Estimated real W2 count: %i", int(sig_integral_rej)));
			tp_W2_rej->AddText(Form("Total fit chi-sq: %.2f", tf1_W2_fullFit_rej->GetChisquare()));
			tp_W2_rej->Draw("same");

	//------MANUALLY PLOTTING FITS/HISTOGRAMS/ETC
			c_W2_fullFit_reject->cd();
			tf1_W2_fullFit_rej->Draw("hist+same");
			h_tf_W2_bg_rej->Draw("hist+same");
			tf1_W2_bg_rej_full->SetLineColor(5);
			tf1_W2_bg_rej_full->Draw("hist+same");

			cout << "-------------------------------------" << endl;
			cout << " WITH REJECT POINTS " << endl;	
			cout << "-------------------------------------" << endl;
			cout << "tf_sig_integral_rej: " << tf_sig_integral_rej << ", sig_integral_rej: " << sig_integral_rej << endl;
			cout << "-------------------------------------" << endl;
		}

		// cout << " NO REJECT " << endl;
		// cout << "-------------------------------------" << endl;
		// cout << "tf_sig_integral: " << tf_sig_integral << ", sig_integral: " << sig_integral << endl;

	//--------------BGcut fits

		h_rej_sig_gaus = (TH1D*) h_tf_W2_sig_rej->Clone("h_rej_sig_gaus");
		h_rej_bg = new TH1D("h_rej_bg", "h_rej_bg", 280, 0, 1.4);

		for( int bin = 0; bin < 280; bin++ ){
			double sub_val = h_W2copy->GetBinContent(bin) - h_rej_sig_gaus->GetBinContent(bin);
			// double sub_val = h_W2copy->GetBinContent(bin) - h_tf_W2_sig->GetBinContent(bin);

			h_rej_bg->SetBinContent(bin, sub_val);
		}

		h_W2_bgCut = (TH1D*) h_W2->Clone("h_W2_bgCut");

		TCanvas *c_W2_bgCut = new TCanvas("c_W2_bgCut", "c_W2_bgCut", 600, 500);
		h_W2_bgCut->Draw();

		tf1_W2_bgCut_pol8 = new TF1("tf1_W2_bgCut_pol8", fitPol8, 0, 1.4, 9);
		tf1_W2_bgCut_pol8->SetNpx(280);
		tf1_W2_bgCut_pol8->SetLineColor(9);
		h_rej_bg->Fit("tf1_W2_bgCut_pol8", "RM0+");
		tf1_W2_bgCut_pol8->Draw("hist+same");

		h_tf_W2_bgCut_bg = (TH1*)tf1_W2_bgCut_pol8->GetHistogram();
		h_tf_W2_bgCut_bg->Scale(1.0/h_tf_W2_bgCut_bg->Integral(), "width");
		h_tf_W2_bgCut_bg->Scale(1.0/h_tf_W2_bgCut_bg->Integral(), "height");
		h_tf_W2_bgCut_bg->Scale(tf1_W2_bgCut_pol8->GetMaximum()/h_tf_W2_bgCut_bg->GetMaximum());

		h_W2_bgCut_bg = (TH1D*)h_tf_W2_bgCut_bg->Clone("h_W2_bgCut_bg");
	
		tf1_W2_bgCut_full = new TF1("tf1_W2_bgCut_full", fullBGCutFitFunction, 0.0, 1.4, 12);
		tf1_W2_bgCut_full->SetNpx(280);
		tf1_W2_bgCut_full->SetParName(0, "bgCut W2 Gaus Norm");
		tf1_W2_bgCut_full->SetParName(1, "bgCut W2 Gaus Mean");
		tf1_W2_bgCut_full->SetParName(2, "bgCut W2 Gaus Sigma");

		for(int param = 3; param < 12; param++ ){
			tf1_W2_bgCut_full->SetParName( param, Form("bgCut pol8 p%i", param-3));
			cout << "Fixing bgCut_full param " << param << " to : " << tf1_W2_bgCut_pol8->GetParameter(param-3) << endl;
			tf1_W2_bgCut_full->FixParameter(param, tf1_W2_bgCut_pol8->GetParameter(param-3));
		}

		tf1_W2_bgCut_full->SetParLimits(0, 0.4*h_W2_bgCut->GetMaximum(), h_W2_bgCut->GetMaximum());
		tf1_W2_bgCut_full->SetParLimits(1, 0.75, 0.95);
		tf1_W2_bgCut_full->SetParLimits(2, 0.0, 0.25);

		h_W2_bgCut->Fit("tf1_W2_bgCut_full", "RM+");

		bgCut_gaus_min = tf1_W2_bgCut_full->GetParameter(1) - 4.0*tf1_W2_bgCut_full->GetParameter(2);
		bgCut_gaus_max = tf1_W2_bgCut_full->GetParameter(1) + 4.0*tf1_W2_bgCut_full->GetParameter(2);

		cout << "bgCut max/min: " << bgCut_gaus_min << ", " << bgCut_gaus_max << endl;

		bgCut_full_integral = tf1_W2_bgCut_full->Integral(bgCut_gaus_min, bgCut_gaus_max)/h_W2_bgCut->GetBinWidth(1);
		sigma_bgCut_full_integral = tf1_W2_bgCut_full->IntegralError(bgCut_gaus_min, bgCut_gaus_max)/h_W2_bgCut->GetBinWidth(1);

		h_W2_bgCut_sig = new TH1D("h_W2_bgCut_sig", "h_W2_bgCut_sig", 280, 0.0, 1.4);
		h_W2_bgCut_sig->SetLineColor(6);

		h_W2_bgCut_bg_sub = new TH1D("h_W2_bgCut_bg_sub", "h_W2_bgCut_bg_sub", 280, 0.0, 1.4);
		h_W2_bgCut_bg_sub->SetLineColor(7);

		TAxis *tax_W2_bgCut_sig = h_W2_bgCut->GetXaxis();

		for(int bin = 0; bin <= 280; bin++ ){
			Int_t min_bin_x = tax_W2_bgCut_sig->FindBin(bgCut_gaus_min);
			Int_t max_bin_x = tax_W2_bgCut_sig->FindBin(bgCut_gaus_max);

			if( bin >= min_bin_x && bin <= max_bin_x ){
				double sub_val = h_W2_bgCut->GetBinContent(bin) - h_W2_bgCut_bg->GetBinContent(bin);

				h_W2_bgCut_sig->SetBinContent(bin, sub_val);
			}

			h_W2_bgCut_bg_sub->SetBinContent( bin, h_W2_bgCut->GetBinContent(bin) -  h_W2_bgCut_sig->GetBinContent(bin) );

		}

		h_W2_bgCut_sig->Draw("same");

		bgCut_sig_cnt = h_W2_bgCut_sig->Integral()/h_W2_bgCut_sig->GetBinWidth(1);

		tf1_W2_bgCut_gaus = new TF1("tf1_W2_bgCut_gaus", fit_gaus, 0.0, 1.4, 3);
		tf1_W2_bgCut_gaus->SetNpx(280);
		tf1_W2_bgCut_gaus->FixParameter(0, tf1_W2_bgCut_full->GetParameter(0));
		tf1_W2_bgCut_gaus->FixParameter(1, tf1_W2_bgCut_full->GetParameter(1));
		tf1_W2_bgCut_gaus->FixParameter(2, tf1_W2_bgCut_full->GetParameter(2));
		tf1_W2_bgCut_gaus->SetLineStyle(3);
		tf1_W2_bgCut_gaus->SetLineColor(15);
		h_W2_bgCut->Fit("tf1_W2_bgCut_gaus", "RM0");

		tf1_W2_bgCut_full->Draw("hist+same");
		tf1_W2_bgCut_gaus->Draw("hist+same");

		TF1 *tf1_bgCut_gaus_err = new TF1("tf1_bgCut_gaus_err", "gaus", bgCut_gaus_min, bgCut_gaus_max);
		tf1_bgCut_gaus_err->SetLineColor(5);
		tf1_bgCut_gaus_err->FixParameter(0, tf1_W2_bgCut_gaus->GetParameter(0));
		tf1_bgCut_gaus_err->FixParameter(1, tf1_W2_bgCut_gaus->GetParameter(1));
		tf1_bgCut_gaus_err->FixParameter(2, tf1_W2_bgCut_gaus->GetParameter(2));

		h_W2_bgCut_sig->Fit(tf1_bgCut_gaus_err, "R0");

		bgCut_sig_integral = tf1_bgCut_gaus_err->Integral(bgCut_gaus_min, bgCut_gaus_max)/h_W2_bgCut_sig->GetBinWidth(1);
		sigma_bgCut_sig_integral = tf1_bgCut_gaus_err->IntegralError(bgCut_gaus_min, bgCut_gaus_max)/h_W2_bgCut_sig->GetBinWidth(1);

		sigma_bgCut_TOTAL = sqrt( pow(sigma_bgCut_full_integral, 2) + pow( sigma_bgCut_sig_integral, 2) );

		cout << endl << endl << endl << endl << endl << "---------------------------------------" << endl;
		cout << "Gaus errors " << endl;

		tf1_W2_bgCut_gaus_err = new TF1("tf1_W2_bgCut_gaus_err", fit_gaus, 0.0, 1.4, 3);
		tfrp_W2_bgCut_gaus = h_W2_bgCut->Fit(tf1_W2_bgCut_gaus_err, "S");

		cov_W2_bgCut_gaus = tfrp_W2_bgCut_gaus->GetCovarianceMatrix();
		W2_bgCut_chi2 = tfrp_W2_bgCut_gaus->Chi2();

		W2_bgCut_gaus_par0 = tfrp_W2_bgCut_gaus->Parameter(0);
		W2_bgCut_gaus_parErr0 = tfrp_W2_bgCut_gaus->ParError(0);
		W2_bgCut_gaus_par1 = tfrp_W2_bgCut_gaus->Parameter(1);
		W2_bgCut_gaus_parErr1 = tfrp_W2_bgCut_gaus->ParError(1);
		W2_bgCut_gaus_par2 = tfrp_W2_bgCut_gaus->Parameter(2);
		W2_bgCut_gaus_parErr2 = tfrp_W2_bgCut_gaus->ParError(2);

		tfrp_W2_bgCut_gaus->Print("V");

		tf1_W2_bgCut_gaus_err->IntegralError(bgCut_gaus_min, bgCut_gaus_max, tfrp_W2_bgCut_gaus->GetParams(), cov_W2_bgCut_gaus.GetMatrixArray() );

		bgCut_sig_gaus_cnt = tf1_W2_bgCut_gaus->GetHistogram()->Integral()/tf1_W2_bgCut_gaus->GetHistogram()->GetBinWidth(1);

		TLegend *tl_W2_bgCut = new TLegend(0.15, 0.70, 0.45, 0.85,"BG-Subtracted Fit");
		tl_W2_bgCut->AddEntry(tf1_W2_bgCut_full, "Total W2 fit");
		tl_W2_bgCut->AddEntry(h_W2_bgCut_sig, "W2 peak from BG sub.");
		tl_W2_bgCut->AddEntry(tf1_W2_bgCut_gaus, "Gaussian fit to W2");
		tl_W2_bgCut->AddEntry(tf1_W2_bgCut_pol8, "BG fit");

		tl_W2_bgCut->Draw("same");

		TPaveText *tp_W2_bgCut = new TPaveText(0.15, 0.25, 0.45, 0.69, "NDCbr");
		tp_W2_bgCut->AddText(Form("BG-sub W2 cnt = %i", int(bgCut_sig_cnt)));
		tp_W2_bgCut->AddText(Form("W2 Gaus-fit cnt = %i", int(bgCut_sig_gaus_cnt)));
		tp_W2_bgCut->AddText(Form("Total fit chi-sq: %.2f", tf1_W2_bgCut_full->GetChisquare()));
		tp_W2_bgCut->AddText("----------------------");
		tp_W2_bgCut->Draw("same");


cout << endl << endl << endl << endl << endl << "---------------------------------------" << endl;
//----------------------------
	//--------INTERPOLATING FITS

		cout << "---------------------------------------" << endl;
		cout << " Plotting WITH Interpolation" << endl;
		cout << "---------------------------------------" << endl;
		double interpol_x_min = 0.0;
		double interpol_x_max = 1.4;
		//-----
		// h_W2elastics = (TH1D*)h_W2_sig->Clone("h_W2elastics");
		TCanvas *c_W2_interpol = new TCanvas("c_W2_interpol", "c_W2_interpol", 600, 500);
		h_W2interpol->Draw();

	//--------------------------------------	
		h_W2elastics_init = (TH1D*)h_W2elastics->Clone("h_W2elastics_init");
		h_W2elastics = (TH1D*)h_W2elastics->Clone("h_W2elastics");

		h_W2elastics_corr = (TH1D*)h_W2elastics->Clone("h_W2elastics_corr");


		// h_W2elastics_init = (TH1D*)h_W2elastics_fcut->Clone("h_W2elastics_init");
		// h_W2elastics = (TH1D*)h_W2elastics_fcut->Clone("h_W2elastics");
	//--------------------------------------

		tf1_W2_elastics_init = new TF1("tf1_W2_elastics_init", fit_gaus, interpol_x_min, interpol_x_max, 3);
		tf1_W2_elastics_init->SetNpx(nBins_fit_W2);
		tf1_W2_elastics_init->SetParLimits(0, 0.99*h_W2elastics_init->GetMaximum(), h_W2elastics_init->GetMaximum());
		tf1_W2_elastics_init->SetParLimits(1, 0.8, 0.825);
		tf1_W2_elastics_init->SetParLimits(2, 0, 0.090);
		h_W2elastics_corr->Fit("tf1_W2_elastics_init", "RM0+");
		h_W2elastics_corr->SetLineColor(3);

		// // // tf1_W2_elastics->GetParameters(&parInterpol[0]);
		tf1_W2_elastics_init->GetParameters(&par[0]);

		h_tf_W2elastics_init = (TH1*)tf1_W2_elastics_init->GetHistogram();
		h_tf_W2elastics_init->Scale(1.0/h_tf_W2elastics_init->Integral(), "width");
		h_tf_W2elastics_init->Scale(1.0/h_tf_W2elastics_init->Integral(), "height");
		h_tf_W2elastics_init->Scale(par[0]/h_tf_W2elastics_init->GetMaximum());
		h_W2elastics_gaus = (TH1D*)h_tf_W2elastics_init->Clone("h_W2elastics_gaus");

		TAxis *W2elastics_init_xaxis = h_W2elastics_gaus->GetXaxis();
		Int_t max_bin_x = W2elastics_init_xaxis->FindBin(0.965);

		for( int bin = max_bin_x; bin < h_W2elastics_gaus->GetXaxis()->GetNbins(); bin++ ){
			h_W2elastics_corr->SetBinContent(bin, h_W2elastics_gaus->GetBinContent(bin));
		}	
		h_W2elastics_init = (TH1D*)h_W2elastics_corr->Clone("h_W2elastics_init");
		// h_W2elastics = (TH1D*)h_W2elastics_corr->Clone("h_W2elastics");

		// h_tf_W2elastics_init = (TH1*)tf1_W2_elastics->GetHistogram();
		// h_tf_W2elastics_init->Scale(1.0/h_tf_W2elastics_init->Integral(), "width");
		// h_tf_W2elastics_init->Scale(1.0/h_tf_W2elastics_init->Integral(), "height");
		// h_tf_W2elastics_init->Scale(tf1_W2_elastics->GetParameter(0)/h_tf_W2elastics_init->GetMaximum());

		// h_W2elastics = (TH1D*)h_tf_W2elastics_init->Clone("h_W2elastics");
		h_W2elastics_thr2->Scale(4.0);
		h_W2elastics = (TH1D*)h_W2elastics_thr2->Clone("h_W2elastics");

		h_W2elastics->GetXaxis()->SetRangeUser(interpol_x_min, interpol_x_max);

		// tf1_W2_interpolate = new TF1("tf1_W2_interpolate", W2FullInterpol, W2_fit_min, W2_fit_max, polN+2);
		tf1_W2_interpolate = new TF1("tf1_W2_interpolate", W2FullInterpol, interpol_x_min, interpol_x_max, interpolNfit);
		tf1_W2_interpolate->SetNpx(nBins_fit_W2);
		tf1_W2_interpolate->SetParName(0, "W2 Elastic Sig Norm" );

		for( int param = 1; param < interpolN+2; param++ ){
			tf1_W2_interpolate->SetParName(param, Form("W2 BG Pol%i p%i (Interpol)", interpolN, param-1));
		}
		h_W2interpol->Fit("tf1_W2_interpolate", "RM0+");
		tf1_W2_interpolate->GetParameters(&parInterpol[0]);
		tf1_W2_interpolate->GetParameters(&par[0]);

		if( interpolN == 99 ){
			tf1_W2_interpol_bg = new TF1("tf1_W2_interpol_bg", fit_expo, interpol_x_min, interpol_x_max, interpolNfit);					
		}

		if( interpolN == 3 ){
			tf1_W2_interpol_bg = new TF1("tf1_W2_interpol_bg", fitPol3, interpol_x_min, interpol_x_max, interpolNfit);					
		}

		if( interpolN == 4 ){
			tf1_W2_interpol_bg = new TF1("tf1_W2_interpol_bg", fitPol4, interpol_x_min, interpol_x_max, interpolNfit);					
		}

		if( interpolN == 5 ){
			tf1_W2_interpol_bg = new TF1("tf1_W2_interpol_bg", fitPol5, interpol_x_min, interpol_x_max, interpolNfit);					
		}

		if( interpolN == 6 ){
			tf1_W2_interpol_bg = new TF1("tf1_W2_interpol_bg", fitPol6, interpol_x_min, interpol_x_max, interpolNfit);					
		}

		if( interpolN == 8 ){
			tf1_W2_interpol_bg = new TF1("tf1_W2_interpol_bg", fitPol8, interpol_x_min, interpol_x_max, interpolNfit);					
		}

		tf1_W2_interpol_bg->SetNpx(nBins_fit_W2);

		for( int param = 0; param < interpolNfit - 1; param++ ){
			cout << "Setting param: " << param << " to parInterpol[" << param + 1 << "] = " << parInterpol[param+1] << endl;
			tf1_W2_interpol_bg->FixParameter(param, parInterpol[param+1]);					
		}	

		tf1_W2_interpol_bg->SetLineColor(9);	

		h_tf_W2_interpol_bg = (TH1*)tf1_W2_interpol_bg->GetHistogram();
		h_tf_W2_interpol_bg->Scale(1.0/h_tf_W2_interpol_bg->Integral(), "width");
		h_tf_W2_interpol_bg->Scale(1.0/h_tf_W2_interpol_bg->Integral(), "height");
		h_tf_W2_interpol_bg->Scale(tf1_W2_interpol_bg->GetMaximum()/h_tf_W2_interpol_bg->GetMaximum());

		h_W2_interpol_bg = (TH1D*)h_tf_W2_interpol_bg->Clone("h_W2_interpol_bg");

		h_W2interpol->Fit("tf1_W2_interpol_bg", "RM+");

		// h_W2elastics_init->Draw("same");

		tf1_W2_interpolate_full = new TF1("tf1_W2_interpolate_full", fullInterpolFitFunction, interpol_x_min, interpol_x_max, interpolNfit+2);
		tf1_W2_interpolate_full->SetNpx( int((W2_fit_max - W2_fit_min)*200.0) );

		tf1_W2_interpolate_full->SetParLimits(0, 0.5*h_W2interpol->GetMaximum(), h_W2interpol->GetMaximum());
		tf1_W2_interpolate_full->SetParLimits(1, 0.81, 0.845);
		tf1_W2_interpolate_full->SetParLimits(2, 0, 0.1);

		tf1_W2_interpolate_full->SetParName(0, "Interpolated Gaus Norm");
		tf1_W2_interpolate_full->SetParName(1, "Interpolated Gaus Mean");
		tf1_W2_interpolate_full->SetParName(2, "Interpolated Gaus Sigma");

		for( int param = 3; param <= (interpolN + 4); param++ ){
			tf1_W2_interpolate_full->FixParameter(param, parInterpol[param - 2]);
			tf1_W2_interpolate_full->SetParName(param, Form("Interpol pol%i p%i", interpolN, param-3));
		}

		h_W2interpol->Fit("tf1_W2_interpolate_full", "R+");

		tf1_W2_interpolate_full->GetParameters(par);

		TF1 *tf1_W2_interpolate_gaus = new TF1("tf1_W2_interpolate_gaus", fit_gaus, interpol_x_min, interpol_x_max, 3);
		tf1_W2_interpolate_gaus->SetLineColor(5);
		tf1_W2_interpolate_gaus->SetNpx(280);
		tf1_W2_interpolate_gaus->FixParameter(0, par[0]);
		tf1_W2_interpolate_gaus->FixParameter(1, par[1]);
		tf1_W2_interpolate_gaus->FixParameter(2, par[2]);

		h_W2interpol->Fit("tf1_W2_interpolate_gaus", "R0+");

		h_W2_interpol_sig = new TH1D("h_W2_interpol_sig", "h_W2_interpol_sig", h_W2interpol->GetXaxis()->GetNbins(), interpol_x_min, interpol_x_max);

		for(int bin = 0; bin <= h_W2interpol->GetXaxis()->GetNbins(); bin++ ){
			double sub_val = h_W2interpol->GetBinContent(bin) - h_W2_interpol_bg->GetBinContent(bin);

			if( sub_val >= 0 ){
				h_W2_interpol_sig->SetBinContent(bin, sub_val);
			}
		}
		h_W2_interpol_sig->SetLineColor(6);
		h_W2_interpol_sig->Draw("hist+same");
		tf1_W2_interpolate_gaus->Draw("hist+same");

		int W2_interpol_sig_min_bin = h_W2_interpol_sig->GetXaxis()->FindBin(0.4);
		int W2_interpol_sig_max_bin = h_W2_interpol_sig->GetXaxis()->FindBin(1.2);

		double W2_interpol_sig_cnt = h_W2_interpol_sig->Integral(W2_interpol_sig_min_bin, W2_interpol_sig_min_bin);
		cout << "W2_interpol_sig_min_bin: " << W2_interpol_sig_min_bin << endl;
		cout << "W2_interpol_sig_max_bin: " << W2_interpol_sig_max_bin << endl;
		cout << "W2_interpol_sig_cnt: " << W2_interpol_sig_cnt << endl;

		TLegend *tl_W2_interpol = new TLegend(0.15, 0.70, 0.45, 0.85,"Full Range Interpolated BG Fit");
		tl_W2_interpol->AddEntry(tf1_W2_interpolate_full, "Total W2 Interpolated fit");
		tl_W2_interpol->AddEntry(h_W2_interpol_sig, "W2 peak after BG sub");
		tl_W2_interpol->AddEntry(tf1_W2_interpol_bg, "BG fit");
		// tl_no_reject->AddEntry((TObject*)0, Form("Fit Chi-Sq = %.2f", tf1_dx_fullFit->GetChisquare()), "");
		tl_W2_interpol->Draw("same");

		TPaveText *tp_W2_interpol = new TPaveText(0.15, 0.58, 0.45, 0.69, "NDCbr");
		tp_W2_interpol->AddText(Form("Interpolated fit with pol%i", interpolN));
		tp_W2_interpol->AddText(Form("Estimated real W2 count: %i", int(W2_interpol_sig_cnt)));
		tp_W2_interpol->AddText(Form("Total fit chi-sq: %.2f", tf1_W2_interpolate_full->GetChisquare()));
		tp_W2_interpol->Draw("same");

//------------------------ dxdy stuff -----------
		cout << "---------------------------------------" << endl;
		cout << " Plotting dxdy" << endl;
		cout << "---------------------------------------" << endl;

		TCanvas *c_dxdy;
		if( plot_dxdy ){
			c_dxdy = new TCanvas("c_dxdy", "c_dxdy", 600, 500);
			h_dxdy->Draw("colz");			
		}


//------------------------------------------------------
//---------------------- DY STUFF ----------------------
//------------------------------------------------------

		cout << "---------------------------------------" << endl;
		cout << " Plotting DY" << endl;
		cout << "---------------------------------------" << endl;

		if( true ){
			int dypolNfit = dypolN+1;
			double dymax = 0.30;
			cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
			cout << "        Working on dy plots   " << endl;

			h_dxdy->GetYaxis()->SetRangeUser(-0.3, 0.3);
			h_dxdy->GetXaxis()->SetRangeUser(-0.4, 0.4);
			
			h_dy_subrange = h_dxdy->ProjectionX();
			h_dxdy->GetYaxis()->SetRangeUser(-2.5, 2.5);
			h_dxdy->GetXaxis()->SetRangeUser(-1.5, 1.5);

			TCanvas *c_dy = new TCanvas("c_dy", "c_dy", 600, 500);
			h_dy_subrange->GetYaxis()->SetRangeUser(0, 1.05*h_dy_subrange->GetMaximum());
			h_dy_subrange->Draw();

			tf1_dy_fullFit = new TF1("tf1_dy_fullFit", fullDYFitFunction, -0.4, 0.4, 6);

			tf1_dy_fullFit->SetNpx(400);
			tf1_dy_fullFit->SetParName(0, "dy Gaus Norm");
			tf1_dy_fullFit->SetParName(1, "dy Gaus Mean");
			tf1_dy_fullFit->SetParName(2, "dy Gaus Sigma");
			tf1_dy_fullFit->SetParName(3, "dy pol3 p0");
			tf1_dy_fullFit->SetParName(4, "dy pol3 p1");
			// tf1_dy_fullFit->SetParName(5, "dy pol3 p2");


			tf1_dy_fullFit->SetParLimits(0, 0.3*h_dy_subrange->GetMaximum(), 0.9*h_dy_subrange->GetMaximum());
			tf1_dy_fullFit->SetParLimits(1, -0.05, 0.05);
			tf1_dy_fullFit->SetParLimits(2, 0.0, 0.11);
			tf1_dy_fullFit->SetParLimits(4, 0.0, 1000);
			// tf1_dy_fullFit->FixParameter(3, h_dy_subrange->GetBinContent(h_dy_subrange->GetMinimumBin()));

			h_dy_subrange->Fit("tf1_dy_fullFit", "R");

			tf1_dy_fullFit->GetParameters(parDY);

		//Signal Gaussian
			tf1_dy_sig = new TF1("tf1_dy_sig", fit_gaus, -0.4, 0.4, 3);
			tf1_dy_sig->SetNpx(400);
			tf1_dy_sig->FixParameter(0, parDY[0]);
			tf1_dy_sig->FixParameter(1, parDY[1]);
			tf1_dy_sig->FixParameter(2, parDY[2]);
			tf1_dy_sig->SetLineColor(6);
			tf1_dy_sig->Draw("hist+same");

			dy_sig_mean = parDY[1];
			dy_sig_sigma = parDY[2];

			dy_sig_min = dy_sig_mean - 4.0*dy_sig_sigma;
			dy_sig_max = dy_sig_mean + 4.0*dy_sig_sigma;

			h_tf_dy_sig = (TH1*)tf1_dy_sig->GetHistogram();
			h_tf_dy_sig->Scale(1.0/h_tf_dy_sig->Integral(), "width");
			h_tf_dy_sig->Scale(1.0/h_tf_dy_sig->Integral(), "height");
			h_tf_dy_sig->Scale(parDY[0]/h_tf_dy_sig->GetMaximum());
			h_tf_dy_sig->SetLineColor(6);

			double dy_sig_integral = h_tf_dy_sig->Integral();

			int dy_bg_rand_cnt = int( h_dy_subrange->GetEntries() - dy_sig_integral);

			cout << "dy_sig_integral: " << dy_sig_integral << ", dy_bg_rand_cnt: " << dy_bg_rand_cnt << endl;

		//BG Pol3

			tf1_dy_bg = new TF1("tf1_dy_bg", fitPol2, -0.4, 0.4, 3);
			tf1_dy_bg->SetNpx(400);
			tf1_dy_bg->SetLineColor(9);
			tf1_dy_bg->SetParName(0, "dy bg pol3 p0");
			tf1_dy_bg->SetParName(1, "dy bg pol3 p1");
			tf1_dy_bg->SetParName(2, "dy bg pol3 p2");
			// tf1_dy_bg->SetParName(3, "dy bg pol3 p3");

			tf1_dy_bg->FixParameter(0, parDY[3]);
			tf1_dy_bg->FixParameter(1, parDY[4]);
			tf1_dy_bg->FixParameter(2, parDY[5]);
			// tf1_dy_bg->FixParameter(3, parDY[6]);

			tf1_dy_bg->Draw("hist+same");

			TLegend *tl_dy = new TLegend(0.60, 0.50, 0.85, 0.70,"Sub-range dy Fit");
			tl_dy->AddEntry(tf1_dy_fullFit, "Total dy fit");
			tl_dy->AddEntry(tf1_dy_sig, "dy peak after BG sub");
			tl_dy->AddEntry(tf1_dy_bg, "BG fit");
			// tl_no_reject->AddEntry((TObject*)0, Form("Fit Chi-Sq = %.2f", tf1_dx_fullFit->GetChisquare()), "");
			tl_dy->Draw("same");

			TPaveText *tp_dy = new TPaveText(0.60, 0.39, 0.85, 0.49, "NDCbr");
			tp_dy->AddText(Form("dy mean: %.3f, dy sigma: %0.3f", dy_sig_mean, dy_sig_sigma));
			tp_dy->AddText(Form("dy_min: %.3f, dy_max: %0.3f", dy_sig_min, dy_sig_max));
			tp_dy->Draw("same");

			TLine *tl_dy_min = new TLine(dy_sig_min, -2.5, dy_sig_min, 2.5);
			tl_dy_min->SetLineColor(2);
			tl_dy_min->SetLineWidth(3);
			tl_dy_min->SetLineStyle(4);

			TLine *tl_dy_max = new TLine(dy_sig_max, -2.5, dy_sig_max, 2.5);
			tl_dy_max->SetLineColor(2);
			tl_dy_max->SetLineWidth(3);
			tl_dy_max->SetLineStyle(4);

			if( plot_dxdy && plot_dxdy_anticut ){
				c_dxdy->cd();
				tl_dy_min->Draw("same");
				tl_dy_max->Draw("same");				
			}


		}

//------------------------------------------------------
//---------------------- DX STUFF ----------------------
//------------------------------------------------------


		cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
		cout << "        Working on dx plots   " << endl;
		TCanvas *c_dx = new TCanvas("c_dx", "c_dx", 600, 500);

		h_dxdy->GetXaxis()->SetRangeUser(dy_sig_min, dy_sig_max);
		h_dxdy->GetYaxis()->SetRangeUser(-0.5, 0.5);
		h_dx_subrange = h_dxdy->ProjectionY();
		h_dxdy->GetXaxis()->SetRangeUser(-1.5, 1.5);
		h_dxdy->GetYaxis()->SetRangeUser(-2.5, 2.5);

		h_dx_subrange->GetYaxis()->SetRangeUser(0, 1.05*h_dx_subrange->GetMaximum());
		h_dx_subrange->Draw();

		tf1_dx_fullFit = new TF1("tf1_dx_fullFit", fullDXFitFunction, -0.5, 0.5, 7);
		tf1_dx_fullFit->SetNpx(200);
		tf1_dx_fullFit->SetParName(0, "dx Gaus Norm (No Rej)" );
		tf1_dx_fullFit->SetParName(1, "dx Gaus Mean (No Rej)" );
		tf1_dx_fullFit->SetParName(2, "dx Gaus Sigma (No Rej)");
		tf1_dx_fullFit->SetParName(3, "dx BG Pol3 p0 (No Rej)");
		tf1_dx_fullFit->SetParName(4, "dx BG Pol3 p1 (No Rej)");
		tf1_dx_fullFit->SetParName(5, "dx BG Pol3 p2 (No Rej)");
		tf1_dx_fullFit->SetParName(6, "dx BG Pol3 p3 (No Rej)");

		tf1_dx_fullFit->SetParLimits(0, 0.3*h_dx->GetMaximum(), h_dx->GetMaximum());
		tf1_dx_fullFit->SetParLimits(1, -0.05, 0.05);
		tf1_dx_fullFit->SetParLimits(2, 0, 0.07);

		h_dx_subrange->Fit("tf1_dx_fullFit", "R+");

		dx_full_integral = tf1_dx_fullFit->Integral(-0.5, 0.5)/tf1_dx_fullFit->GetHistogram()->GetBinWidth(1);
		sigma_dx_full_integral = tf1_dx_fullFit->IntegralError(-0.5, 0.5)/tf1_dx_fullFit->GetHistogram()->GetBinWidth(1);
		
		tf1_dx_fullFit->GetParameters(&par[0]);		

		tf1_dx_sig = new TF1("tf1_dx_sig", fit_gaus, -0.5, 0.5, 3);
		tf1_dx_sig->SetNpx(200);
		tf1_dx_sig->FixParameter(0, par[0]);
		tf1_dx_sig->FixParameter(1, par[1]);
		tf1_dx_sig->FixParameter(2, par[2]);

		h_tf_dx_sig = (TH1*)tf1_dx_sig->GetHistogram();
		h_tf_dx_sig->SetLineColor(6);
		h_tf_dx_sig->Scale(1./h_tf_dx_sig->Integral(), "width");
		h_tf_dx_sig->Scale(1./h_tf_dx_sig->Integral(), "height");
		h_tf_dx_sig->Scale(par[0]/(h_tf_dx_sig->GetMaximum()));

		// h_tf_dx_sig->Draw("hist+same");

		tf1_dx_bg = new TF1("tf1_dx_bg", fitPol3, -0.5, 0.5, 4);
		tf1_dx_bg->SetNpx(200);
		tf1_dx_bg->SetLineColor(9);
		tf1_dx_bg->FixParameter(0, par[3]);
		tf1_dx_bg->FixParameter(1, par[4]);
		tf1_dx_bg->FixParameter(2, par[5]);
		tf1_dx_bg->FixParameter(3, par[6]);
		tf1_dx_bg->Draw("hist+same");

		h_dx_subrange->Fit("tf1_dx_bg", "R0+");

		TF1 *tf1_dx_bg_pol3_err = new TF1("tf1_dx_bg_pol3_err", "pol3", -0.5, 0.5);
		tf1_dx_bg_pol3_err->SetNpx(200);
		tf1_dx_bg_pol3_err->FixParameter(0, tf1_dx_bg->GetParameter(0));
		tf1_dx_bg_pol3_err->FixParameter(1, tf1_dx_bg->GetParameter(1));
		tf1_dx_bg_pol3_err->FixParameter(2, tf1_dx_bg->GetParameter(2));
		tf1_dx_bg_pol3_err->FixParameter(3, tf1_dx_bg->GetParameter(3));

		h_dx_subrange->Fit(tf1_dx_bg_pol3_err, "RS");

		dx_bg_integral = tf1_dx_bg_pol3_err->Integral(-0.5, 0.5);

		sigma_dx_bg_integral = tf1_dx_bg_pol3_err->IntegralError(-0.5, 0.5);
		// sigma_dx_bg_integral = 151.17716;


		double tf1_dx_bg_par0 = tf1_dx_bg->GetParameter(0);
		double tf1_dx_bg_par1 = tf1_dx_bg->GetParameter(1);
		double tf1_dx_bg_par2 = tf1_dx_bg->GetParameter(2);
		double tf1_dx_bg_par3 = tf1_dx_bg->GetParameter(3);

		h_tf_dx_bg = (TH1*)tf1_dx_bg->GetHistogram();
		h_tf_dx_bg->Scale(1.0/h_tf_dx_bg->Integral(), "width");
		h_tf_dx_bg->Scale(1.0/h_tf_dx_bg->Integral(), "height");
		h_tf_dx_bg->Scale(tf1_dx_bg->GetMaximum()/h_tf_dx_bg->GetMaximum());

		h_dx_bg = (TH1D*)h_tf_dx_bg->Clone("h_dx_bg");

		h_dx_sig = new TH1D("h_dx_sig", "h_dx_sig", 200, -0.5, 0.5);
		h_dx_sig->SetLineColor(6);
		for( int bin = 0; bin <= h_dx_sig->GetXaxis()->GetNbins(); bin++ ){
			double sub_val = h_dx_subrange->GetBinContent(bin) - h_dx_bg->GetBinContent(bin);

			if(sub_val >= 0){
				h_dx_sig->SetBinContent(bin, sub_val);
			}
			else{
				h_dx_sig->SetBinContent(bin, 0);
			}
		}
		h_dx_sig->Draw("hist+same");

		// tf1_dx_bg->Scale(1.0/tf1_dx_bg->Integral(), "width");
		// tf1_dx_bg->Scale(1.0/tf1_dx_bg->Integral(), "height");

		h_dx_subrange_sig = new TH1D("h_dx_subrange_sig", "h_dx_subrange_sig", 200, -0.5, 0.5);

		dx_sig_gaus_min = tf1_dx_sig->GetParameter(1) - 4.0*tf1_dx_sig->GetParameter(2);
		dx_sig_gaus_max = tf1_dx_sig->GetParameter(1) + 4.0*tf1_dx_sig->GetParameter(2);

		for( int bin = 0; bin <= 200; bin++ ){
			Int_t min_bin_x = h_dx_subrange_sig->FindBin(dx_sig_gaus_min);
			Int_t max_bin_x = h_dx_subrange_sig->FindBin(dx_sig_gaus_max);

			if( bin >= min_bin_x && bin <= max_bin_x ){
				double sub_val = h_dx_subrange->GetBinContent(bin) - h_dx_bg->GetBinContent(bin);
				
				h_dx_subrange_sig->SetBinContent(bin, sub_val);
			}
		}

		TF1 *tf1_dx_sig_err = new TF1("tf1_dx_sig_err", "gaus", dx_sig_gaus_min, dx_sig_gaus_max);
		tf1_dx_sig_err->SetNpx(int(200.0*(dx_sig_gaus_max - dx_sig_gaus_min) ));
		h_dx_subrange_sig->Fit(tf1_dx_sig_err, "RS0+");



		dx_sig_integral = tf1_dx_sig_err->Integral(dx_sig_gaus_min, dx_sig_gaus_max);
		sigma_dx_sig_integral = tf1_dx_sig_err->IntegralError(dx_sig_gaus_min, dx_sig_gaus_max)/tf1_dx_sig_err->GetHistogram()->GetBinWidth(1);

		sigma_dx_TOTAL = sqrt( pow( sigma_dx_full_integral, 2) + pow( sigma_dx_sig_integral, 2) );


		dx_sig_integral = h_dx_sig->Integral()/h_dx_sig->GetBinWidth(1);

		cout << "dx_sig_integral: " << dx_sig_integral << endl;

	// //------------------------------------------------------------
	// 	tf1_dx_bg_pol8 = new TF1("tf1_dx_bg_pol8", fitPol8, -2.5, 2.5, 9);
	// 	tf1_dx_bg_pol8->SetNpx(1000);
	// 	tf1_dx_bg_pol8->FixParameter(0, par[3]);
	// 	tf1_dx_bg_pol8->FixParameter(1, par[4]);
	// 	tf1_dx_bg_pol8->FixParameter(2, par[5]);
	// 	tf1_dx_bg_pol8->FixParameter(3, par[6]);
	// 	tf1_dx_bg_pol8->FixParameter(4, par[7]);
	// 	tf1_dx_bg_pol8->FixParameter(5, par[8]);
	// 	tf1_dx_bg_pol8->FixParameter(6, par[9]);
	// 	tf1_dx_bg_pol8->FixParameter(7, par[10]);
	// 	tf1_dx_bg_pol8->FixParameter(8, par[11]);

	// 	tf1_dx_bg_pol8->SetLineColor(9);
	// 	h_dx->Fit("tf1_dx_bg_pol8", "R0+");

	// 	tf_dx_sig_integral = h_tf_dx_sig->Integral();
	// 	int dx_minus_bg_cnt = h_dx->GetEntries() - int(tf_dx_sig_integral);
	// 	TH1D *h_dx_bg_rand = new TH1D("h_dx_bg_rand", "dx: Randomly filled BG from tf for BG", 1000, -2.5, 2.5);
	// 	cout << "Filling dx bg_rand with: " << dx_minus_bg_cnt << " events." << endl;
	// 	h_dx_bg_rand->FillRandom("tf1_dx_bg_pol8", dx_minus_bg_cnt);

	// 	h_dx->GetXaxis()->SetRangeUser(-2.5, -0.3);
	// 	// dx_bg_max = h_dx->GetMaximum();
	// 	dx_bg_max = tf1_dx_bg_pol8->GetMaximum();
	// 	h_dx->GetXaxis()->SetRangeUser(-2.5, 2.5);

	// 	h_tf_dx_bg = (TH1*)tf1_dx_bg_pol8->GetHistogram();
	// 	h_tf_dx_bg->SetLineColor(7);
	// 	h_tf_dx_bg->Scale(1.0/h_tf_dx_bg->Integral(), "width");
	// 	h_tf_dx_bg->Scale(1.0/h_tf_dx_bg->Integral(), "height");
	// 	h_tf_dx_bg->Scale( (dx_bg_max)/(h_tf_dx_bg->GetMaximum()) );

	// 	TH1D *h_dx_sig = (TH1D*)h_dx->Clone("h_dx_sig");
	// 	h_dx_sig->SetLineColor(6);
	// 	for( int bin = 0; bin <= 1000; bin++ ){
			
	// 		TAxis *dx_xaxis = h_dx->GetXaxis();
	// 		// double x = bin*0.005;
	// 		double min_x = par[1] - 4*par[2];
	// 		double max_x = par[1] + 4*par[2];

	// 		Int_t min_bin_x = dx_xaxis->FindBin(min_x);
	// 		Int_t max_bin_x = dx_xaxis->FindBin(max_x);

	// 		if( bin >= min_bin_x && bin <= max_bin_x ){
	// 			h_dx_sig->SetBinContent( bin, (h_dx->GetBinContent(bin)) - (h_tf_dx_bg->GetBinContent(bin)) );	
	// 		}
	// 		else{
	// 			h_dx_sig->SetBinContent( bin, 0 );
	// 		}			
	// 	}

	// 	h_dx_sig->Draw("hist+same");
	// 	sig_integral = h_dx_sig->Integral();

		TLegend *tl_dx = new TLegend(0.125, 0.70, 0.425, 0.85,"Full Range dx Fit");
		tl_dx->AddEntry(tf1_dx_fullFit, "Total dx fit");
		tl_dx->AddEntry(h_dx_sig, "dx peak after BG sub");
		tl_dx->AddEntry(tf1_dx_bg, "BG fit");
		// tl_no_reject->AddEntry((TObject*)0, Form("Fit Chi-Sq = %.2f", tf1_dx_fullFit->GetChisquare()), "");
		tl_dx->Draw("same");

		TPaveText *tp_dx = new TPaveText(0.60, 0.55, 0.90, 0.69, "NDCbr");
		tp_dx->AddText(Form("dx total fit chi-sq: %.2f", tf1_dx_fullFit->GetChisquare()));
		tp_dx->AddText(Form("Estimated real dx count: %i", int(dx_sig_integral)));
		// tp_dx->AddText("---------");
		// tp_dx->AddText(Form("BG-sub. W2 cnt = %i", int(bgCut_sig_cnt)));
		// tp_dx->AddText(Form("--> efficiency = %.2f%%", 100.0*dx_sig_integral/bgCut_sig_cnt));
		// tp_dx->AddText("---------");
		// tp_dx->AddText(Form("W2 Gaus fit count: %i", int(bgCut_sig_gaus_cnt)));
		// tp_dx->AddText(Form("--> efficiency: %.2f%%", 100.0*dx_sig_integral/bgCut_sig_gaus_cnt));		
		tp_dx->Draw("same");

		//-----ADD SOME OF THIS INFO TO THE W2_BGCUT PLOT---------
		tp_W2_bgCut->AddText(Form("Estimated real dx count: %i", int(dx_sig_integral)));
		tp_W2_bgCut->AddText("----------------------");
		tp_W2_bgCut->AddText(Form("Bg-sub W2 (max.) efficiency = %.2f%%", 100.0*dx_sig_integral/bgCut_sig_cnt));
		tp_W2_bgCut->AddText("----------------------");
		tp_W2_bgCut->AddText(Form("Gaus-fit W2 (min.) efficiency = %.2f%%", 100.0*dx_sig_integral/bgCut_sig_gaus_cnt));

		cout << "dx_sig_integral before calc: " << dx_sig_integral << endl;
		det_eff_FINAL = dx_sig_integral/bgCut_sig_cnt;

		sigma_det_eff_FINAL = 100.0*det_eff_FINAL*sqrt( pow( (100.0*sigma_dx_TOTAL/dx_sig_integral), 2) + pow( (100.0*sigma_bgCut_TOTAL/bgCut_sig_cnt) , 2) );

		double basic_W2_err = sqrt(bgCut_sig_integral);
		double basic_W2_err_pct = sqrt(bgCut_sig_integral)/bgCut_full_integral;

		double basic_dx_err = sqrt(dx_sig_integral);
		double basic_dx_err_pct = sqrt(dx_sig_integral)/dx_full_integral;

		double basic_TOTAL_err = sqrt( pow( basic_W2_err_pct, 2 ) + pow( basic_dx_err_pct, 2)  );
		double basic_TOTAL_sigma = basic_TOTAL_err*100.0;

		// cout << "Detector efficiency = " << 100.0*det_eff_FINAL << "% +/- " << sigma_det_eff_FINAL << "%" << endl;
		cout << "Detector efficiency = " << 100.0*det_eff_FINAL << "% +/- " << basic_TOTAL_sigma << "%" << endl;

	}
		auto total_time_end = high_resolution_clock::now();
		auto total_time_duration = duration_cast<minutes>(total_time_end - total_time_start);
		cout << "Total time for analysis: " << total_time_duration.count() << " minutes. " << endl;
		cout << endl << "Outfile: " << outfile_name.Data() << endl;
}