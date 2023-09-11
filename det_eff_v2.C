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

double W2_fullFit_par[14], W2_fullFit_par_errors[14], W2_bgSub_polN_par[9], dx_fullFit_par[14], dx_fullFit_par_errors[14], dy_fullFit_par[7];
double W2_bg_reject_par[5];

//Run info and lookups
TString run_target = "LH2";
int kine = 4;
int sbsfieldscale = 0;

//Choose from 3 options:
// "rebuild_only" to only rebuild without fitting/plotting
// "fit_only" to only run fits only
// "rebuild_and_fit" to rebuild then subsequently fit. 

TString rebuild_andOr_plot = "rebuild_and_fit";
int rebuild_andOr_plot_iter;
vector<bool> fit_or_build = {false, true};
bool fit_only;

//Calc Method for kinematic-based variables:
//1) Use four-momentum functions
//2) Use tree kinematic variables (e.kine.**)
//3) Use reconstructed variable functions

int calc_method = 0; //1 isn't good

bool fit_with_expo = false;

int polN = 4;
int interpolN = 99;
int dypolN, dxpolN;
bool use_parsed = true;

double dx_sig_multiplier = 3.5; //SBS4: 3.5 SBS8: 2.75
double dy_sig_multiplier = 2.5; //SBS4: 2.5 SBS8: 2.5
double W2_sig_multiplier = 4.0;

double dx_bin_factor = 200.0;
double dxdy_min_y = -2.5; //this is the dx range
double dxdy_max_y = 2.5;  //y on hcal is dx
int dxdy_y_nbins = int(dx_bin_factor*(dxdy_max_y - dxdy_min_y )); //number of dx bins

double dxdy_min_x = -1.5; //this is the dy range
double dxdy_max_x = 1.5;  //x on hcal is dy
int dxdy_x_nbins = int(dx_bin_factor*(dxdy_max_x - dxdy_min_x )); //number of dy bins

double dx_min_x;
double dx_max_x;
int dx_nbins;

double dy_min_y;
double dy_max_y;
int dy_nbins;

double W2_bin_factor = dx_bin_factor;
double W2_fullScale_min_x = 0.0;
double W2_fullScale_max_x = 3;
int W2_fullScale_nbins = int(W2_bin_factor*(W2_fullScale_max_x - W2_fullScale_min_x));
double W2_min_x = 0.0;
double W2_max_x = 1.4; //SBS4 1.4, SBS8, 1.5 maybe 1.2?
int W2_nbins = int(W2_bin_factor*(W2_max_x - W2_min_x));

double reject_min;
double reject_max;

bool use_bbcal_cuts = false;
bool fiducial_cut = true;
bool correct_beam_energy = true;
bool apply_fcut = false;

bool calc_W_only = false;
bool fit_with_reject = true;

bool acceptance_match = true;

//hcal_cluster_minimize
//Which method to use for "Best cluster selection"
//Choices are: "coin_time", "theta_pq", "dxdy"
TString hcal_cluster_minimize = "theta_pq";

TString hcal_cluster_PID = "proton";
bool use_hcal_cluster_PID = false;
int hcal_cluster_PID_particle = 0; //0: None, 1: Proton, 2: Neutron;
bool PID_clus_match_bool = false;

bool sort_hcal_cluster_energy = true;
bool use_best_cluster = true;
bool use_scoring = true;
int num_best_cluster_types = 1;
//Choose the best cluster type to prioritize:
// 1 --> Highest energy cluster
// 2 --> Cluster to minimize "dx"
// 3 --> Cluster to minimize ADC time different

int priority_best_cluster_type = 3;
bool theta_pq_cut = true;

//best cluster counters:
int e_clus_in_match = 0, ADC_clus_in_match = 0, dx_clus_in_match = 0;
int num_score_0 = 0, num_score_1 = 0, num_score_2 = 0;
int e_ADC_clus_match = 0, e_dx_clus_match = 0, ADC_dx_clus_match = 0;
int e_clus_selections = 0, ADC_clus_selections = 0, dx_clus_selections = 0;
int total_clusters_selected = 0;

bool plot_dxdy = true;
bool plot_dxdy_anticut = true;

//FIT FUNCTIONS FOR DETECTOR EFFICIENCY

//----------
//W2 fullFit (standard) fits
TF1 *tf1_W2_fullFit, *tf1_W2_gaus, *tf1_W2_gaus_err, *tf1_W2_bg_pol4, *tf1_W2_bg_reject;

//W2 fullFit (standard) fits
TH1 *h_tf1_W2_gaus, *h_tf1_W2_bg_pol4, *h_tf1_W2_bg_reject;
TH1D *h_W2_gaus, *h_W2_bg_pol4, *h_W2_bg_pol8;
TH1D *h_W2_bg_from_sub, *h_W2_sig_from_sub, *h_W2_bg_reject;

//W2 pol8 fit to Bg subtracted W2
TH1 *h_tf1_W2_bgSub_polN, *h_tf1_W2_bgSub_gaus;
TH1D *h_W2_bgSub_gaus, *h_W2_bgSub_polN;
//----------

//dy fits
TF1 *tf1_dy_fullFit, *tf1_dy_pol, *tf1_dy_gaus;

//dx fits
TF1 *tf1_dx_fullFit, *tf1_dx_pol, *tf1_dx_gaus, *tf1_dx_bgSub_gaus;

//dx histograms
TH1 *h_tf1_dx_pol;
TH1D *h_dx_gaus, *h_dx_pol, *h_dx_sig_from_sub;

//----------

///INFile histograms
TH1D *h_W2_in, *h_dx_in, *h_dy_in, *h_p_Nucleon_in, *h_p_proton_in, *h_p_neutron_in;
TH2D *h_dxdy_in;

TH1D *h_W2elastics, *h_W2elastics_thr2, *h_W2elastics_fcut, *h_W2elasticsFull_fcut;

double parInterpol[14];

//Integral Counts and values
double W2_gaus_integral, W2_bgSub_pol8_integral, W2_bgSub_fullFit_integral;
double dx_fullFit_integral, dx_pol_integral, dx_gaus_integral;
double p_Nucleon_mean, p_proton_mean, p_neutron_mean;

//Errors for fits and propagation
double W2_bgSub_fullFit_integral_error, W2_gaus_integral_error, W2_bgSub_pol8_integral_error, W2_sig_from_sub_integral, W2_sig_from_sub_integral_error;
double W2_sig_from_sub_par0, W2_sig_from_sub_par1, W2_sig_from_sub_par2;
double W2_sig_from_sub_par0err, W2_sig_from_sub_par1err, W2_sig_from_sub_par2err, W2_sig_from_sub_IntegralError;
double dx_fullFit_integral_error, dx_pol_integral_error, dx_gaus_integral_error, dx_bgSub_gaus_integral_error;
double W2_addsub_error, W2_multdiv_error, W2_total_error;
double dx_addsub_error, dx_multdiv_error, dx_total_error;
double det_eff_total_error, det_eff_FINAL;
double dy_fullFit_integral_error, dy_gaus_integral_error;

//norms and means for dx, dy, etc.
double dy_sig_norm, dy_sig_mean, dy_sig_sigma, dy_sig_min, dy_sig_max, dy_integral, dy_fit_min, dy_fit_max;
double W2_gaus_min_x, W2_gaus_max_x, W2_gaus_min_bin, W2_gaus_max_bin;
double dx_gaus_min_x, dx_gaus_max_x, dx_gaus_min_bin, dx_gaus_max_bin, dx_gaus_norm, dx_gaus_mean, dx_gaus_sigma;
double dx_sig_mean, dx_sig_sigma;

#include "/w/halla-scshelf2102/sbs/jboyd/include/fit_functions.h"


vector<int> runnum_vec;

TString experiment = "gmn";

int pass;
double W2_yield, W2_yield_reject;
double dx_scale = dx_sig_multiplier;
double dy_scale = dy_sig_multiplier;
double W2_scale = W2_sig_multiplier;
double dx_min, dx_max, dy_min, dy_max, W2_min, W2_max;
bool global_cut_fail = false;
bool hit_on_HCal = false;

double x_min_W2 = 0.0;
double x_max_W2 = W2_max_x;
int nBins_W2 = int((x_max_W2 - x_min_W2)*W2_bin_factor);

double W2_fit_min = 0.0;
double W2_fit_max = W2_max_x;
int nBins_fit_W2 = int((W2_fit_max - W2_fit_min)*W2_bin_factor);

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
double hcal_y_fmin, hcal_y_fmax, hcal_x_fmin, hcal_x_fmax, ADC_time_min, ADC_time_max, ADC_time_mean;

//Scattered kinematics
double e_prime_theta; //Scattered electron theta angle
double e_prime_phi; //Scattered electron phi angle
double p_el, nu, p_Nucleon, nucleon_theta, nucleon_phi, E_ep, p_ep, Q2, W, W2, E_pp, E_nucleon, KE_p, dx, dy;

TLorentzVector pN; //Scattered Nucleon momentum
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
const Int_t max_clus = 100;
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
Double_t hcal_clus_e[max_clus], hcal_clus_x[max_clus], hcal_clus_y[max_clus], hcal_clus_atime[max_clus], hcal_clus_tdctime[max_clus];
Array1DValueWithIndex hcal_clus_e_sorted[maxTracks];
// Array1DScoredWithIndex hcal_clus_scored[maxTracks];
vector<vector<int>> hcal_clus_scored;
//{ SCORE, HIGHEST_E_CLUS_INDEX, BEST_TIMING_INDEX, BEST_DXDY_INDEX, BEST_THETA_PQ_INDEX}
vector<int> hcal_clus_row = {0, -1, -1, -1, -1};

Double_t dx_bestcluster, dy_bestcluster, HCal_ADC_time_bestcluster;
Double_t dx_cluster_final, dy_cluster_final, HCal_ADC_time_final;
Double_t hcal_clus_id, hcal_clus_nblk;

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
Int_t TDCTndata, Nhcal_clus_id;

Long64_t nTCevents, Nevents;

//INITIALIZE ALL HISTOGRAMS:
TH1D *hin_bb_gem_Ntracks, *h_W2_full_rand;
TH1D *h_atime, *h_W, *h_W2, *h_W2copy, *h_W2interpol, *h_W2recon, *h_KE_p, *h_KE_low, *h_Diff, *h_X, *h_Y, *h_E_eloss, *h_hcal_clusblk_ADC_time;
TH1D *h_Wfull, *h_W2full, *h_W2_residual;
TH1D *h_W_cut, *h_W_fcut, *h_vz_cut, *h_theta_pq_p, *h_theta_pq_n, *h_theta_pq_p_pNhat, *h_theta_pq_n_pNhat, *h_hcal_e;
TH1D *h_Q2, *h_E_ep, *h_E_pp, *h_Ep, *h_p_Nucleon, *h_p_proton, *h_p_neutron;
TH1D *h_dy, *h_dy_cut, *h_dy_wcut, *h_dx, *h_dx_cut, *h_dx_wcut, *h_dx_fcut, *h_dx_wcut_fcut, *h_dy_wcut_fcut;
TH1D *h_dx_w2cut, *h_dx_w2cut_fcut;
TH1D *h_dy_subrange, *h_dx_subrange;
TH1D *h_Nevents;

TH1D *h_rej_bg, *h_rej_sig_gaus;

//Interpolation "tight" histograms
TH1D *h_W2antiCut, *h_W2fullCut, *h_W2elasticsFull, *h_W2antiCutFull, *h_dx_allCuts;
TH1D *h_W2antiCut_fcut, *h_W2fullCut_fcut, *h_W2antiCutFull_fcut, *h_dx_allCuts_fcut;

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

int polNfit, interpolNfit;
int min_max_iter_cnt;
int hcal_missed_cnt = 0;
int hcal_hit_cnt = 0;
int Nclusters = 0;

int badHcalESort = 0, badHcalEsortElastics = 0, numScore3 = 0, clusterWithScore3 = -1;
bool foundclusterWithScore3 = false;

double defficiency;

void det_eff_v2(){

	auto total_time_start = chrono::high_resolution_clock::now();
	TStopwatch *StopWatch = new TStopwatch();

	gStyle->SetPalette(55);
	cout << "--------------------------------------" << endl;
	cout << "Analysis started. " << endl;
	cout << "--------------------------------------" << endl;
	if( use_parsed == true ){
		cout << "----------Using Parsed Rootfiles----------" << endl;
	}
	if( use_parsed == false){
		cout << "----------Using NON-Parsed Rootfiles----------" << endl;
	}
	cout << "--------------------------------------" << endl;

	cout << "Run parameters: " << endl;
	cout << "Kinematic: SBS" << kine << endl;
	cout << "Target: " << run_target.Data() << endl;
	cout << "Beam Energy: " << E_beam << endl;
	cout << "SBS Field: " << sbsfieldscale << "%" << endl;
	cout << "-----------------------------------" << endl;
	cout << "BB angle [deg]: " << lookup_BB_angle_by_kine(kine, "deg") << endl;
	cout << "SBS angle: " << lookup_SBS_angle_by_kine(kine, "deg") << endl;
	cout << "HCal angle [deg]: " << lookup_HCal_angle_by_kine(kine, "deg") << endl;
	cout << "HCal distance: " << HCal_dist << endl;
	cout << "-----------------------------------" << endl << endl;

	if( run_target == "LH2" ){
		Nucleon_mass = Mp;
	}
	if( run_target == "LD2" ){
		Nucleon_mass = 0.5*(Mp + Mn);
	}

	if( use_hcal_cluster_PID && hcal_cluster_PID == "proton" ){
		hcal_cluster_PID_particle = 1; 
	}
	else if( use_hcal_cluster_PID && hcal_cluster_PID == "neutron" ){
		hcal_cluster_PID_particle = 2; 
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

		dx_min_x = -0.5;
		dx_max_x = 0.5;

		dy_min_y = -0.4;
		dy_max_y = 0.4;

	}

	if( kine == 8 ){
		pass = 1;
		hcalheight = -0.450; // Height of HCal above beamline

		W_mean = 8.81088e-01;
		W_sigma = 9.95799e-02;
		// W2_mean = 7.28445e-01;
		// W2_sigma = 1.16788e-01;
		
		W2_mean = 0.72903487;
		W2_sigma = 0.12343966;

		//offset: -1.37521e-01
		dx_min = -6.98838e-02 - (dx_scale)*7.28711e-02;
		dx_max = -6.98838e-02 + (dx_scale)*7.28711e-02;
		dx_p = -0.050000000;
		dx_p_sigma = 0.080900537;
		dx_n = dx_p;
		dx_n_sigma = dx_p_sigma;

		dy_min = -0.0374881 - (dy_scale)*0.056;
		dy_min = -0.0374881 + (dy_scale)*0.056;

		dy_p = -0.0374881;
		dy_p_sigma = 0.056;
		dy_n = dy_p;
		dy_n_sigma = dy_p_sigma;

		dx_pn_max = 0;

		dx_min_x = -1.0;
		dx_max_x = 0.5;

		dy_min_y = -0.4;
		dy_max_y = 0.2;

	}

	W2_min = W2_mean - (W2_scale)*W2_sigma;
	W2_max = W2_mean + (W2_scale)*W2_sigma;

	dx_nbins = int(dx_bin_factor*(dx_max_x - dx_min_x));
	dy_nbins = int(dx_bin_factor*(dy_max_y - dy_min_y));

	TString dy_sig_mult_str = Form("%.1f", dy_sig_multiplier);
	dy_sig_mult_str.ReplaceAll(".", "_");

	TString W2_sig_mult_str = Form("%.1f", W2_sig_multiplier);
	W2_sig_mult_str.ReplaceAll(".", "_");
	
	outfile_name = Form("rootfiles/det_eff_histos_SBS%i_%s_mag%i_bgSub_W2BinFactor%i_dxBinFactor%i_dySigMult%s_w2SigMult%s_scoreTest_%i.root", kine, run_target.Data(), sbsfieldscale, int(W2_bin_factor), int(dx_bin_factor), dy_sig_mult_str.Data(), W2_sig_mult_str.Data(), priority_best_cluster_type);

	cout << "Outfile name: " << outfile_name.Data() << endl;

	if( rebuild_andOr_plot == "rebuild_and_fit"){
		cout << "---------------------------------" << endl << endl;
		cout << "		Running in Rebuild AND Fit Mode " << endl;
		cout << "---------------------------------" << endl << endl;
		rebuild_andOr_plot_iter = 2;

	}
	if( rebuild_andOr_plot == "rebuild_only" ){
		cout << "---------------------------------" << endl << endl;
		cout << "		Running in Rebuild Only Mode " << endl;
		cout << "---------------------------------" << endl << endl;

		rebuild_andOr_plot_iter = 1;
		fit_only = false;
	}
	if( rebuild_andOr_plot == "fit_only" ){
		cout << "---------------------------------" << endl << endl;
		cout << "		Running in Fit Only Mode " << endl;
		cout << "---------------------------------" << endl << endl;

		rebuild_andOr_plot_iter = 1;
		fit_only = true;
	}

	for( int rebuild_plot_sel = 0; rebuild_plot_sel < rebuild_andOr_plot_iter; rebuild_plot_sel++ ){
	
		if( rebuild_andOr_plot == "rebuild_and_fit"){
			fit_only = fit_or_build[rebuild_plot_sel];
			cout << "---------------------------------" << endl << endl;
			if( fit_only ){
				cout << "		Running in Rebuild AND Fit Mode " << endl;
				cout << "             On Fit Only portion" << endl;
			}
			if( !fit_only ){
				cout << "		Running in Rebuild AND Fit Mode " << endl;
				cout << "             On Rebuild Only portion" << endl;
			}			
			cout << "---------------------------------" << endl << endl;
		}

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
		ADC_time_mean = lookup_ADC_time_cut(run_target, kine, sbsfieldscale, "ADC_time_mean");

		if( use_parsed ){
			cout << endl << "-------- Parsed run mode --------" << endl;
			// rootfile_dir = "/lustre19/expphy/volatile/halla/sbs/jboyd/analysis_rootfiles/jboyd_parsed";
			rootfile_dir = Form("/volatile/halla/sbs/jboyd/analysis_rootfiles/jboyd_parsed/SBS%i/%s/mag%i", kine, run_target.Data(), sbsfieldscale);
			infile_name = Form("gmn_parsed_%s_SBS%i_mag%i_binFactor.root", run_target.Data(), kine, sbsfieldscale);
			cout << endl << "Adding files to TChain from..." << endl << endl;
			cout << rootfile_dir.Data() << endl << endl;
			cout << "                ------- " << endl;
			// TC->Add(Form("%s/%s", rootfile_dir.Data(), infile_name.Data()));
			TC->Add(Form("%s/*", rootfile_dir.Data()));
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

		// HClas Cluster Tree Variables
		TC->SetBranchStatus( "sbs.hcal.clus.e", 1);
		TC->SetBranchStatus( "sbs.hcal.clus.x", 1);
		TC->SetBranchStatus( "sbs.hcal.clus.y", 1);
		TC->SetBranchStatus( "sbs.hcal.clus.atime", 1);
		TC->SetBranchStatus( "sbs.hcal.clus.tdctime", 1);
		TC->SetBranchStatus( "sbs.hcal.clus.id", 1);
		TC->SetBranchStatus( "sbs.hcal.clus.nblk", 1);

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
		TC->SetBranchStatus( "Ndata.sbs.hcal.clus.id", 1 );

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
		TC->SetBranchAddress( "Ndata.sbs.hcal.clus.id", &Nhcal_clus_id );

		// HCal Cluster Tree Variables
		TC->SetBranchAddress( "sbs.hcal.clus.e", &hcal_clus_e);
		TC->SetBranchAddress( "sbs.hcal.clus.x", &hcal_clus_x);
		TC->SetBranchAddress( "sbs.hcal.clus.y", &hcal_clus_y);
		TC->SetBranchAddress( "sbs.hcal.clus.atime", &hcal_clus_atime);
		TC->SetBranchAddress( "sbs.hcal.clus.tdctime", &hcal_clus_tdctime);
		TC->SetBranchAddress( "sbs.hcal.clus.id", &hcal_clus_id);
		TC->SetBranchAddress( "sbs.hcal.clus.nblk", &hcal_clus_nblk);

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

		h_Ep = new TH1D("h_Ep", "Total E over P", 400, 0.0, 2.0);
		h_p_Nucleon = new TH1D("h_p_Nucleon", Form("Scattered Nucleon momentum (pN) - SBS%i %i, %s", kine, sbsfieldscale, run_target.Data()), 1200, 1.0, 7.0);
		h_p_proton = new TH1D("h_p_proton", Form("Scattered proton momentum (pp) - SBS%i %i, %s", kine, sbsfieldscale, run_target.Data()), 1200, 1.0, 7.0);
		h_p_neutron = new TH1D("h_p_neutron", Form("Scattered neutron momentum (pn) - SBS%i %i, %s", kine, sbsfieldscale, run_target.Data()), 1200, 1.0, 7.0);		
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

		h_dx = new TH1D("h_dx",Form("dx (NO CUTS) - SBS%i %i, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), dxdy_y_nbins, dxdy_min_y, dxdy_max_y);

		h_dx_cut = new TH1D("h_dx_cut",Form("dx (Basic CUTS) - SBS%i %i, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), dxdy_y_nbins, dxdy_min_y, dxdy_max_y);
		h_dx_wcut = new TH1D("h_dx_wcut",Form("dx (W cut) - SBS%i %i, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), dxdy_y_nbins, dxdy_min_y, dxdy_max_y);
		h_dx_w2cut = new TH1D("h_dx_w2cut",Form("dx (W^2 cut) - SBS%i %i, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), dxdy_y_nbins, dxdy_min_y, dxdy_max_y);
		h_dx_fcut = new TH1D("h_dx_fcut",Form("dx (f cut) - SBS%i %i, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), dxdy_y_nbins, dxdy_min_y, dxdy_max_y);
		h_dx_wcut_fcut = new TH1D("h_dx_wcut_fcut",Form("dx (W & Fiduc. Cuts) - SBS%i %i, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), dxdy_y_nbins, dxdy_min_y, dxdy_max_y);
		h_dx_w2cut_fcut = new TH1D("h_dx_w2cut_fcut",Form("dx (W^2 & Fiduc. Cuts) - SBS%i %i, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), dxdy_y_nbins, dxdy_min_y, dxdy_max_y);
		
		h_dy = new TH1D("h_dy",Form("dy (NO CUTS) - SBS%i %i, %s; y_{HCal} - y_{exp} (m);", kine, sbsfieldscale, run_target.Data()), dxdy_x_nbins, dxdy_min_x, dxdy_max_x);
		h_dy_cut = new TH1D("h_dy_cut",Form("dy (Basic Cuts) - SBS%i %i, %s; y_{HCal} - y_{exp} (m);", kine, sbsfieldscale, run_target.Data()), dxdy_x_nbins, dxdy_min_x, dxdy_max_x);  
		h_dy_wcut = new TH1D("h_dy_wcut",Form("dy (W Cuts) - SBS%i %i, %s; y_{HCal} - y_{exp} (m);", kine, sbsfieldscale, run_target.Data()), dxdy_x_nbins, dxdy_min_x, dxdy_max_x);
		h_dy_wcut_fcut = new TH1D("h_dy_wcut_fcut",Form("dy (W & Fiduc. Cuts) - SBS%i %i, %s; y_{HCal} - y_{exp} (m);", kine, sbsfieldscale, run_target.Data()), dxdy_x_nbins, dxdy_min_x, dxdy_max_x);    

		h_dxdy = new TH2D("h_dxdy", Form("Hadron Spot(s) on HCal (NO CUTS) - SBS%i %i, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), dxdy_x_nbins, dxdy_min_x, dxdy_max_x, dxdy_y_nbins, dxdy_min_y, dxdy_max_y );
		h_dxdy_wcut = new TH2D("h_dxdy_wcut", Form("Hadron Spot(s) on HCal (W cut) - SBS%i %i, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), dxdy_x_nbins, dxdy_min_x, dxdy_max_x, dxdy_y_nbins, dxdy_min_y, dxdy_max_y );
		h_dxdy_wcut_fcut = new TH2D("h_dxdy_wcut_fcut", Form("Hadron Spot(s) on HCal (W & Fiduc. Cuts) - SBS%i %i, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), dxdy_x_nbins, dxdy_min_x, dxdy_max_x, dxdy_y_nbins, dxdy_min_y, dxdy_max_y );
		h_dxdy_cut = new TH2D("h_dxdy_cut", Form("Hadron Spot(s) on HCal (Basic cuts) - SBS%i %i, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), dxdy_x_nbins, dxdy_min_x, dxdy_max_x, dxdy_y_nbins, dxdy_min_y, dxdy_max_y );
		h_dxdy_ncut = new TH2D("h_dxdy_ncut", Form("Hadron Spot(s) on HCal (n cut) - SBS%i %i, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), 250, -1.25, 1.25, 250, -2.5, 2.5 );
		h_dxdy_pcut = new TH2D("h_dxdy_pcut", Form("Hadron Spot(s) on HCal (p cut) - SBS%i %i, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), 250, -1.25, 1.25, 250, -2.5, 2.5 );
		h_dxdy_fcut = new TH2D("h_dxdy_fcut", Form("Hadron Spot(s) on HCal (f cut) - SBS%i %i, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), dxdy_x_nbins, dxdy_min_x, dxdy_max_x, dxdy_y_nbins, dxdy_min_y, dxdy_max_y );
		
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

		h_hcal_e = new TH1D("h_hcal_e", Form("HCal cluster energy - SBS%i %i, %s; theta_pq_n (rad);", kine, sbsfieldscale, run_target.Data()), 120, 0.0, 0.6);

		cout << "--------------------------------------" << endl;
		cout << "Number of events to analyze: " << Nevents << endl;
		cout << "--------------------------------------" << endl;
		cout << "--------------------------------------" << endl;
		cout << "Starting analysis loop on events..... " << endl;

	//Basic Energy calcs

		E_loss_outgoing = cell_diameter/2.0/sin(BB_theta)*rho_tgt*dEdx_tgt; //Should be about 1 MeV
		if( useAlshield !=0 ) E_loss_outgoing += Alshieldthick*rho_Al*dEdx_Al;

	//FIDUCIAL CUT:
		//dimensions for HCal edges:
		//x_i = -2.268095, x_f = 1.538095..... y_i = -0.931545, y_f = 0.931545
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

			if( sort_hcal_cluster_energy ){
			//SORT THE HCAL CLUSTER ENERGY ARRAY BY DESCENDING ORDER
				for(int clus = 0; clus < max_clus; clus++){
					hcal_clus_e_sorted[clus].ArrayValue = hcal_clus_e[clus];
					hcal_clus_e_sorted[clus].index = clus;
				}

				sort( hcal_clus_e_sorted, hcal_clus_e_sorted + max_clus, CompareArrayValueWithIndex );

				hcal_x = hcal_clus_x[hcal_clus_e_sorted[0].index];
				hcal_y = hcal_clus_y[hcal_clus_e_sorted[0].index];			
			}


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
			h_hcal_clusblk_ADC_time->Fill(hcal_clusblk_ADC_time[0]);


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

			h_hcal_e->Fill(hcal_e);

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
				p_Nucleon = sqrt( pow(nu, 2.0) + 2.0*(Nucleon_mass*nu) );
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
			// cout << "p_Nucleon = " << p_Nucleon << endl;
			h_p_Nucleon->Fill( p_Nucleon );
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
			TVector3 qvec = q.Vect();

		//Expected neutron direction
			TVector3 Neutron_Direction = (HCAL_pos - vertex).Unit();

		//Expected proton direction
			//Need to incorporate deflection due to SBS magnet
			double Bdl = sbsfieldscale*maxSBSfield*Dgap/100.0;
			double Proton_Deflection = tan( 0.3*Bdl/q.Mag() )*(HCal_dist - (SBSdist + Dgap/2.0) ); 

			TVector3 Proton_Direction = (HCAL_pos + Proton_Deflection*HCAL_xaxis - vertex).Unit();

			// theta_pq_n = acos( Neutron_Direction.Dot( qvec_recon.Unit() ) );
			// theta_pq_p = acos( Proton_Direction.Dot( qvec_recon.Unit() ) );

		// THETA_PQ
			theta_pq_p = acos( Proton_Direction.Dot( pNhat ));
			theta_pq_n = acos( Neutron_Direction.Dot( pNhat ));

			h_theta_pq_n->Fill(theta_pq_n);
			h_theta_pq_p->Fill(theta_pq_p);

	//-----------------------------------------------------
	// Search for best cluster on HCal (Not just highest energy cluster)
	//-----------------------------------------------------

			Int_t clus_sel_dx_atime = 0;
			Double_t clus_sel_atimediff = 1000.0;

			Int_t clus_sel_dx_theta_pq = 0;
			Double_t clus_sel_theta_pq_diff = 1000.0;

			Int_t clus_sel_dx_dxdy = 0;
			Double_t clus_sel_dxdy_diff = 1000.0;

			Double_t clus_sel_theta_pq;

			//Add number of elastic clusters to: Nclusters
			Nclusters+=Nhcal_clus_id;
			// cout << "Nhcal_clus_id to search through: " << Nhcal_clus_id << endl;

			for( Int_t clus_sel = 0; clus_sel < Nhcal_clus_id; clus_sel++){
			// for( Int_t clus_sel = 0; clus_sel < nclus; clus_sel++){
				Double_t clus_sel_dx = hcal_clus_x[clus_sel] - x_expected_HCal;
				Double_t clus_sel_dy = hcal_clus_y[clus_sel] - y_expected_HCal;

				TVector3 clus_sel_HCal_pos = HCAL_origin + hcal_clus_x[clus_sel]*HCAL_xaxis + hcal_clus_y[clus_sel]*HCAL_yaxis;
				TVector3 clus_sel_proton_direction = ( clus_sel_HCal_pos + Proton_Deflection*HCAL_xaxis - vertex ).Unit();
				TVector3 clus_sel_neutron_direction = (clus_sel_HCal_pos - vertex ).Unit();

				Double_t clus_sel_theta_pq_p = acos( clus_sel_proton_direction.Dot( qvec.Unit() ) );
				Double_t clus_sel_theta_pq_n = acos( clus_sel_neutron_direction.Dot( qvec.Unit() ) );


			// Best cluster particle ID

				if( use_hcal_cluster_PID ){
					PID_clus_match_bool = false;
				}
				if( !use_hcal_cluster_PID ){
					PID_clus_match_bool = true;
				}

				Int_t clus_sel_PID; //0 is NEITHER, 1 is PROTON, 2 is NEUTRON

				bool clus_sel_PID_p = (abs(clus_sel_dx - dx_p)<(dx_sig_multiplier*dx_p_sigma)) && (abs(clus_sel_dy - dy_p )<dx_sig_multiplier*dy_p_sigma);
				bool clus_sel_PID_n = (abs(clus_sel_dx - dx_n)<(dx_sig_multiplier*dx_n_sigma)) && (abs(clus_sel_dy - dy_p )<dx_sig_multiplier*dy_p_sigma);

				if( clus_sel_PID_p && !clus_sel_PID_n){
					clus_sel_PID = 1;

					if( use_hcal_cluster_PID && hcal_cluster_PID == "proton" ){
						hcal_cluster_PID_particle = 1;
						PID_clus_match_bool = true;
					}
				}
				else if( clus_sel_PID_n && !clus_sel_PID_p){
					clus_sel_PID = 2;

					if( use_hcal_cluster_PID && hcal_cluster_PID == "neutron" ){
						hcal_cluster_PID_particle = 2;
						PID_clus_match_bool = true;
					}
				}
				else if( clus_sel_PID_p && clus_sel_PID_n){
					clus_sel_PID = 0;

					if( use_hcal_cluster_PID && ( hcal_cluster_PID == "proton" || hcal_cluster_PID == "neutron" ) ){
						hcal_cluster_PID_particle = 0;
						PID_clus_match_bool = false;
					}
				}
				else{
					clus_sel_PID = -1;
				}

			//After PID we can pull indices of best clusters
			//We want cluster which minimized difference between coincidence times
				if( abs(hcal_clus_atime[clus_sel] - ADC_time_mean) < clus_sel_atimediff ){
					clus_sel_dx_atime = clus_sel;
					clus_sel_atimediff = hcal_clus_atime[clus_sel];
				}

			//Cluster to minize theta pq for proton
			//Can select neutron. Should study effects of selection/bias	
				if( clus_sel_theta_pq < clus_sel_theta_pq_diff ){
					clus_sel_dx_theta_pq = clus_sel;
					clus_sel_theta_pq_diff = clus_sel_theta_pq_p;
				}

			//Cluster to minize distance btween expected dx/dy lcoation and center of cluster
				Double_t dxdy_dist = sqrt( pow(clus_sel_dx - dx_p, 2) + pow(clus_sel_dy - dy_p, 2) );
				if( dxdy_dist < clus_sel_dxdy_diff ){
					clus_sel_dx_dxdy = clus_sel;
					clus_sel_dxdy_diff = dxdy_dist;
					clus_sel_theta_pq = clus_sel_theta_pq_p; //Saving proton here
				}

			}

		//Now we can apply this and use it to select "best clusters" to use:
			if( hcal_cluster_minimize == "coin_time" || ( PID_clus_match_bool ) ){
				dx_bestcluster = hcal_clus_x[clus_sel_dx_atime] - x_expected_HCal;
				dy_bestcluster = hcal_clus_y[clus_sel_dx_atime] - y_expected_HCal;
				HCal_ADC_time_bestcluster = hcal_clus_atime[clus_sel_dx_atime];				
			}
			if( hcal_cluster_minimize == "theta_pq" ){
				dx_bestcluster = hcal_clus_x[clus_sel_dx_theta_pq] - x_expected_HCal;
				dy_bestcluster = hcal_clus_y[clus_sel_dx_theta_pq] - y_expected_HCal;
				HCal_ADC_time_bestcluster = hcal_clus_atime[clus_sel_dx_theta_pq];				
			}
			if( hcal_cluster_minimize == "dxdy" ){
				dx_bestcluster = hcal_clus_x[clus_sel_dx_dxdy] - x_expected_HCal;
				dy_bestcluster = hcal_clus_y[clus_sel_dx_dxdy] - y_expected_HCal;
				HCal_ADC_time_bestcluster = hcal_clus_atime[clus_sel_dx_dxdy];				
			}


			//Count bad hcal clus e sorts:
			if( FindIndexOfMaxElement( hcal_clus_e, max_clus ) != 0 ){
				badHcalESort++;
			}
	//////////////////////////////////////
	//////////////////////////////////////
		//Start of best cluster scoring
			hcal_clus_scored.push_back(hcal_clus_row);

			// //DECLARE CLUSTER VALUE INFO IN THE SCORED ARRAY:
			//Index 0 is the score of the cluster:
			//[0] = 2 --> Three matching clusters
			//[0] = 1 --> Two matching clusters
			//[0] = 0 --> All clusters are different
			
			hcal_clus_scored[nevent][1] = FindIndexOfMaxElement( hcal_clus_e, max_clus );
			hcal_clus_scored[nevent][2] = clus_sel_dx_dxdy;
			hcal_clus_scored[nevent][3] = clus_sel_dx_atime;

		//CASES FOR VARIETIES OF MATCHES ON THE INDICES. 
		//------------------------------------
		//Loop through the hcal_cluse_scored array and count the matches
			int maxMatchCount = 0;
			int currentMatchCount = 0;
			for( int i = 1; i <= num_best_cluster_types; i ++ ){
				// Reset the current matching count for each new index
				currentMatchCount = 0;

				// Compare the element at the current index with elements at indices [1] through [last]
				for( int j = 1; j <= num_best_cluster_types; j++ ){
					if( (hcal_clus_scored[nevent][i] == hcal_clus_scored[nevent][j]) && (i != j) ){
						currentMatchCount++;
					}
				}
				// Update the maximum matching count if necessary
		        if (currentMatchCount > maxMatchCount) {
		            maxMatchCount = currentMatchCount; 
		        }
			}

		//Put the maximum number of matches into element [0]
			//minimum number of matches should be 0 and max should be num_best_cluster_types - 1.
			hcal_clus_scored[nevent][0] = maxMatchCount; 

// --------------------------------------

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

			h_Ep->Fill( (bb_sh_e+bb_ps_e)/(bb_tr_p[0]) );

			//Get invariant mass transfer W from the four-momentum of the scattered nucleon
		
			W = pow(W2, 0.5);
			// W = P_gammaN.M();
			// W2 = pow(W, 2);				

		//IS NUCLEON ON HCAL FACE?
			if( acceptance_match ){
				hit_on_HCal = false;
				if( hcal_y>hcal_y_fmin && hcal_y<hcal_y_fmax && hcal_x>hcal_x_fmin && hcal_x<hcal_x_fmax ){
					hit_on_HCal = true;
					hcal_hit_cnt++;
				}
			}
			if( !acceptance_match ){
				hit_on_HCal = true;
			}

			if( hit_on_HCal == false ){
				hcal_missed_cnt++;
				continue;
			}

		//fill work histograms

			h_W->Fill( W );
			h_Wfull->Fill( W );

	//For detection efficiency we should have some general cuts on 
			int bbcal_min_nhits = 3;

			if( false ){
				if( hit_on_HCal && (nclus>0) && (PS_nclus>0) && (SH_nclus>0) && (abs(bb_tr_vz[0])<=0.075) && (bb_nhits[0]>bbcal_min_nhits) && (bb_tr_n==1) && (bb_ps_e>0.15) && (hcal_e>0.025) ){
					h_W2->Fill( W2 );	
				}				
			}
			// bb_tgt_th[0]>-0.2 && 
			bool det_eff_cut;

			if( kine == 4 ){
				det_eff_cut =  	
								hit_on_HCal && 
								SH_nclus>0 && 
								PS_nclus>0 && (abs(bb_tr_vz[0])<=0.075) 
								&& (bb_ps_e>0.150) 
								&& (bb_nhits[0]>3) 
								&& (bb_tr_n==1);				
			}

			if( kine == 8 ){

				det_eff_cut =  	
								hit_on_HCal && 
								SH_nclus>0 && 
								PS_nclus>0 && (abs(bb_tr_vz[0])<=0.075) 
								&& (bb_ps_e>0.150) 
								&& (bb_nhits[0]>4) 
								&& (bb_tr_n==1);
			}
			
			if( det_eff_cut ){
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

				if( use_best_cluster ){
					dx = dx_bestcluster;
					dy = dy_bestcluster;
				}
				if( !use_best_cluster ){
					dx = hcal_x - x_expected_HCal;
					dy = hcal_y - y_expected_HCal;				
				}

			//If using the scoring method for best cluster we need to sort through the scores and matches and pick a best cluster
			//This is all based on the original priority_best_cluster.
			//I think we should work chronogically from score 0 to score 3
				if( use_scoring ){
					// cout << "Using score3 cluster " << numScore3 << endl;

				//CASE 1: No matches. If no elements match we should default to the priority_best_cluster_type selection	
					if( hcal_clus_scored[nevent][0] == 0 ){
						dx = hcal_clus_x[hcal_clus_scored[nevent][priority_best_cluster_type]] - x_expected_HCal;
						dy = hcal_clus_y[hcal_clus_scored[nevent][priority_best_cluster_type]] - y_expected_HCal;

						num_score_0++;
						if( priority_best_cluster_type == 1 ){e_clus_selections++;}
						if( priority_best_cluster_type == 2 ){dx_clus_selections++;}
						if( priority_best_cluster_type == 3 ){ADC_clus_selections++;}
						total_clusters_selected++;
					}

			//These next cases is a bit more interesting. We want to look and see if the priority_best_cluster_type is part of a match.
			//IF THE priority_best_cluster_type IS NOT PART OF THE MAX MATCH THEN WE DON'T WANT IT. 
			//IF THE priority_best_cluster_type IS PART OF THE MAX MATCH THEN WE USE IT. 
					
				//CASE 2: ONE MATCHING PAIR
					if( hcal_clus_scored[nevent][0] == 1 ){
						bool priority_in_match = false;
						vector<int> matching_indices = {};
						//Lets default the cluster choice to the originally chosen priority_best_cluster_type
						int one_match_best_cluster = priority_best_cluster_type;

						for( int i = 1; i <= num_best_cluster_types; i++ ){
							if( (hcal_clus_scored[nevent][i] == hcal_clus_scored[nevent][priority_best_cluster_type]) && (i !=  priority_best_cluster_type) ){
								priority_in_match = true;
								//Priority cluster is in a match so we will pick it!
							}
							//NOW, the priority cluster is not part of a match so.... we pick the second best option. 
							//Sort through and find the indices of the matching pair
							matching_indices = {};
							for( int j = 1; j <= num_best_cluster_types; j++ ){
								if( hcal_clus_scored[nevent][i] == hcal_clus_scored[nevent][j] && (i != j) ){
									matching_indices.push_back(i);
								}
							}
							
							//Counting stuff for metrics
							if( (hcal_clus_scored[nevent][i] == hcal_clus_scored[nevent][1]) && (i !=  1) ){
								e_clus_in_match++;
							}
							if( (hcal_clus_scored[nevent][i] == hcal_clus_scored[nevent][2]) && (i !=  2) ){
								dx_clus_in_match++;
							}
							if( (hcal_clus_scored[nevent][i] == hcal_clus_scored[nevent][3]) && (i !=  3) ){
								ADC_clus_in_match++;
							}

						}
						//Counting stuff for metrics
						if( (hcal_clus_scored[nevent][1]) == hcal_clus_scored[nevent][2] ){
							e_dx_clus_match++;
						}
						if( (hcal_clus_scored[nevent][1]) == hcal_clus_scored[nevent][3] ){
							e_ADC_clus_match++;
						}
						if( (hcal_clus_scored[nevent][2]) == hcal_clus_scored[nevent][3] ){
							ADC_dx_clus_match++;
						}
						total_clusters_selected++;

						if( !priority_in_match ){
							//These are ordered by importance. So, matching_indices[0] will be the next best choice after the priority_best_cluster_type
							one_match_best_cluster = matching_indices[0];
						}
						if( priority_in_match ){
							one_match_best_cluster = priority_best_cluster_type; //A Bit redundant but that's fine
						}
					//The best cluster should have been sorted out by now. 
						dx = hcal_clus_x[hcal_clus_scored[nevent][one_match_best_cluster]] - x_expected_HCal;
						dy = hcal_clus_y[hcal_clus_scored[nevent][one_match_best_cluster]] - y_expected_HCal;						

						num_score_1++;

						if( one_match_best_cluster == 1 ){e_clus_selections++;}
						if( one_match_best_cluster == 2 ){dx_clus_selections++;}
						if( one_match_best_cluster == 3 ){ADC_clus_selections++;}
					}

				//CASE 3: All element match: Select index of priority_best_cluster_type. Doesn't matter, but let us be diligent and select it. 
					if( hcal_clus_scored[nevent][0] == 2 ){
						dx = hcal_clus_x[hcal_clus_scored[nevent][priority_best_cluster_type]] - x_expected_HCal;
						dy = hcal_clus_y[hcal_clus_scored[nevent][priority_best_cluster_type]] - y_expected_HCal;

						//Counting stuff for metrics
						num_score_2++;
						e_clus_in_match++;
						ADC_clus_in_match++;
						dx_clus_in_match++;
						e_ADC_clus_match++;
						e_dx_clus_match++;
						ADC_dx_clus_match++;
						total_clusters_selected++;

						if( priority_best_cluster_type == 1 ){e_clus_selections++;}
						if( priority_best_cluster_type == 2 ){dx_clus_selections++;}
						if( priority_best_cluster_type == 3 ){ADC_clus_selections++;}
					}

					// dx = dx_bestcluster;
					// dy = dy_bestcluster;					
				}
				//Resolve the hadron spots without cuts
				h_dx->Fill( dx );
				h_dy->Fill( dy );
				
				if( det_eff_cut && (nclus>0) ){ //
				// if( hit_on_HCal  && (nclus>0)){ //&& (hcal_e>0.025)
					h_dxdy->Fill( dy, dx );					
				}
				
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
				// if( (theta_pq_p > theta_pq_p_thresh) && (dx<dx_min || dx>dx_max) && (dy<dy_min || dy>dy_max) && (global_cut_fail) && hit_on_HCal ){
				// 	h_W2antiCut->Fill( W2 );
				// 	h_W2antiCutFull->Fill( W2 );
				// }

				// //Looking to cut out all the background and only focus in on the "real elastics" signal here
				// if( (theta_pq_p < theta_pq_p_thresh) && (dx>dx_min || dx<dx_max) && (dy>dy_min || dy<dy_max) && (!global_cut_fail) && hit_on_HCal ){

				// 	h_W2elastics->Fill( W2 );
				// 	h_W2elasticsFull->Fill( W2 );
				// }
				// if( (theta_pq_p < (0.30*theta_pq_p_thresh)) && (dx>dx_min || dx<dx_max) && (dy>dy_min || dy<dy_max) && (!global_cut_fail) && hit_on_HCal ){

				// 	h_W2elastics_thr2->Fill( W2 );

				// }
				if( (W2>W2_min) && (W2<W2_max) && (dy>dy_min) && (dy<dy_max) && (!global_cut_fail) && hit_on_HCal ){
					cout << "filling h_dx_allCuts... " << endl;
					h_dx_allCuts->Fill(dx);
				}

			// Coincidence timing cut and vertex cut to resolve W well
				// if( fabs(diff - tdiff)<tdiff_max && fabs(bb_tr_vz[0])<=0.075 ){
				// 	h_W_cut->Fill( W );
				// } 

				// Preliminary HCal projections with single cut on W
				if( fabs(W - W_mean) < W2_sig_multiplier*W_sigma ){
				// if( ( W > (W_mean - 0.75*W_sigma) ) && ( W < (W_mean + 0.75*W_sigma) ) ){
				// if( W2 < 1.05 && W2 > 0.85){
					h_dx_wcut->Fill( dx );
					h_dy_wcut->Fill ( dy );
					h_dxdy_wcut->Fill( dy, dx );
					h_W_cut->Fill( W );
				}
				if( fabs(W2 - W2_mean) < W2_sig_multiplier*W2_sigma ){
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

					if( is_p ){
						h_p_proton->Fill(p_Nucleon);
					}
					if( is_n ){
						h_p_neutron->Fill(p_Nucleon);
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

							//count bad hcal clus e sorts for elastics:
							if( FindIndexOfMaxElement( hcal_clus_e, max_clus ) != 0 ){
								badHcalEsortElastics++;
							}
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
							
							//count bad hcal clus e sorts for elastics:
							if( FindIndexOfMaxElement( hcal_clus_e, max_clus ) != 0 ){
								badHcalEsortElastics++;
							}

						}
					}
				}
			//END OF FIDUCIAL Cut
				//Still should count elastic yields if we got this far.....
				if( !fiducial_cut ){
					elastic_yield++;

					//count bad hcal clus e sorts for elastics:
					if( FindIndexOfMaxElement( hcal_clus_e, max_clus ) != 0 ){
						badHcalEsortElastics++;
					}
				}
			}

			if( calc_W_only ){
				elastic_yield++;

				//count bad hcal clus e sorts for elastics:
				if( FindIndexOfMaxElement( hcal_clus_e, max_clus ) != 0 ){
					badHcalEsortElastics++;
				}
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

			h_W2_in = static_cast<TH1D*>(histo_infile->Get("h_W2"));
			h_W2copy = static_cast<TH1D*>(histo_infile->Get("h_W2"));

			h_dxdy_in = static_cast<TH2D*>(histo_infile->Get("h_dxdy"));
			h_p_Nucleon_in = static_cast<TH1D*>(histo_infile->Get("h_p_Nucleon"));
			p_Nucleon_mean = h_p_Nucleon_in->GetMean();

			h_p_proton_in = static_cast<TH1D*>(histo_infile->Get("h_p_proton"));
			p_proton_mean = h_p_proton_in->GetMean();

			h_p_neutron_in = static_cast<TH1D*>(histo_infile->Get("h_p_neutron"));
			p_neutron_mean = h_p_neutron_in->GetMean();

			//Reduce the range of dxdy to limit the excess plotting range and excess backround in areas far away from elastics
			h_dxdy_in->GetXaxis()->SetRangeUser(dy_min_y, dy_max_y);  //X-axis of Hcal corresponds to dx plot
			h_dxdy_in->GetYaxis()->SetRangeUser(dx_min_x, dx_max_x);	//Y-axis of HCal corresponds to dy plot

			h_dy_in = (TH1D*)(h_dxdy_in->ProjectionX())->Clone("h_dy_in");

			//Gotta turn the range back to normal
			h_dxdy_in->GetXaxis()->SetRangeUser(dxdy_min_y, dxdy_max_y);
			h_dxdy_in->GetYaxis()->SetRangeUser(dxdy_min_x, dxdy_max_x);


			h_W2_bg_from_sub = new TH1D("h_W2_bg_from_sub", "h_W2_bg_from_sub", W2_nbins, W2_min_x, W2_max_x);
			h_W2_sig_from_sub = new TH1D("h_W2_sig_from_sub", "h_W2_sig_from_sub", W2_nbins, W2_min_x, W2_max_x);
			h_dx_sig_from_sub = new TH1D("h_dx_sig_from_sub", "h_dx_sig_from_sub", dx_nbins, dx_min_x, dx_max_x);


	//----------------------------------------------------------
	//------------------ W2 Fits and Extractions ---------------
	//----------------------------------------------------------


			TCanvas *c_W2 = new TCanvas("c_W2", "c_W2", 600, 500);
			h_W2_in->Draw();

		// BG from rejecting middle "signal" points

			if( kine == 4 ){
				reject_min = 0.35;
				reject_max = 1.16;
			}
			if( kine == 8 ){
				reject_min = 0.30;
				reject_max = 1.07;
			}

			tf1_W2_bg_reject = new TF1("tf1_W2_bg_reject", fitPol4_with_reject, W2_min_x, W2_max_x, 5);
			tf1_W2_bg_reject->SetNpx(W2_nbins);
			tf1_W2_bg_reject->SetLineColor(9);
			h_W2_in->Fit("tf1_W2_bg_reject", "RMSE0+");
			tf1_W2_bg_reject->GetParameters(W2_bg_reject_par);

			tf1_W2_bg_pol4 = new TF1("tf1_W2_bg_pol4", "pol4", 0.0, W2_max_x);
			tf1_W2_bg_pol4->SetNpx(W2_nbins);
			tf1_W2_bg_pol4->SetLineColor(9);
			tf1_W2_bg_pol4->SetParameters(W2_bg_reject_par);

			h_tf1_W2_bg_pol4 = (TH1*)tf1_W2_bg_pol4->GetHistogram();
			h_tf1_W2_bg_pol4->Scale(1.0/h_tf1_W2_bg_pol4->Integral(), "width");
			h_tf1_W2_bg_pol4->Scale(1.0/h_tf1_W2_bg_pol4->Integral(), "height");
			h_tf1_W2_bg_pol4->Scale(tf1_W2_bg_pol4->GetMaximum()/h_tf1_W2_bg_pol4->GetMaximum() );

			h_W2_bg_pol4 = (TH1D*)h_tf1_W2_bg_pol4->Clone("h_W2_bg_pol4");
			h_W2_bg_pol4->SetLineColor(9);
			h_W2_bg_pol4->Draw("hist+same");

		// FULL FIT
			tf1_W2_fullFit = new TF1("tf1_W2_fullFit", fullFitFunction, W2_min_x, W2_max_x, polN+4);
			tf1_W2_fullFit->SetNpx(W2_nbins);

			//Set tf1_W2_fullFit parameter names:
			tf1_W2_fullFit->SetParName(0, "W2_fullFit Gaus Norm");
			tf1_W2_fullFit->SetParName(1, "W2_fullFit Gaus Mean");
			tf1_W2_fullFit->SetParName(2, "W2_fullFit Gaus Sigma");
			for(int param = 3; param < (polN + 4); param++ ){
				tf1_W2_fullFit->SetParName(param, Form("W2_fullFit BG pol%i p%i", polN, param-3) );
			}
			// double W2_gaus_norm_fit_min = 0.9*h_W2_in->GetMaximum();
			double W2_gaus_norm_fit_min = 0.0;
			//Set parameter limits for signal gaussian:

			if( kine == 4 ){
				tf1_W2_fullFit->SetParLimits(0, W2_gaus_norm_fit_min, h_W2_in->GetMaximum());
				tf1_W2_fullFit->SetParLimits(1, 0.775, 0.9);
				tf1_W2_fullFit->SetParLimits(2, 0.0, 0.2);			
			}
			if( kine == 8 ){
				tf1_W2_fullFit->SetParLimits(0, W2_gaus_norm_fit_min, h_W2_in->GetMaximum());
				tf1_W2_fullFit->SetParLimits(1, 0.7, 0.775);
				tf1_W2_fullFit->SetParLimits(2, 0.0, 0.2);			
			}

			//Pull parameters from the BG reject fit
			tf1_W2_fullFit->FixParameter(3, W2_bg_reject_par[0]);
			tf1_W2_fullFit->FixParameter(4, W2_bg_reject_par[1]); //168
			tf1_W2_fullFit->FixParameter(5, W2_bg_reject_par[2]);
			tf1_W2_fullFit->FixParameter(6, W2_bg_reject_par[3]);
			tf1_W2_fullFit->FixParameter(7, W2_bg_reject_par[4]); //1932

			// tf1_W2_fullFit->FixParameter(3, 50);
			// tf1_W2_fullFit->FixParameter(4, 188); //168
			// tf1_W2_fullFit->FixParameter(5, 1405);
			// tf1_W2_fullFit->FixParameter(6, -2865);
			// tf1_W2_fullFit->FixParameter(7, 1915); //1932

			// tf1_W2_fullFit->FixParameter(3, 52.47);
			// tf1_W2_fullFit->FixParameter(4, -93.50); //168
			// tf1_W2_fullFit->FixParameter(5, 970);
			// tf1_W2_fullFit->FixParameter(6, -1352.5);
			// tf1_W2_fullFit->FixParameter(7, 905); //1932

			// tf1_W2_fullFit->SetParLimits(3, 47.5, 48.5);// 48, 52.5
			// tf1_W2_fullFit->SetParLimits(4, -92.5, -91.5); //-94, -92
			// tf1_W2_fullFit->SetParLimits(5, 989, 991); // 988, 992
			// tf1_W2_fullFit->SetParLimits(6, -1359, -1357);//-1362, -1359
			// tf1_W2_fullFit->SetParLimits(7, 894, 896); //923, 925

			h_W2_in->Fit("tf1_W2_fullFit", "RMSE+");
			W2_sig_from_sub_integral_error = tf1_W2_fullFit->IntegralError(W2_min_x, W2_max_x)/tf1_W2_fullFit->GetHistogram()->GetBinWidth(1);

			tf1_W2_fullFit->SetLineColor(3);
			tf1_W2_fullFit->SetLineStyle(3);

			tf1_W2_fullFit->GetParameters(par);
			tf1_W2_fullFit->GetParameters(W2_fullFit_par);
			W2_fullFit_par_errors[0] = tf1_W2_fullFit->GetParError(0);
			W2_fullFit_par_errors[1] = tf1_W2_fullFit->GetParError(1);
			W2_fullFit_par_errors[2] = tf1_W2_fullFit->GetParError(2);

			W2_gaus_min_x = W2_fullFit_par[1] - W2_sig_multiplier*fabs(W2_fullFit_par[2]);
			W2_gaus_max_x = W2_fullFit_par[1] + W2_sig_multiplier*fabs(W2_fullFit_par[2]);

			W2_gaus_min_bin = h_W2_in->GetXaxis()->FindBin(W2_gaus_min_x);
			W2_gaus_max_bin = h_W2_in->GetXaxis()->FindBin(W2_gaus_max_x);

		//------ W2 Sig from Sub
			for(int bin = 0; bin < W2_nbins; bin++ ){
				double sub_val = h_W2_in->GetBinContent(bin) - h_W2_bg_pol4->GetBinContent(bin);

				if( (bin >= W2_gaus_min_bin) && ( bin <= W2_gaus_max_bin ) ){
					h_W2_sig_from_sub->SetBinContent(bin, sub_val);			
				}

			}
			h_W2_sig_from_sub->SetLineColor(6);
			h_W2_sig_from_sub->Draw("same");

			tf1_W2_gaus = new TF1("tf1_W2_gaus", "gaus", 0.0, W2_max_x);
			tf1_W2_gaus->SetNpx(W2_nbins);
			tf1_W2_gaus->SetLineColor(6);
			tf1_W2_gaus->SetLineStyle(6);
			tf1_W2_gaus->FixParameter(0, par[0]);
			tf1_W2_gaus->FixParameter(1, par[1]);
			tf1_W2_gaus->FixParameter(2, par[2]);

			W2_gaus_integral = tf1_W2_gaus->Integral(0.0, W2_max_x)/tf1_W2_gaus->GetHistogram()->GetBinWidth(1);
			W2_gaus_integral_error = tf1_W2_gaus->IntegralError(W2_gaus_min_x, W2_gaus_max_x)/tf1_W2_gaus->GetHistogram()->GetBinWidth(1);

			tf1_W2_gaus_err = new TF1("tf1_W2_gaus_err", "gaus", W2_gaus_min_x, W2_gaus_max_x);
			tf1_W2_gaus_err->SetNpx(W2_nbins);
			tf1_W2_gaus_err->SetParLimits(0, 0.999*par[0], 1.001*par[0]);
			tf1_W2_gaus_err->SetParLimits(1, 0.999*par[1], 1.001*par[1]);
			tf1_W2_gaus_err->SetParLimits(2, 0.999*par[2], 1.001*par[2]);
			h_W2_sig_from_sub->Fit("tf1_W2_gaus_err", "RMSE0");

			W2_sig_from_sub_par0 = tf1_W2_gaus_err->GetParameter(0);
			W2_sig_from_sub_par1 = tf1_W2_gaus_err->GetParameter(1);
			W2_sig_from_sub_par2 = tf1_W2_gaus_err->GetParameter(2);

			W2_sig_from_sub_par0err = tf1_W2_gaus_err->GetParError(0);
			W2_sig_from_sub_par1err = tf1_W2_gaus_err->GetParError(1);
			W2_sig_from_sub_par2err = tf1_W2_gaus_err->GetParError(2);
			W2_sig_from_sub_integral_error = tf1_W2_gaus_err->IntegralError(W2_gaus_min_x, W2_gaus_max_x)/tf1_W2_gaus_err->GetHistogram()->GetBinWidth(1);

			W2_sig_from_sub_integral = h_W2_sig_from_sub->Integral();

			if( W2_gaus_integral_error == 0.0 ){
				W2_gaus_integral_error = calc_gaus_error(W2_sig_from_sub_par0, W2_sig_from_sub_par0err, W2_sig_from_sub_par2, W2_sig_from_sub_par2err, h_W2_sig_from_sub->GetBinWidth(1) );
				W2_sig_from_sub_integral_error = W2_gaus_integral_error;
			}

			TLegend *tl_W2_bg_sub = new TLegend(0.15, 0.7, 0.45, 0.85);
			tl_W2_bg_sub->AddEntry(tf1_W2_fullFit, "Total W2 fit");
			tl_W2_bg_sub->AddEntry(h_W2_sig_from_sub, "W2 peak from BG sub.");
			tl_W2_bg_sub->AddEntry(tf1_W2_bg_pol4, "BG fit");
			tl_W2_bg_sub->Draw("Same");

			TPaveText *tp_W2_bg_sub = new TPaveText(0.15, 0.60, 0.45, 0.69, "NDCbr");
			tp_W2_bg_sub->AddText(Form("W2 count = %i +/- %i", int(W2_sig_from_sub_integral), int(W2_sig_from_sub_integral_error) ));
			tp_W2_bg_sub->AddText(Form("Total fit #chi^{2} = %.1f", tf1_W2_fullFit->GetChisquare() ));
			tp_W2_bg_sub->Draw("same");

			cout << endl << "-------------------------------" << endl;
			cout << "W2 gaus_min: " << W2_gaus_min_x << endl;
			cout << "W2 gaus_max: " << W2_gaus_max_x << endl << endl;
			cout << "---" << endl;
			cout << "W2 fullFit gaus integral: " << W2_gaus_integral << " +/- " << W2_gaus_integral_error <<  endl;
			cout << "---" << endl;
			cout << "W2 Signal from BG sub: " << W2_sig_from_sub_integral << " +/- " << W2_sig_from_sub_integral_error << endl;
			cout << "W2 signal from sub error: " << W2_sig_from_sub_integral_error << endl;
			cout << "-------------------------------" << endl << endl;


	//----------------------------------------------------------
	//------------------ dy Fits and Extractions ---------------
	//----------------------------------------------------------
			dypolN = 2;

			TCanvas *c_dy = new TCanvas("c_dy", "c_dy", 600, 500);
			h_dy_in->GetXaxis()->SetRangeUser(-0.4, 0.4);
			h_dy_in->Draw();

			tf1_dy_fullFit = new TF1("tf1_dy_fullFit", fullDYFitFunction, dy_min_y, dy_max_y, 7);
			tf1_dy_fullFit->SetNpx(dy_nbins);

			tf1_dy_fullFit->SetParName(0, "dy Gaus Norm");
			tf1_dy_fullFit->SetParName(1, "dy Gaus Mean");
			tf1_dy_fullFit->SetParName(2, "dy Gaus Sigma");
			tf1_dy_fullFit->SetParName(3, "dy BG pol3 p0");
			tf1_dy_fullFit->SetParName(4, "dy BG pol3 p1");
			tf1_dy_fullFit->SetParName(5, "dy BG pol3 p2");
			if(dypolN == 3){
				tf1_dy_fullFit->SetParName(6, "dy BG pol3 p3");		
			}

			tf1_dy_fullFit->SetParLimits(0, 0.5*h_dy_in->GetMaximum(), 0.93*h_dy_in->GetMaximum());
			tf1_dy_fullFit->SetParLimits(1, -0.05, 0.05);
			tf1_dy_fullFit->SetParLimits(2, 0.0, 0.20);

			h_dy_in->Fit("tf1_dy_fullFit", "RMSE+");
			tf1_dy_fullFit->GetParameters(dy_fullFit_par);
			dy_fullFit_integral_error = tf1_dy_fullFit->IntegralError(dy_min_y, dy_max_y)/tf1_dy_fullFit->GetHistogram()->GetBinWidth(1);

			if( dypolN == 2){
				tf1_dy_pol = new TF1("tf1_dy_pol", fitPol2, dy_min_y, dy_max_y, 3);			
			}
			if( dypolN == 3){
				tf1_dy_pol = new TF1("tf1_dy_pol", fitPol3, dy_min_y, dy_max_y, 4);			
			}
			tf1_dy_pol->SetNpx(dy_nbins);
			tf1_dy_pol->SetLineColor(9);

			for( int param = 0; param < (dypolN+1); param++ ){
				tf1_dy_pol->FixParameter(param, dy_fullFit_par[param+3]);
			}

			tf1_dy_pol->Draw("hist+same");

			tf1_dy_gaus = new TF1("tf1_dy_gaus", "gaus", dy_min_y, dy_max_y);
			tf1_dy_gaus->SetNpx(2.0*dy_nbins);
			h_dy_in->Fit("tf1_dy_gaus", "RMSE0+");
			dy_gaus_integral_error = tf1_dy_gaus->Integral(dy_min_y, dy_max_y)/tf1_dy_gaus->GetHistogram()->GetBinWidth(1);
			tf1_dy_gaus->SetLineColor(6);

			tf1_dy_gaus->FixParameter(0, dy_fullFit_par[0]);
			tf1_dy_gaus->FixParameter(1, dy_fullFit_par[1]);
			tf1_dy_gaus->FixParameter(2, dy_fullFit_par[2]);
			tf1_dy_gaus->Draw("hist+same");

			dy_fit_min = dy_fullFit_par[1] - dy_sig_multiplier*fabs(dy_fullFit_par[2]);
			dy_fit_max = dy_fullFit_par[1] + dy_sig_multiplier*fabs(dy_fullFit_par[2]);

			TLine *tl_dy_min = new TLine(dy_fit_min, 0, dy_fit_min, 1.05*h_dy_in->GetMaximum() );
			tl_dy_min->SetLineStyle(6);
			tl_dy_min->Draw("same");

			TLine *tl_dy_max = new TLine(dy_fit_max, 0, dy_fit_max, 1.05*h_dy_in->GetMaximum() );
			tl_dy_max->SetLineStyle(6);
			tl_dy_max->Draw("same");

			dy_integral = tf1_dy_gaus->Integral(dy_min_y, dy_max_y)/tf1_dy_gaus->GetHistogram()->GetBinWidth(1);

			TLegend *tl_dy;
			if( kine == 4 ){
				tl_dy = new TLegend(0.6, 0.6, 0.9, 0.75);
			}
			if( kine == 8 ){
				tl_dy = new TLegend(0.15, 0.65, 0.45, 0.80);
			}

			tl_dy->AddEntry(tf1_dy_fullFit, "Totaly dy fit");
			tl_dy->AddEntry(tf1_dy_gaus, "dy fit");
			tl_dy->AddEntry(tf1_dy_pol, "BG fit");
			tl_dy->Draw("same");

			TPaveText *tp_dy;
			if( kine == 4 ){
				tp_dy = new TPaveText(0.6, 0.5, 0.9, 0.59, "NDCbr");
			}
			if( kine == 8 ){
				tp_dy = new TPaveText(0.15, 0.65, 0.45, 0.80, "NDCbr");				
			}
			tp_dy->AddText(Form("Total fit #chi^{2} = %.1f", tf1_dy_fullFit->GetChisquare()));
			tp_dy->Draw("same");


			cout << endl << "-----------------------------------" << endl;
			cout << "dy_fullFit chiSquare: " << tf1_dy_fullFit->GetChisquare() << endl;
			cout << "dy_mean: " << dy_fullFit_par[1] << endl;
			cout << "dy_sigma: " << dy_fullFit_par[2] << endl;
			cout << "dy_sig_multiplier: " << dy_sig_multiplier << endl;
			cout << "dy_min: " << dy_fit_min << ", dy_max: " << dy_fit_max << endl;
			cout << "---" << endl;
			cout << "dy_integral: " << dy_integral << endl;
			cout << "---" << endl;

		//----------------------------------------------------------
		//------------------ dx Fits and Extractions ---------------
		//----------------------------------------------------------
			h_dxdy_in->GetXaxis()->SetRangeUser(dy_fit_min, dy_fit_max);
			h_dxdy_in->GetYaxis()->SetRangeUser(dx_min_x, dx_max_x);
			h_dx_in = (TH1D*)(h_dxdy_in->ProjectionY())->Clone("h_dx_in");
			
			dxpolN = 2;

			TCanvas *c_dx = new TCanvas("c_dx", "c_dx", 600, 500);
			h_dx_in->SetMinimum(0);
			h_dx_in->Draw();

			tf1_dx_fullFit = new TF1("tf1_dx_fullFit", fullDXFitFunction, dx_min_x, dx_max_x, dxpolN + 4);
			tf1_dx_fullFit->SetNpx(dx_nbins);

			tf1_dx_fullFit->SetParName(0, "dx Gaus Norm");
			tf1_dx_fullFit->SetParName(1, "dx Gaus Mean");
			tf1_dx_fullFit->SetParName(2, "dx Gaus Sigma");
			tf1_dx_fullFit->SetParName(3, Form("dx BG pol%i p0", dxpolN));
			tf1_dx_fullFit->SetParName(4, Form("dx BG pol%i p1", dxpolN));
			tf1_dx_fullFit->SetParName(5, Form("dx BG pol%i p2", dxpolN));

			if( dxpolN == 2 ){
				tf1_dx_fullFit->SetParLimits(5, -10000.0, 0.0);
			}

			if( dxpolN == 3){
				tf1_dx_fullFit->SetParName(6, Form("dx BG pol%i p3", dxpolN));
			}

			tf1_dx_fullFit->SetParLimits(0, 0.5*h_dx_in->GetMaximum(), 0.95*h_dx_in->GetMaximum());
			tf1_dx_fullFit->SetParLimits(1, -0.5, 0.5);

			if( kine == 4 ){
				tf1_dx_fullFit->SetParLimits(2, 0.0, 0.15);				
			}
			if( kine == 8 ){
				tf1_dx_fullFit->SetParLimits(2, 0.0, 0.085);
			}
			
			h_dx_in->Fit("tf1_dx_fullFit", "RMSE+");
			tf1_dx_fullFit->GetParameters(dx_fullFit_par);
			dx_fullFit_par_errors[0] = tf1_dx_fullFit->GetParError(0);
			dx_fullFit_par_errors[1] = tf1_dx_fullFit->GetParError(1);
			dx_fullFit_par_errors[2] = tf1_dx_fullFit->GetParError(2);

			dx_fullFit_integral_error = tf1_dx_fullFit->IntegralError(dx_min_x, dx_max_x)/tf1_dx_fullFit->GetHistogram()->GetBinWidth(1);

			tf1_dx_gaus = new TF1("tf1_dx_gaus", "gaus", dx_min_x, dx_max_x);
			tf1_dx_gaus->SetNpx(dx_nbins);
			h_dx_in->Fit("tf1_dx_gaus", "RMSE0+");
			dx_gaus_integral_error = tf1_dx_gaus->IntegralError(dx_min_x, dx_max_x)/tf1_dx_gaus->GetHistogram()->GetBinWidth(1);

			tf1_dx_gaus->FixParameter(0, dx_fullFit_par[0]);
			tf1_dx_gaus->FixParameter(1, dx_fullFit_par[1]);
			tf1_dx_gaus->FixParameter(2, dx_fullFit_par[2]);

			dx_gaus_norm = dx_fullFit_par[0];
			dx_gaus_mean = dx_fullFit_par[1];
			dx_gaus_sigma = dx_fullFit_par[2];

			if( dxpolN == 2 ){
				tf1_dx_pol = new TF1("tf1_dx_pol", fitPol2, dx_min_x, dx_max_x, 3);		
			}
			if( dxpolN == 3 ){
				tf1_dx_pol = new TF1("tf1_dx_pol", fitPol3, dx_min_x, dx_max_x, 4);		
			}

			tf1_dx_pol->SetNpx(dx_nbins);
			tf1_dx_pol->SetLineColor(9);

			for(int param = 0; param < (dxpolN + 1); param++ ){
				tf1_dx_pol->FixParameter(param, dx_fullFit_par[param+3]);
			}
			tf1_dx_pol->Draw("hist+same");

			h_tf1_dx_pol = (TH1*)tf1_dx_pol->GetHistogram();
			h_tf1_dx_pol->Scale(1.0/h_tf1_dx_pol->Integral(), "width");
			h_tf1_dx_pol->Scale(1.0/h_tf1_dx_pol->Integral(), "height");
			h_tf1_dx_pol->Scale(tf1_dx_pol->GetMaximum()/h_tf1_dx_pol->GetMaximum());

			h_dx_pol = (TH1D*)h_tf1_dx_pol->Clone("h_dx_pol");


		//Now lets get the "signal" Gaussian by subracting the polN BG histogram from h_dx_in
		//We base this off of the gaussian determined from the full fit 
			dx_gaus_min_x = dx_gaus_mean - dx_sig_multiplier*fabs(dx_gaus_sigma);
			dx_gaus_max_x = dx_gaus_mean + dx_sig_multiplier*fabs(dx_gaus_sigma);

			dx_gaus_min_bin = h_dx_in->GetXaxis()->FindBin(dx_gaus_min_x);
			dx_gaus_max_bin = h_dx_in->GetXaxis()->FindBin(dx_gaus_max_x);

			cout << "dx_gaus_mean: " << dx_gaus_mean << endl;
			cout << "dx_gaus_sigma: " << dx_gaus_sigma << endl;

			cout << endl << "dx_gaus --- min x: " << dx_gaus_min_x << ", max x: " << dx_gaus_max_x << endl;
			cout << "dx_gaus --- min bin: " << dx_gaus_min_bin << ", max bin: " << dx_gaus_max_bin << endl << endl;


			for( int bin = 0; bin <= dx_nbins; bin++ ){
				double sub_val = h_dx_in->GetBinContent(bin) - h_dx_pol->GetBinContent(bin);

				if( bin >= dx_gaus_min_bin && bin <= dx_gaus_max_bin ){
					h_dx_sig_from_sub->SetBinContent(bin, sub_val);
				}
			}

			h_dx_sig_from_sub->SetLineColor(6);
			h_dx_sig_from_sub->Draw("same");

			tf1_dx_bgSub_gaus = new TF1("tf1_dx_bgSub_gaus", "gaus", dx_min_x, dx_max_x);
			tf1_dx_bgSub_gaus->SetNpx(dx_nbins);
			tf1_dx_bgSub_gaus->FixParameter(0, dx_fullFit_par[0]);
			tf1_dx_bgSub_gaus->FixParameter(1, dx_fullFit_par[1]);
			tf1_dx_bgSub_gaus->FixParameter(2, dx_fullFit_par[2]);

			h_dx_sig_from_sub->Fit("tf1_dx_bgSub_gaus", "RMSE0+");

			dx_bgSub_gaus_integral_error = tf1_dx_bgSub_gaus->IntegralError(dx_min_x, dx_max_x)/tf1_dx_bgSub_gaus->GetHistogram()->GetBinWidth(1);

			cout << "dx_gaus_integral: " << h_dx_sig_from_sub->Integral() << endl;
			cout << "dx_gaus_integral error: " << dx_bgSub_gaus_integral_error << endl;
			// dx_addsub_error = sqrt( pow( dx_fullFit_integral_error , 2) + pow(dx_gaus_integral_error, 2) );
			dx_addsub_error = sqrt( pow( dx_fullFit_integral_error , 2) + pow(dx_bgSub_gaus_integral_error, 2) );

			TLegend *tl_dx = new TLegend(0.15, 0.7, 0.45, 0.85);
			tl_dx->AddEntry(tf1_dx_fullFit, "Total dx fit");
			tl_dx->AddEntry(tf1_dx_bgSub_gaus, "dx peak from BG sub");
			tl_dx->AddEntry(tf1_dx_pol, "BG fit");
			tl_dx->Draw("same");

			TPaveText *tp_dx = new TPaveText(0.15, 0.55, 0.45, 0.69, "NDCbr");
			tp_dx->AddText(Form("dx count = %i +/- %i", int(h_dx_sig_from_sub->Integral()), int(dx_bgSub_gaus_integral_error) ));
			tp_dx->AddText(Form("Total fit #chi^{2} = %.1f", tf1_dx_fullFit->GetChisquare()));
			tp_dx->Draw("same");

		/// -----------------------------------------------------
		/////	---------------------print stuff

			dx_fullFit_integral = tf1_dx_fullFit->Integral(dx_min_x, dx_max_x)/tf1_dx_fullFit->GetHistogram()->GetBinWidth(1);
			dx_pol_integral = tf1_dx_pol->Integral(dx_min_x, dx_max_x)/tf1_dx_pol->GetHistogram()->GetBinWidth(1);
			dx_gaus_integral = h_dx_sig_from_sub->Integral();

			det_eff_FINAL = 100.0*dx_gaus_integral/W2_sig_from_sub_integral;

			//Only including error from ratio for efficiency  W2_gaus_integral_error
			// det_eff_total_error = det_eff_FINAL*sqrt( pow( W2_gaus_integral_error/W2_gaus_integral, 2) + pow( dx_bgSub_gaus_integral_error/dx_gaus_integral, 2) );

			//OLD
			// det_eff_total_error = det_eff_FINAL*sqrt( 100.0*(pow( (dx_gaus_integral_error/dx_gaus_integral), 2) + pow( (W2_sig_from_sub_integral_error/W2_sig_from_sub_integral) , 2) ));
			//NEW
			det_eff_total_error = det_eff_FINAL*sqrt( 100.0*(pow( (dx_bgSub_gaus_integral_error/dx_gaus_integral), 2) + pow( (W2_sig_from_sub_integral_error/W2_sig_from_sub_integral) , 2) ));

			//calculated with sqrt(N)
			double sqrt_N_err = det_eff_FINAL*sqrt( 100.0*(pow( sqrt(dx_gaus_integral)/dx_gaus_integral, 2 ) + pow( sqrt(W2_gaus_integral)/W2_gaus_integral, 2) ));

			cout << endl << "-------------------------------" << endl;
			cout << "Nucleon Momentum = " << p_Nucleon_mean << endl;
			cout << "-------------------------------" << endl;
			cout << "dx fullFit ChiSquare: " << tf1_dx_fullFit->GetChisquare() << endl;
			cout << "Full Integral: " << dx_fullFit_integral << endl;
			cout << "BG integral: " << dx_pol_integral << endl;
			cout << "---" << endl;
			cout << "dx gaus fit integral: " << tf1_dx_gaus->Integral(dx_min_x, dx_max_x)/tf1_dx_gaus->GetHistogram()->GetBinWidth(1) << endl;
			// cout << "dx gaus integral from function: " << gaus_integral(dx_gaus_norm, dx_gaus_sigma, h_dx_in->GetBinWidth(1)) << endl << endl;
			cout << "-------------------------------" << endl;
			// cout << "dx signal integral: " << dx_gaus_integral << " +/- " << dx_gaus_integral_error << endl;
			cout << "dx signal integral: " << dx_gaus_integral << " +/- " << dx_bgSub_gaus_integral_error << endl;
			cout << "---" << endl;
			cout << "W2 signal integral: " << W2_sig_from_sub_integral << " +/- " << W2_sig_from_sub_integral_error << endl;
			cout << "-------------------------------" << endl;
			// cout << "HCal Detector Efficiency = " << det_eff_FINAL << "% +/- " << det_eff_total_error << "% " << endl << endl;
	
			defficiency = 100.0*sqrt( fabs(det_eff_FINAL*( 100.0 - det_eff_FINAL )/(1.0*dx_gaus_integral)));
			cout << "defficiency: " << defficiency << endl;

			//Calculated with sqrt(N)
			cout << "HCal Detector Efficiency = " << det_eff_FINAL << "% +/- " << det_eff_total_error << "% " << endl << endl;
			cout << "---" << endl;


		//----------------------------------------------------------
		//------------------ dxdy Fits and Extractions ---------------
		//----------------------------------------------------------
			if( true ){

				h_dxdy_in->GetXaxis()->SetRangeUser(dxdy_min_x, dxdy_max_x);
				h_dxdy_in->GetYaxis()->SetRangeUser(dxdy_min_y, dxdy_max_y);

				TCanvas *c_dxdy = new TCanvas("c_dxdy", "c_dxdy", 600, 500);
				h_dxdy_in->Draw("colz");
				TLine *tl_dxdy_dy_min = new TLine(dy_fit_min, dxdy_min_y, dy_fit_min, dxdy_max_y);
				tl_dxdy_dy_min->SetLineColor(41);
				tl_dxdy_dy_min->SetLineStyle(6);
				tl_dxdy_dy_min->SetLineWidth(3);
				tl_dxdy_dy_min->Draw("same");

				TLine *tl_dxdy_dy_max = new TLine(dy_fit_max, dxdy_min_y, dy_fit_max, dxdy_max_y);
				tl_dxdy_dy_max->SetLineColor(41);
				tl_dxdy_dy_max->SetLineStyle(6);
				tl_dxdy_dy_max->SetLineWidth(3);
				tl_dxdy_dy_max->Draw("same");

				// TLine *tl_dxdy_dx_min = new TLine(dxdy_min_x, dx_min_x, dxdy_max_x, dx_min_x);
				// tl_dxdy_dx_min->SetLineColor(29);
				// tl_dxdy_dx_min->SetLineStyle(5);
				// tl_dxdy_dx_min->SetLineWidth(3);
				// tl_dxdy_dx_min->Draw("same");

				// TLine *tl_dxdy_dx_max = new TLine(dxdy_min_x, dx_max_x, dxdy_max_x, dx_max_x);
				// tl_dxdy_dx_max->SetLineColor(29);
				// tl_dxdy_dx_max->SetLineStyle(5);
				// tl_dxdy_dx_max->SetLineWidth(3);
				// tl_dxdy_dx_max->Draw("same");		
			}


		}
	}

		auto total_time_end = chrono::high_resolution_clock::now();
		auto total_time_duration = chrono::duration_cast<chrono::minutes>(total_time_end - total_time_start);
		cout << endl << "---------------------------------------------------" << endl;
		cout << "Total time for analysis: " << total_time_duration.count() << " minutes. " << endl;
		cout << endl << "Outfile: " << outfile_name.Data() << endl;
		cout << endl << "Infile: " << infile_name.Data() << endl;
}