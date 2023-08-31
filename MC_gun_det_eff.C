//Script to calculate the detection efficiency for neutrons and protons for MC output
//Different
//First need to calibrate so that we can pull mean values for HCalE and p_Nucleon
//Written by J. Boyd - May 22, 2023

TH1D* MakeHisto(int, int, double, double, std::string);

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
#include "/w/halla-scshelf2102/sbs/jboyd/include/calc_functions.h"

double par_n[5], par_p[5];
double par_n_error[5], par_p_error[5];
double det_eff_err_p = 0.0;
double det_eff_err_n = 0.0;

double calc_pol_error( double *par_err, int polN){

	vector<double> vec_errors(polN, 0.0);
	for(int i = 0; i < polN; i++ ){
		vec_errors[i] = par_err[i];
	}

	double error;
	double temp_error = 0.0;

	for( int elem = 0; elem < vec_errors.size(); elem++ ){
		temp_error += pow( vec_errors[elem], 2 );
	}

	error = sqrt( temp_error );
	return error;
}

Double_t fit_det_eff( Double_t *x, Double_t *par){
	double det_eff_fit = 0.0;

	det_eff_fit = par[0] + (par[1]/(pow( par[2], x[0])));

	return det_eff_fit;
}

Double_t calc_det_eff(double p_central, double *par ){

	double det_eff = 0.0;

	det_eff = par[0] + par[1]*p_central + par[2]*pow( p_central, 2 ) + par[3]*pow( p_central, 3 ) + par[4]*pow( p_central, 4 );

	return det_eff;

}

int kine = 4;
int sbsfieldscale = 0;

int nBins = 80;
double p_min = 1.0;
double p_max = 9.0;

const int pN_bins = 80;
const double pN_min = 1.0;
const double pN_max = 2.2;

const int HCal_bins = 220;
const double HCal_min = 0.0;
const double HCal_max = 2.2;

const double det_eff_yaxis_min = 80.0;
const double n_xaxis_min = 1.4;
const double p_xaxis_min = 1.3;

double p_central_sbs4 = 2.35;
double p_central_sbs7 = 6.20;
double p_central_sbs8 = 3.22;
double p_central_sbs9 = 3.21;
double p_central_sbs14 = 4.84;
double p_central_kine = 0.0;

double p_sbs4_neut_min, p_sbs4_neut_max;
double p_sbs7_neut_min, p_sbs7_neut_max;
double p_sbs8_neut_min, p_sbs8_neut_max;
double p_sbs9_neut_min, p_sbs9_neut_max;
double p_sbs14_neut_min, p_sbs14_neut_max;

double p_sbs4_prot_min, p_sbs4_prot_max;
double p_sbs7_prot_min, p_sbs7_prot_max;
double p_sbs8_prot_min, p_sbs8_prot_max;
double p_sbs9_prot_min, p_sbs9_prot_max;
double p_sbs14_prot_min, p_sbs14_prot_max;

vector<double> p_sbs4_neut, p_sbs4_prot, p_sbs7_neut, p_sbs7_prot, p_sbs8_neut, p_sbs8_prot;
vector<double> p_sbs9_neut, p_sbs9_prot, p_sbs14_neut, p_sbs14_prot;

double histo_max = HCal_max;

//Select Nucleon to calculate from
//"neutron", "proton", "both"
TString Nucleon_Select = "both";
TString nucl_type = "";

Double_t HCal_theta_rad = lookup_HCal_angle_by_kine(kine, "rad");
Double_t HCal_dist = lookup_HCal_dist_by_kine(kine);

Double_t sbs_theta_rad = lookup_SBS_angle_by_kine(kine, "rad");

//Filenames and related
TString outfile_name, rootfile_dir, proton_infile_name, neutron_infile_name;
vector<TString> proton_infilenames_vec, neutron_infilenames_vec;
TFile *outfile, *proton_infile, *neutron_infile;
TChain *TC;

//Histograms
TH1D *h_dx_p, *h_dy_p, *h_dx_n, *h_dy_n;
TH1D *h_HCalE_p, *h_HCalE_n, *h_det_eff_p, *h_det_eff_n, *h_det_eff_pn;
TH1D *h_n_det_eff, *h_p_det_eff;

TH2D *h2_xy_p, *h2_xy_exp_p, *h2_HCal_pN_p, *h2_HCal_pN_cut_p;
TH2D *h2_xy_n, *h2_xy_exp_n, *h2_HCal_pN_n, *h2_HCal_pN_cut_n;

TF1 *tf_det_eff_n, *tf_det_eff_p;

TAxis *tax_h_det_eff_n, *tax_h_det_eff_p;
TPaveText *tpt_det_eff_fits, *tpt_det_eff_bins;

TProfile *tprof_h2_HCalE_pN_n, *tprof_h2_HCalE_pN_p;

auto gr_det_eff_n = new TGraph();
auto gr_det_eff_p = new TGraph();

//FITS
TF1 *tf_gaus_n, *tf_gaus_p;

double n_efficiency = 0.0;
double p_efficiency = 0.0;

//Efficiency calculations:
double prot_det_eff_sbs4 = 0.0, prot_det_eff_sbs7 = 0.0, prot_det_eff_sbs8 = 0.0, prot_det_eff_sbs9 = 0.0, prot_det_eff_sbs14 = 0.0;
double neut_det_eff_sbs4 = 0.0, neut_det_eff_sbs7 = 0.0, neut_det_eff_sbs8 = 0.0, neut_det_eff_sbs9 = 0.0, neut_det_eff_sbs14 = 0.0;

double prot_bin_det_eff_sbs4 = 0.0, prot_bin_det_eff_sbs7 = 0.0, prot_bin_det_eff_sbs8 = 0.0, prot_bin_det_eff_sbs9 = 0.0, prot_bin_det_eff_sbs14 = 0.0;
double neut_bin_det_eff_sbs4 = 0.0, neut_bin_det_eff_sbs7 = 0.0, neut_bin_det_eff_sbs8 = 0.0, neut_bin_det_eff_sbs9 = 0.0, neut_bin_det_eff_sbs14 = 0.0;

//TBranch variables
//HCal
Double_t HCal_id, HCal_e, HCal_x, HCal_y, HCal_row, HCal_col, HCal_tdc, HCal_atime;
Double_t mc_omega, mc_sigma, mc_fnucl;

//MC
Double_t mc_p, mc_px, mc_py, mc_pz, mc_vx, mc_vy, mc_vz, mc_nucl, mc_posx, mc_posy;

long nEvents = 0;
int max_iter = 1000;

bool hit_on_HCal = false;
int HCal_miss_cnt_p = 0;
int HCal_miss_cnt_n = 0;

int HCal_hit_cnt_p = 0;
int HCal_hit_cnt_n = 0;

TH1D *h_HCalE_array_n[pN_bins];
TH1D *h_HCalE_array_cut_n[pN_bins];

TH1D *h_HCalE_array_p[pN_bins];
TH1D *h_HCalE_array_cut_p[pN_bins];

vector<double> threshold_n, threshold_p;

void MC_gun_det_eff(){

	auto total_time_start = high_resolution_clock::now();
	TStopwatch *StopWatch = new TStopwatch();

	gStyle->SetPalette(55);
	cout << "--------------------------------------" << endl;
	cout << "Analysis started for MC-based Detection Efficiencies. " << endl;
	cout << "--------------------------------------" << endl;

	cout << "Run parameters: " << endl;
	cout << "Kinematic: SBS" << kine << endl;
	cout << "SBS Field Scale: " << sbsfieldscale << endl;

	switch( kine ){
		case 4: p_central_kine = 2.35; break;
		case 7: p_central_kine = 6.20; break;
		case 8: p_central_kine = 3.22; break;
		case 9: p_central_kine = 3.21; break;
		case 14: p_central_kine = 4.84; break;
		default: 
			p_central_kine = 0.0;
	}

	rootfile_dir = "/lustre19/expphy/volatile/halla/sbs/jboyd/simulation/out_dir/MC_REPLAY_OUT_DIR";

	if( Nucleon_Select == "proton" || Nucleon_Select == "both"){
		list_files(rootfile_dir, proton_infilenames_vec, "replayed_jb_gmn_SBS4_LD2_mag0_100uA_pGun_Px_300k_v2_job");
		// proton_infile_name = Form("%s/replayed_jb_gmn_SBS4_LD2_mag0_100uA_pGun_Px_300k_v2_job0.root", rootfile_dir.Data());
		// proton_infile = new TFile(proton_infile_name.Data(), "READ");	
	}

	if( Nucleon_Select == "neutron" || Nucleon_Select == "both"){
		list_files(rootfile_dir, neutron_infilenames_vec, "replayed_jb_gmn_SBS4_LD2_mag0_100uA_nGun_Px_300k_v2_job");
	// 	neutron_infile_name = Form("%s/replayed_jb_gmn_SBS%i_LD2_mag%i_100uA_nGun_Px_100k_v2_job0.root", rootfile_dir.Data(), kine, sbsfieldscale);
	// 	neutron_infile = new TFile(proton_infile_name.Data(), "READ");		
	}
	if( Nucleon_Select == "both" ){
		h_n_det_eff = new TH1D("h_n_det_eff", "Neutron detection efficiency", nBins, p_min, p_max);
		h_p_det_eff = new TH1D("h_p_det_eff", "Proton detection efficiency", nBins, p_min, p_max);
	}

	outfile_name = Form("rootfiles/mc_gun_det_eff_SBS%i.root", kine);
	outfile = new TFile(outfile_name.Data(), "RECREATE");

	//Histogram creations
	if( Nucleon_Select == "proton" || Nucleon_Select == "both" ){
		h_dx_p = new TH1D("h_dx_p", Form("Proton dx - SBS%i; x_{HCal}-x_{Expect} (m)", kine), 400, -2, 2);
		h_dy_p = new TH1D("h_dy_p", Form("Proton dy - SBS%i;y_{HCal}-y_{expect} (m)", kine), 400, 3.8, 7.8);		
		h_HCalE_p = new TH1D("h_HCalE_p", Form("Proton HCal Cluster Energy - SBS%i; HCal Cluster Energy (GeV)", kine), HCal_bins, HCal_min, HCal_max);

		h_det_eff_p = new TH1D("h_det_eff_p", Form("Proton Detection Efficiency - SBS%i; Nucleon Momentum (GeV/c); Efficiency (%%)", kine), 80, 1, 9);

		h2_xy_p = new TH2D("h2_xy_p", Form("Proton xy - SBS%i", kine), 12, -0.92837, 0.92837, 24, -2.35183, 1.45182);
		h2_xy_exp_p = new TH2D("h2_xy_exp_p", Form("Proton Expected xy - SBS%i", kine), 12, -0.92837, 0.92837, 24, -2.35183, 1.45182);
		h2_HCal_pN_p = new TH2D("h2_HCal_pN_p", Form("Proton Momentum - SBS%i; Momentum (GeV/c); HCal Cluster Energy (GeV)", kine), 80, 1, 9, HCal_bins, HCal_min, HCal_max);
		h2_HCal_pN_cut_p = new TH2D("h2_HCal_pN_cut_p", Form("Proton momentum with threshold = HCal E_{mean} / (4.0 each bin) - SBS%i; Proton Momentum (GeV/c); HCal Cluster Energy (GeV)", kine), 80, 1, 9, HCal_bins, HCal_min, HCal_max);

		tprof_h2_HCalE_pN_p = new TProfile("tprof_h2_HCalE_pN_p", Form("Proton Momentum Energy Profile - SBS%i; Proton Momentum (GeV/c); HCal Cluster Energy (GeV)", kine), pN_bins, pN_min, pN_max, HCal_min, HCal_max);
	}

	if( Nucleon_Select == "neutron" || Nucleon_Select == "both" ){
		h_dx_n = new TH1D("h_dx_n","Neutron dx;x_{HCal}-x_{expect} (m)", 400, -2, 2);
		h_dy_n = new TH1D("h_dy_n","Neutron dy;y_{HCal}-y_{expect} (m)", 400, 3.8, 7.8);

		h_HCalE_n = new TH1D("h_HCalE_n", Form("Neutron HCal Cluster Energy - SBS%i; HCal Cluster Energy (GeV)", kine), HCal_bins, HCal_min, HCal_max);

		h_det_eff_n = new TH1D("h_det_eff_n", Form("Neutron Detection Efficiency - SBS%i; Proton Momentum (GeV/c); Efficiency (%%)", kine), 80, 1, 9);

		h2_xy_n = new TH2D("h2_xy_n", Form("Neutron xy - SBS%i", kine), 12, -0.92837, 0.92837, 24, -2.35183, 1.45182);
		h2_xy_exp_n = new TH2D("h2_xy_exp_n", Form("Neutron Expected xy - SBS%i", kine), 12, -0.92837, 0.92837, 24, -2.35183, 1.45182);
		h2_HCal_pN_n = new TH2D("h2_HCal_pN_n", Form("Neutron Momentum - SBS%i; Momentum (GeV/c); HCal Cluster Energy (GeV)", kine), 80, 1, 9, HCal_bins, HCal_min, HCal_max);
		h2_HCal_pN_cut_n = new TH2D("h2_HCal_pN_cut_n", Form("Neutron momentum with threshold = HCal E_{mean} / (4.0 each bin) - SBS%i; Proton Momentum (GeV/c); HCal Cluster Energy (GeV)", kine), 80, 1, 9, HCal_bins, HCal_min, HCal_max);

		tprof_h2_HCalE_pN_n = new TProfile("tprof_h2_HCalE_pN_n", Form("Neutron Momentum Energy Profile - SBS%i; Nucleon Momentum (GeV/c); HCal Cluster Energy (GeV)", kine), pN_bins, pN_min, pN_max, HCal_min, HCal_max);

	}

	TVector3 HCal_xaxis(0, -1, 0);
	TVector3 HCal_zaxis( sin( -sbs_theta_rad ), 0, cos( -sbs_theta_rad ) );
	TVector3 HCal_yaxis = HCal_zaxis.Cross( HCal_xaxis ).Unit();

	vector<TVector3> HCal_axes = {HCal_xaxis, HCal_yaxis, HCal_zaxis};

	TVector3 HCal_origin = HCal_dist*HCal_zaxis;

	for(int ibin = 0; ibin < pN_bins; ibin++ ){


		if( ibin < (pN_bins/1.6) ){
			histo_max = HCal_max/1.2;
		}

		if( ibin < (pN_bins/4.0) ){
			histo_max = HCal_max/2.5;
		}

		if( ibin < (pN_bins/8.0) ){
			histo_max = HCal_max/4.0;
		}

		if( Nucleon_Select == "neutron" || Nucleon_Select == "both" ){
			h_HCalE_array_n[ibin] = new TH1D(Form("h_HCalE_array_n_%i", ibin), Form("h_HCalE_array_pbins_n_%i", ibin), pN_bins, 0.0, histo_max);
			h_HCalE_array_cut_n[ibin] = new TH1D(Form("h_HCalE_array_cut_n_%i", ibin), Form("h_HCalE_array_pbins_cut_n_%i", ibin), pN_bins, 0.0, histo_max);
		}

		if( Nucleon_Select == "proton" || Nucleon_Select == "both" ){
			h_HCalE_array_p[ ibin ] = new TH1D(Form("h_HCalE_array_p_%i", ibin), Form("h_HCalE_array_pbins_p_%i", ibin), pN_bins, 0.0, histo_max);
			h_HCalE_array_cut_p[ ibin ] = new TH1D(Form("h_HCalE_array_cut_p_%i", ibin), Form("h_HCalE_array_pbins_cut_p_%i", ibin), pN_bins, 0.0, histo_max);
		}

	}

//We can now loop over the two nucleons
//neutron: nucl = 0; proton: nucl = 1

	for(Int_t nucl = 0; nucl < 2; nucl++){

	//Forcing TChain to reset so that we don't use old values.
		TC = nullptr;
		TC = new TChain("T");

		if( nucl == 0 ){

			nucl_type = "neutron";

			for( size_t n_file = 0; n_file < neutron_infilenames_vec.size(); n_file++ ){
				TC->Add(neutron_infilenames_vec[n_file].Data());
			}
		}

		if( nucl == 1 ){

			nucl_type = "proton";

			for( size_t p_file = 0; p_file < proton_infilenames_vec.size(); p_file++ ){
				TC->Add(proton_infilenames_vec[p_file].Data());
			}
		}

			
	//Set up Tree Branches....
		TC->SetBranchStatus("*", 0);

		// HCal
		// TC->SetBranchStatus( "sbs.hcal.x", 1 );
		// TC->SetBranchStatus( "sbs.hcal.y", 1 );
		// TC->SetBranchStatus( "sbs.hcal.e", 1 );
		// TC->SetBranchStatus( "sbs.hcal.nclus", 1);

		//MC values for normalization
		TC->SetBranchStatus( "MC.mc_omega", 1 );
		TC->SetBranchStatus( "MC.mc_sigma", 1 );
		TC->SetBranchStatus( "MC.mc_nucl", 1);

		//MC normalization variables
		TC->SetBranchAddress( "MC.mc_omega", &mc_omega );
		TC->SetBranchAddress( "MC.mc_sigma", &mc_sigma );
		TC->SetBranchAddress( "MC.mc_fnucl", &mc_fnucl );	

		// HCal Cluster Tree Variables
		TC->SetBranchStatus( "sbs.hcal.e", 1);
		TC->SetBranchStatus( "sbs.hcal.x", 1);
		TC->SetBranchStatus( "sbs.hcal.y", 1);
		// TC->SetBranchStatus( "sbs.hcal.rowblk", 1);
		// TC->SetBranchStatus( "sbs.hcal.colblk", 1);
		// TC->SetBranchStatus( "sbs.hcal.tdctimeblk", 1);
		// TC->SetBranchStatus( "sbs.hcal.atimeblk", 1);
		// TC->SetBranchStatus( "sbs.hcal.idblk", 1);
		
		// HCal Cluster Tree Variables
		TC->SetBranchAddress( "sbs.hcal.e", &HCal_e);
		TC->SetBranchAddress( "sbs.hcal.x", &HCal_x);
		TC->SetBranchAddress( "sbs.hcal.y", &HCal_y);
		// TC->SetBranchAddress( "sbs.hcal.rowblk", &HCal_row);
		// TC->SetBranchAddress( "sbs.hcal.colblk", &HCal_col);
		// TC->SetBranchAddress( "sbs.hcal.atimeblk", &HCal_atime);
		// TC->SetBranchAddress( "sbs.hcal.tdctimeblk", &HCal_tdc);
		// TC->SetBranchAddress( "sbs.hcal.idblk", &HCal_id);

//--------------------------------------------------

		// MC Tree Variables
		TC->SetBranchStatus( "MC.mc_ep", 1);
		TC->SetBranchStatus( "MC.mc_epx", 1);
		TC->SetBranchStatus( "MC.mc_epy", 1);
		TC->SetBranchStatus( "MC.mc_epz", 1);
		TC->SetBranchStatus( "MC.mc_vx", 1);
		TC->SetBranchStatus( "MC.mc_vy", 1);
		TC->SetBranchStatus( "MC.mc_vz", 1);
		TC->SetBranchStatus( "MC.mc_nucl", 1);

		// MC Tree Variables
		TC->SetBranchAddress( "MC.mc_ep", &mc_p);
		TC->SetBranchAddress( "MC.mc_epx", &mc_px);
		TC->SetBranchAddress( "MC.mc_epy", &mc_py);
		TC->SetBranchAddress( "MC.mc_epz", &mc_pz);
		TC->SetBranchAddress( "MC.mc_vx", &mc_vx);
		TC->SetBranchAddress( "MC.mc_vy", &mc_vy);
		TC->SetBranchAddress( "MC.mc_vz", &mc_vz);
		TC->SetBranchAddress( "MC.mc_nucl", &mc_nucl);

		nEvents = TC->GetEntries();

		long nEvent = 0;

		while( TC->GetEntry(nEvent++) ){
			if( nEvent%5000 == 0 ){
				cout << "Thresh calcs... " << nucl_type.Data() << " Event " << nEvent << " of " << nEvents << endl;
				std::cout.flush();
			}

			if( nucl_type == "neutron" ){
				if( int(mc_fnucl) == 1 ){  
					continue;
				}
			}
			if( nucl_type == "proton" ){
				if( int(mc_fnucl) == 0 ){ 
					continue;
				}				
			}

		//Q-Vector
			TVector3 vertex(mc_vx, mc_vy, mc_vz); //Target Vertex vector in Hall Coordinates
			Double_t theta_nucleon_exp = acos( mc_pz/mc_p ); //Scattered Nucleon Theta
			Double_t phi_nucleon_exp = atan2( mc_py, mc_px ); //Scattered Nucleon Phi
			// Projected q vector:
			TVector3 pN_hat( sin(theta_nucleon_exp)*cos(phi_nucleon_exp), sin(theta_nucleon_exp)*sin(phi_nucleon_exp), cos(theta_nucleon_exp) );

		//Expected values at HCal
			double HCal_x_exp, HCal_y_exp;

		//Intersection of q ray at plane of HCal
			double s_intersect = (HCal_origin - vertex).Dot(HCal_zaxis)/(pN_hat.Dot(HCal_zaxis));

			TVector3 HCal_intersect = vertex + s_intersect*pN_hat;

			HCal_x_exp = (HCal_intersect - HCal_origin).Dot(HCal_xaxis);
			HCal_y_exp = (HCal_intersect - HCal_origin).Dot(HCal_yaxis);
			

			vector<double> HCal_xy_exp = {HCal_x_exp, HCal_y_exp};

		//Active area on HCal 
			//HCal_active_area[0] = top;
			//HCal_active_area[1] = bottom;
			//HCal_active_area[2] = right;
			//HCal_active_area[3] = left;

			vector<double> HCal_active_area = calc_HCal_active_area_MC();

			if( HCal_x_exp > HCal_active_area[0] && HCal_x_exp < HCal_active_area[1] &&
				HCal_y_exp > HCal_active_area[2] && HCal_y_exp < HCal_active_area[3] ){
				hit_on_HCal = true;
			}

			if( !hit_on_HCal ){
				if( nucl_type == "neutron" ){
					HCal_miss_cnt_n++;
				}
				if( nucl_type == "proton" ){
					HCal_miss_cnt_p++;			
				}

				continue;
			}

			if( nucl_type == "neutron" ){
				h_dx_n->Fill(HCal_x - HCal_x_exp);
				h_dy_n->Fill(HCal_y - HCal_y_exp);

				h2_xy_n->Fill(HCal_y, HCal_x);
				h2_xy_exp_n->Fill(HCal_y_exp, HCal_x_exp);

				h_HCalE_n->Fill(HCal_e);
				h2_HCal_pN_n->Fill(mc_p, HCal_e);
				tprof_h2_HCalE_pN_n->Fill(mc_p, HCal_e, 1.0);

				//Fill the histogram for various pN values
				// double bin_width = (pN_max - pN_min)/(pN_bins);
				double bin_width = h2_HCal_pN_n->GetXaxis()->GetBinWidth(1);
				int bin_i = int((mc_p - pN_min)/bin_width );

				if( bin_i >= pN_bins ){
					cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
					cout << "Error!!! Trying to access a bin that is out of range of pN_bins" << endl;
					cout << "Accessing neutron bin: " << bin_i << ".  Max bins = " << pN_bins << endl;
					cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
					continue;
				}
				// cout << "neutron bin_i is: " << bin_i << endl;
				// if( !(h_HCalE_array_n[bin_i]) ){
				// 	continue;
				// }
				h_HCalE_array_n[bin_i]->Fill(HCal_e);
			}

			if( nucl_type == "proton" ){
				h_dx_p->Fill(HCal_x - HCal_x_exp);
				h_dy_p->Fill(HCal_y - HCal_y_exp);

				h2_xy_p->Fill(HCal_y, HCal_x);
				h2_xy_exp_p->Fill(HCal_y_exp, HCal_x_exp);

				h_HCalE_p->Fill(HCal_e);
				h2_HCal_pN_p->Fill(mc_p, HCal_e);
				tprof_h2_HCalE_pN_p->Fill(mc_p, HCal_e, 1.0);

				//Fill the histogram for various pN values
				// double bin_width = (pN_max - pN_min)/(pN_bins);
				double bin_width = h2_HCal_pN_n->GetXaxis()->GetBinWidth(1);
				int bin_i = (mc_p - pN_min)/bin_width;

				if( bin_i >= pN_bins ){
					cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
					cout << "Error!!! Trying to access a bin that is out of range of pN_bins" << endl;
					cout << "Accessing proton bin: " << bin_i << ".  Max bins = " << pN_bins << endl;
					cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
					continue;
				}

				// if( !(h_HCalE_array_p[bin_i]) ){
				// 	continue;
				// }
				h_HCalE_array_p[bin_i]->Fill(HCal_e);					
			}	

		}

		if( nucl_type == "proton" ){
			tprof_h2_HCalE_pN_p->SetMarkerStyle(7);
			tprof_h2_HCalE_pN_p->SetMarkerColor(2);
		}

		if( nucl_type == "neutron" ){
			tprof_h2_HCalE_pN_n->SetMarkerStyle(7);
			tprof_h2_HCalE_pN_n->SetMarkerColor(2);
		}


//After that first loop we now need to determine bin-wise threshold (E_peak/4.0);
		cout << endl << "----------------------------" << endl;
		cout << "Calculating thresholds..." << endl;		
		cout << "----------------------------" << endl << endl;

		if( nucl_type == "neutron" ){
			tf_gaus_n = new TF1("tf_gaus_n", "gaus");

			for( int bin = 0; bin < pN_bins; bin++ ){
			//Pull E_mean from the histograms
				int max_bin = h_HCalE_array_n[bin]->GetMaximumBin();
				double max_bin_center = h_HCalE_array_n[bin]->GetXaxis()->GetBinCenter(max_bin);
				double max_count = h_HCalE_array_n[bin]->GetMaximum();
				double bin_width = h_HCalE_array_n[bin]->GetBinWidth(max_bin);
				double stdDev = h_HCalE_array_n[bin]->GetStdDev();

				// cout << "Neutron: " << endl;
				// cout << "max_bin: " << max_bin << ", max_bin_center: " << max_bin_center << endl;
				// cout << "bin_width: " << bin_width << ", stdDev: " << stdDev << endl;

			//Reject low energy peaks
				if( max_bin < 2 ){
					while( h_HCalE_array_n[bin]->GetBinContent(max_bin + 1) < h_HCalE_array_n[bin]->GetBinContent(max_bin) ||
						h_HCalE_array_n[bin]->GetBinContent(max_bin + 1) == h_HCalE_array_n[bin]->GetBinContent(max_bin) ){
						max_bin++;
						// if( max_bin >= pN_bins ){
						// 	break;
						// }	
					}
					h_HCalE_array_n[bin]->GetXaxis()->SetRange(max_bin + 1, h_HCalE_array_n[bin]->GetNbinsX());
					max_bin = h_HCalE_array_n[bin]->GetMaximumBin();
					max_bin_center = h_HCalE_array_n[bin]->GetXaxis()->GetBinCenter(max_bin);
					max_count = h_HCalE_array_n[bin]->GetMaximum();
					bin_width = h_HCalE_array_n[bin]->GetBinWidth(max_bin);
					stdDev = h_HCalE_array_n[bin]->GetStdDev();
				}

			//Fit the histograms
				double lower_binC = HCal_min + max_bin*bin_width - stdDev;
				double upper_binC = HCal_min + max_bin*bin_width + stdDev;

				tf_gaus_n->SetParameters(max_count, max_bin_center, stdDev);
				tf_gaus_n->SetRange(lower_binC, upper_binC);
				h_HCalE_array_n[bin]->Fit(tf_gaus_n, "RQ0");

			//Threshold = E_peak/4;
				threshold_n.push_back(tf_gaus_n->GetParameter(1)/4.0); 
			}
		}

		if( nucl_type == "proton" ){
			tf_gaus_p = new TF1("tf_gaus_p", "gaus");

			for( int bin = 0; bin < pN_bins; bin++ ){
			//Pull E_mean from the histograms
				int max_bin = h_HCalE_array_p[bin]->GetMaximumBin();
				double max_bin_center = h_HCalE_array_p[bin]->GetXaxis()->GetBinCenter(max_bin);
				double max_count = h_HCalE_array_p[bin]->GetMaximum();
				double bin_width = h_HCalE_array_p[bin]->GetBinWidth(max_bin);
				double stdDev = h_HCalE_array_p[bin]->GetStdDev();

				// cout << "Neutron: " << endl;
				// cout << "max_bin: " << max_bin << ", max_bin_center: " << max_bin_center << endl;
				// cout << "bin_width: " << bin_width << ", stdDev: " << stdDev << endl;

			//Reject low energy peaks
				if( max_bin < 2 ){
					while( h_HCalE_array_p[bin]->GetBinContent(max_bin + 1) < h_HCalE_array_p[bin]->GetBinContent(max_bin) ||
						h_HCalE_array_p[bin]->GetBinContent(max_bin + 1) == h_HCalE_array_p[bin]->GetBinContent(max_bin) ){
						max_bin++;
						// if( max_bin >= pN_bins ){
						// 	break;
						// }	
					}
					h_HCalE_array_p[bin]->GetXaxis()->SetRange(max_bin + 1, h_HCalE_array_p[bin]->GetNbinsX());
					max_bin = h_HCalE_array_p[bin]->GetMaximumBin();
					max_bin_center = h_HCalE_array_p[bin]->GetXaxis()->GetBinCenter(max_bin);
					max_count = h_HCalE_array_p[bin]->GetMaximum();
					bin_width = h_HCalE_array_p[bin]->GetBinWidth(max_bin);
					stdDev = h_HCalE_array_p[bin]->GetStdDev();
				}

			//Fit the histograms
				double lower_binC = HCal_min + max_bin*bin_width - stdDev;
				double upper_binC = HCal_min + max_bin*bin_width + stdDev;

				tf_gaus_p->SetParameters(max_count, max_bin_center, stdDev);
				tf_gaus_p->SetRange(lower_binC, upper_binC);
				h_HCalE_array_p[bin]->Fit(tf_gaus_p, "RQ0");

			//Threshold = E_peak/4;
				threshold_p.push_back(tf_gaus_p->GetParameter(1)/4.0); 
			}
		}

	//Loop again and now apply threshold cuts:
		cout << "---------------------------------" << endl;
		cout << "Applying threshold cuts..." << endl;
		cout << "---------------------------------" << endl;

		nEvent = 0;

		while( TC->GetEntry(nEvent++) ){
			if( nEvent%5000 == 0 ){
				cout << "Applying thresh cuts... " << nucl_type.Data() << " Event " << nEvent << " of " << nEvents << endl;
				std::cout.flush();
			}

			if( nucl_type == "proton" ){
				if( int(mc_fnucl) == 0 ){ 
					continue;
				}				
			}

	//Kinematics Again
		//Q-Vector
			TVector3 vertex(mc_vx, mc_vy, mc_vz); //Target Vertex vector in Hall Coordinates
			Double_t theta_nucleon_exp = acos( mc_pz/mc_p ); //Scattered Nucleon Theta
			Double_t phi_nucleon_exp = atan2( mc_py, mc_px ); //Scattered Nucleon Phi
			// Projected q vector:
			TVector3 pN_hat( sin(theta_nucleon_exp)*cos(phi_nucleon_exp), sin(theta_nucleon_exp)*sin(phi_nucleon_exp), cos(theta_nucleon_exp) );

		//Expected values at HCal
			double HCal_x_exp, HCal_y_exp;

		//Intersection of q ray at plane of HCal
			double s_intersect = (HCal_origin - vertex).Dot(HCal_zaxis)/(pN_hat.Dot(HCal_zaxis));

			TVector3 HCal_intersect = vertex + s_intersect*pN_hat;

			HCal_x_exp = (HCal_intersect - HCal_origin).Dot(HCal_xaxis);
			HCal_y_exp = (HCal_intersect - HCal_origin).Dot(HCal_yaxis);
			

			vector<double> HCal_xy_exp = {HCal_x_exp, HCal_y_exp};

		//Active area on HCal 
			//HCal_active_area[0] = top;
			//HCal_active_area[1] = bottom;
			//HCal_active_area[2] = right;
			//HCal_active_area[3] = left;

			vector<double> HCal_active_area = calc_HCal_active_area_MC();

			if( HCal_x_exp > HCal_active_area[0] && HCal_x_exp < HCal_active_area[1] &&
				HCal_y_exp > HCal_active_area[2] && HCal_y_exp < HCal_active_area[3] ){
				hit_on_HCal = true;
			}

			if( nucl_type == "neutron" ){
				
				if( int(mc_fnucl) == 1 ){  
					continue;
				}

			//Fill histograms for various nucleon momentum values.....
				// double bin_width = (pN_max - pN_min)/pN_bins;
				double bin_width = h2_HCal_pN_n->GetXaxis()->GetBinWidth(1);
				int bin = (mc_p - pN_min)/bin_width;

				//Apply threshold cut
				if( HCal_e > threshold_n[bin] ){
					h_HCalE_array_cut_n[bin]->Fill(HCal_e);

					double bin_width_y = (HCal_max - HCal_min)/HCal_bins;
					int bin_y = (HCal_e - HCal_min)/bin_width_y;

					h2_HCal_pN_cut_n->SetBinContent(bin, bin_y, 1);
				}

				if( mc_p > ( floor(p_central_sbs4*10.0)/10.0 ) && mc_p < ( ceil(p_central_sbs4*10.0)/10.0 ) ){
					p_sbs4_neut.push_back(mc_p);
				}
				if( mc_p > ( floor(p_central_sbs7*10.0)/10.0 ) && mc_p < ( ceil(p_central_sbs7*10.0)/10.0 ) ){
					p_sbs7_neut.push_back(mc_p);
				}
				if( mc_p > ( floor(p_central_sbs8*10.0)/10.0 ) && mc_p < ( ceil(p_central_sbs8*10.0)/10.0 ) ){
					p_sbs8_neut.push_back(mc_p);
				}
				if( mc_p > ( floor(p_central_sbs9*10.0)/10.0 ) && mc_p < ( ceil(p_central_sbs9*10.0)/10.0 ) ){
					p_sbs9_neut.push_back(mc_p);
				}
				if( mc_p > ( floor(p_central_sbs14*10.0)/10.0 ) && mc_p < ( ceil(p_central_sbs14*10.0)/10.0 ) ){
					p_sbs14_neut.push_back(mc_p);
				}

			}

			if( nucl_type == "proton" ){
				
				if( int(mc_fnucl) == 0 ){  
					continue;
				}

			//Fill histograms for various nucleon momentum values.....
				// double bin_width = (pN_max - pN_min)/pN_bins;
				double bin_width = h2_HCal_pN_p->GetXaxis()->GetBinWidth(1);
				int bin = (mc_p - pN_min)/bin_width;

				//Apply threshold cut
				if( HCal_e > threshold_p[bin] ){
					h_HCalE_array_cut_p[bin]->Fill(HCal_e);

					double bin_width_y = (HCal_max - HCal_min)/HCal_bins;
					int bin_y = (HCal_e - HCal_min)/bin_width_y;

					h2_HCal_pN_cut_p->SetBinContent(bin, bin_y, 1);
				}

				if( mc_p > ( floor(p_central_sbs4*10.0)/10.0 ) && mc_p < ( ceil(p_central_sbs4*10.0)/10.0 ) ){
					p_sbs4_prot.push_back(mc_p);
				}
				if( mc_p > ( floor(p_central_sbs7*10.0)/10.0 ) && mc_p < ( ceil(p_central_sbs7*10.0)/10.0 ) ){
					p_sbs7_prot.push_back(mc_p);
				}
				if( mc_p > ( floor(p_central_sbs8*10.0)/10.0 ) && mc_p < ( ceil(p_central_sbs8*10.0)/10.0 ) ){
					p_sbs8_prot.push_back(mc_p);
				}
				if( mc_p > ( floor(p_central_sbs9*10.0)/10.0 ) && mc_p < ( ceil(p_central_sbs9*10.0)/10.0 ) ){
					p_sbs9_prot.push_back(mc_p);
				}
				if( mc_p > ( floor(p_central_sbs14*10.0)/10.0 ) && mc_p < ( ceil(p_central_sbs14*10.0)/10.0 ) ){
					p_sbs14_prot.push_back(mc_p);
				}					
			}

		}

		// p_sbs4_neut_min = *min_element(p_sbs4_neut.begin(), p_sbs4_neut.end());
		// p_sbs4_neut_max = *max_element(p_sbs4_neut.begin(), p_sbs4_neut.end());

		// p_sbs4_prot_min = *min_element(p_sbs4_prot.begin(), p_sbs4_prot.end());
		// p_sbs4_prot_max = *max_element(p_sbs4_prot.begin(), p_sbs4_prot.end());

		// p_sbs7_neut_min = *min_element(p_sbs7_neut.begin(), p_sbs7_neut.end());
		// p_sbs7_neut_max = *max_element(p_sbs7_neut.begin(), p_sbs7_neut.end());

		// p_sbs7_prot_min = *min_element(p_sbs7_prot.begin(), p_sbs7_prot.end());
		// p_sbs7_prot_max = *max_element(p_sbs7_prot.begin(), p_sbs7_prot.end());

		// p_sbs8_neut_min = *min_element(p_sbs8_neut.begin(), p_sbs8_neut.end());
		// p_sbs8_neut_max = *max_element(p_sbs8_neut.begin(), p_sbs8_neut.end());

		// p_sbs8_prot_min = *min_element(p_sbs8_prot.begin(), p_sbs8_prot.end());
		// p_sbs8_prot_max = *max_element(p_sbs8_prot.begin(), p_sbs8_prot.end());

		// p_sbs9_neut_min = *min_element(p_sbs9_neut.begin(), p_sbs9_neut.end());
		// p_sbs9_neut_max = *max_element(p_sbs9_neut.begin(), p_sbs9_neut.end());

		// p_sbs9_prot_min = *min_element(p_sbs9_prot.begin(), p_sbs9_prot.end());
		// p_sbs9_prot_max = *max_element(p_sbs9_prot.begin(), p_sbs9_prot.end());				

		// p_sbs14_neut_min = *min_element(p_sbs14_neut.begin(), p_sbs14_neut.end());
		// p_sbs14_neut_max = *max_element(p_sbs14_neut.begin(), p_sbs14_neut.end());

		// p_sbs14_prot_min = *min_element(p_sbs14_prot.begin(), p_sbs14_prot.end());
		// p_sbs14_prot_max = *max_element(p_sbs14_prot.begin(), p_sbs14_prot.end());

		cout << "------------------------------------------------" << endl;
		cout << "------------------------------------------------" << endl;
		cout << "             Calculating Efficiencies" << endl;
		cout << "------------------------------------------------" << endl;
		cout << "------------------------------------------------" << endl << endl;

		if( nucl_type == "neutron" ){

			double gr_n_x[pN_bins];
			double gr_n_y[pN_bins];

			for(int n_bin = 0; n_bin < pN_bins; n_bin++ ){
				n_efficiency = 0.0;
				n_efficiency = (h_HCalE_array_cut_n[n_bin]->Integral()/h_HCalE_array_n[n_bin]->Integral())*100.0;
				double bin_width = h_HCalE_array_n[n_bin]->GetBinWidth(1);
				double bin_x = pN_min + (n_bin)*(bin_width) + bin_width/2.0;

				h_det_eff_n->SetBinContent(n_bin, n_efficiency);

				double x_val = ((p_max - p_min)/pN_bins)*n_bin + p_min;
				gr_n_x[n_bin] = x_val;
				gr_n_y[n_bin] = n_efficiency;

			}

			TCanvas *c_det_eff_n = new TCanvas("c_det_eff_n", "c_det_eff_n", 1200, 500);
			h_det_eff_n->GetYaxis()->SetRangeUser(det_eff_yaxis_min, 100.0);
			h_det_eff_n->Draw("same");

			TCanvas *c_gr_det_eff_n = new TCanvas("c_gr_det_eff_n", "c_gr_det_eff_n", 1200, 500);
			gr_det_eff_n = new TGraph(pN_bins, gr_n_x, gr_n_y);
			gr_det_eff_n->GetYaxis()->SetRangeUser(det_eff_yaxis_min, 100);
			gr_det_eff_n->GetXaxis()->SetRangeUser(1.2, 9.0);
			gr_det_eff_n->SetMarkerColor(4);
			gr_det_eff_n->Draw("AP0*");

			tf_det_eff_n = new TF1("tf_det_eff_n", "pol4", 1.2 , 9.0);
			// tf_det_eff_n->SetParLimits(1, -1000, -0.0001);
			// tf_det_eff_n->SetParLimits(2, 0.0001, 1000);
			tf_det_eff_n->SetLineColor(4);
			tf_det_eff_n->SetLineStyle(5);
			gr_det_eff_n->Fit("tf_det_eff_n", "RMSEF+");

			tf_det_eff_n->GetParameters( par_n );

			par_n_error[0] = tf_det_eff_n->GetParError(0);
			par_n_error[1] = tf_det_eff_n->GetParError(1);
			par_n_error[2] = tf_det_eff_n->GetParError(2);
			par_n_error[3] = tf_det_eff_n->GetParError(3);
			par_n_error[4] = tf_det_eff_n->GetParError(4);

			det_eff_err_n = calc_pol_error(par_n_error, 4);

			neut_det_eff_sbs4 = calc_det_eff( p_central_sbs4, par_n );
			neut_det_eff_sbs7 = calc_det_eff( p_central_sbs7, par_n );
			neut_det_eff_sbs8 = calc_det_eff( p_central_sbs8, par_n );
			neut_det_eff_sbs9 = calc_det_eff( p_central_sbs9, par_n );
			neut_det_eff_sbs14 = calc_det_eff( p_central_sbs14, par_n );


		}

		if( nucl_type == "proton"){

			double gr_p_x[pN_bins];
			double gr_p_y[pN_bins];

			for(int p_bin = 0; p_bin < pN_bins; p_bin++ ){
				p_efficiency = 0.0;
				p_efficiency = (h_HCalE_array_cut_p[p_bin]->Integral()/h_HCalE_array_p[p_bin]->Integral())*100.0;

				double bin_width = h_HCalE_array_p[p_bin]->GetBinWidth(1);
				double bin_x = p_min + (p_bin)*(bin_width) + bin_width/2.0;
				// double bin_x = pN_min + (p_bin)*(bin_width) + bin_width/2.0

				h_det_eff_p->SetBinContent(p_bin, p_efficiency);
								
				double x_val = ((p_max - p_min)/pN_bins)*p_bin + p_min;
				// cout << "x val: " << x_val << endl;
				gr_p_x[p_bin] = x_val;
				gr_p_y[p_bin] = p_efficiency;

			}

			TCanvas *c_det_eff_p = new TCanvas("c_det_eff_p", "c_det_eff_p", 1200, 500);
			h_det_eff_p->GetYaxis()->SetRangeUser(det_eff_yaxis_min, 100.0);
			h_det_eff_p->Draw("same");

			TCanvas *c_gr_det_ef_p = new TCanvas("c_gr_det_eff_p", "c_gr_det_eff_p", 1200, 500);
			gr_det_eff_p = new TGraph(pN_bins, gr_p_x, gr_p_y);
			gr_det_eff_p->GetYaxis()->SetRangeUser(det_eff_yaxis_min, 100);
			gr_det_eff_p->GetXaxis()->SetRangeUser(1.2, 9.0);
			gr_det_eff_p->SetMarkerColor(3);
			gr_det_eff_p->Draw("AP0*");

			tf_det_eff_p = new TF1("tf_det_eff_p", "pol4", 1.1, 9.0);
			// tf_det_eff_p->SetParLimits(1, -1000, -0.0001);
			// tf_det_eff_p->SetParLimits(2, 0.0001, 1000);
			tf_det_eff_p->SetLineStyle(5);
			tf_det_eff_p->SetLineColor(3);
			gr_det_eff_p->Fit("tf_det_eff_p", "RMSEF+");

			tf_det_eff_p->GetParameters( par_p );

			par_p_error[0] = tf_det_eff_p->GetParError(0);
			par_p_error[1] = tf_det_eff_p->GetParError(1);
			par_p_error[2] = tf_det_eff_p->GetParError(2);
			par_p_error[3] = tf_det_eff_p->GetParError(3);
			par_p_error[4] = tf_det_eff_p->GetParError(4);

			det_eff_err_p = calc_pol_error(par_p_error, 4);

			prot_det_eff_sbs4 = calc_det_eff( p_central_sbs4, par_p );
			prot_det_eff_sbs7 = calc_det_eff( p_central_sbs7, par_p );
			prot_det_eff_sbs8 = calc_det_eff( p_central_sbs8, par_p );
			prot_det_eff_sbs9 = calc_det_eff( p_central_sbs9, par_p );
			prot_det_eff_sbs14 = calc_det_eff( p_central_sbs14, par_p );



		}

		if( Nucleon_Select != "proton" && Nucleon_Select != "neutron" && Nucleon_Select == "both" ){
			TCanvas *c_det_eff_np = new TCanvas("c_det_eff_np", "c_det_eff_np", 1200, 500);

			h_det_eff_n->SetMarkerStyle(20);
			h_det_eff_n->SetMarkerSize(0.5);
			h_det_eff_n->SetMarkerColor(4);

			h_det_eff_p->SetMarkerStyle(21);
			h_det_eff_p->SetMarkerSize(0.5);
			h_det_eff_p->SetMarkerColor(3);

			h_det_eff_p->SetStats(false);
			h_det_eff_n->SetStats(false);

			h_det_eff_n->Draw("P0");
			tf_det_eff_n->Draw("same");

			h_det_eff_p->Draw("same+P0");
			if( nucl_type == "proton" ){
				tf_det_eff_p->Draw("same");
			}
			

			TLegend *tleg_pn_det_eff = new TLegend(0.12, 0.79, 0.22, 0.89);
			tleg_pn_det_eff->AddEntry(h_det_eff_p, "Proton");
			tleg_pn_det_eff->AddEntry(h_det_eff_n, "Neutron");
			tleg_pn_det_eff->Draw("SAME");

			TLine *tl_sbs4 = new TLine(2.35, det_eff_yaxis_min, 2.35, 100);
			tl_sbs4->SetLineStyle(3);
			tl_sbs4->Draw("SAME");

			TLine *tl_sbs7 = new TLine(6.20, det_eff_yaxis_min, 6.20, 100);
			tl_sbs7->SetLineStyle(3);
			tl_sbs7->Draw("SAME");

			TLine *tl_sbs8 = new TLine(3.22, det_eff_yaxis_min, 3.22, 100);
			tl_sbs8->SetLineStyle(3);
			tl_sbs8->Draw("SAME");

			TLine *tl_sbs9 = new TLine(3.21, det_eff_yaxis_min, 3.21, 100);
			tl_sbs9->SetLineStyle(3);
			tl_sbs9->Draw("SAME");

			TLine *tl_sbs14 = new TLine(4.84, det_eff_yaxis_min, 4.84, 100);
			tl_sbs14->SetLineStyle(3);
			tl_sbs14->Draw("SAME");

		//EFFICIENCIES FROM THE BIN VALUES

			tax_h_det_eff_n = h_det_eff_n->GetXaxis();
			tax_h_det_eff_p = h_det_eff_p->GetXaxis();

			prot_bin_det_eff_sbs4 = h_det_eff_p->GetBinContent(tax_h_det_eff_p->FindBin(p_central_sbs4));
			prot_bin_det_eff_sbs7 = h_det_eff_p->GetBinContent(tax_h_det_eff_p->FindBin(p_central_sbs7));
			prot_bin_det_eff_sbs8 = h_det_eff_p->GetBinContent(tax_h_det_eff_p->FindBin(p_central_sbs8));
			prot_bin_det_eff_sbs9 = h_det_eff_p->GetBinContent(tax_h_det_eff_p->FindBin(p_central_sbs9));
			prot_bin_det_eff_sbs14 = h_det_eff_p->GetBinContent(tax_h_det_eff_p->FindBin(p_central_sbs14));

			neut_bin_det_eff_sbs4 = h_det_eff_n->GetBinContent(tax_h_det_eff_n->FindBin(p_central_sbs4));
			neut_bin_det_eff_sbs7 = h_det_eff_n->GetBinContent(tax_h_det_eff_n->FindBin(p_central_sbs7));
			neut_bin_det_eff_sbs8 = h_det_eff_n->GetBinContent(tax_h_det_eff_n->FindBin(p_central_sbs8));
			neut_bin_det_eff_sbs9 = h_det_eff_n->GetBinContent(tax_h_det_eff_n->FindBin(p_central_sbs9));
			neut_bin_det_eff_sbs14 = h_det_eff_n->GetBinContent(tax_h_det_eff_n->FindBin(p_central_sbs14));
		//------------------------------------		

		//Labels for efficiency values	

			tpt_det_eff_bins = new TPaveText(6.4, 82, 7.475, 92.5);
			tpt_det_eff_bins->SetFillColor(kCyan);
			tpt_det_eff_bins->SetBorderSize(1);

			tpt_det_eff_bins->AddText("Efficiencies from p Bin");
			tpt_det_eff_bins->AddText("");
			tpt_det_eff_bins->AddText("SBS4:");
			tpt_det_eff_bins->AddText(Form("p: %0.1f%%, n: %0.1f%%", prot_bin_det_eff_sbs4, neut_bin_det_eff_sbs4));
			tpt_det_eff_bins->AddText("");

			tpt_det_eff_bins->AddText("SBS7:");
			tpt_det_eff_bins->AddText(Form("p: %0.1f%%, n: %0.1f%%", prot_bin_det_eff_sbs7, neut_bin_det_eff_sbs7));
			tpt_det_eff_bins->AddText("");

			tpt_det_eff_bins->AddText("SBS8:");
			tpt_det_eff_bins->AddText(Form("p: %0.1f%%, n: %0.1f%%", prot_bin_det_eff_sbs8, neut_bin_det_eff_sbs8));
			tpt_det_eff_bins->AddText("");

			tpt_det_eff_bins->AddText("SBS9:");
			tpt_det_eff_bins->AddText(Form("p: %0.1f%%, n: %0.1f%%", prot_bin_det_eff_sbs9, neut_bin_det_eff_sbs9));
			tpt_det_eff_bins->AddText("");

			tpt_det_eff_bins->AddText("SBS14:");
			tpt_det_eff_bins->AddText(Form("p: %0.1f%%, n: %0.1f%%", prot_bin_det_eff_sbs14, neut_bin_det_eff_sbs14));

			tpt_det_eff_bins->Draw("same");

		//-----------------------------------
 
			tpt_det_eff_fits = new TPaveText(7.525, 82, 8.6, 92.5);
			tpt_det_eff_fits->SetFillColor(kAzure+10);
			tpt_det_eff_fits->SetBorderSize(1);

			tpt_det_eff_fits->AddText("Efficiencies from Fit");
			tpt_det_eff_fits->AddText("");
			tpt_det_eff_fits->AddText("SBS4:");
			tpt_det_eff_fits->AddText(Form("p: %0.1f%%, n: %0.1f%%", prot_det_eff_sbs4, neut_det_eff_sbs4));
			tpt_det_eff_fits->AddText("");

			tpt_det_eff_fits->AddText("SBS7:");
			tpt_det_eff_fits->AddText(Form("p: %0.1f%%, n: %0.1f%%", prot_det_eff_sbs7, neut_det_eff_sbs7));
			tpt_det_eff_fits->AddText("");

			tpt_det_eff_fits->AddText("SBS8:");
			tpt_det_eff_fits->AddText(Form("p: %0.1f%%, n: %0.1f%%", prot_det_eff_sbs8, neut_det_eff_sbs8));
			tpt_det_eff_fits->AddText("");

			tpt_det_eff_fits->AddText("SBS9:");
			tpt_det_eff_fits->AddText(Form("p: %0.1f%%, n: %0.1f%%", prot_det_eff_sbs9, neut_det_eff_sbs9));
			tpt_det_eff_fits->AddText("");

			tpt_det_eff_fits->AddText("SBS14:");
			tpt_det_eff_fits->AddText(Form("p: %0.1f%%, n: %0.1f%%", prot_det_eff_sbs14, neut_det_eff_sbs14));

			tpt_det_eff_fits->Draw("same");


			TText *txt_sbs4 = new TText(p_central_sbs4 - 0.001, 82.5, "SBS4");
			txt_sbs4->SetTextAngle(90);
			txt_sbs4->SetTextSize(0.035f);
			txt_sbs4->Draw("same");

			TText *txt_sbs7 = new TText(p_central_sbs7 - 0.001, 82.5, "SBS7");
			txt_sbs7->SetTextAngle(90);
			txt_sbs7->SetTextSize(0.035f);
			txt_sbs7->Draw("same");

			TText *txt_sbs8 = new TText(p_central_sbs8 - 0.05, 82.5, "SBS8");
			txt_sbs8->SetTextAngle(90);
			txt_sbs8->SetTextSize(0.035f);
			txt_sbs8->Draw("same");

			TText *txt_sbs9 = new TText(p_central_sbs9 + 0.15, 82.5, "SBS9");
			txt_sbs9->SetTextAngle(90);
			txt_sbs9->SetTextSize(0.035f);
			txt_sbs9->Draw("same");

			TText *txt_sbs14 = new TText(p_central_sbs14 - 0.001, 82.5, "SBS14");
			txt_sbs14->SetTextAngle(90);
			txt_sbs14->SetTextSize(0.035f);
			txt_sbs14->Draw("same");


			TPad *tpad_det_eff = new TPad("tpad_det_eff", "Pad to hold the title for MC Det. Eff.", 0.04, 0.92, 0.96, 1.0);
			tpad_det_eff->Draw();
			tpad_det_eff->cd();

			TText *ttxt_det_eff_title = new TText(0.3, 0.04, "MC HCal Detection Efficiencies");
			ttxt_det_eff_title->SetTextSize(0.8f);
			ttxt_det_eff_title->Draw();
			
		}

		cout << "Writing to outfile...." << endl;
		outfile->Write();

		if( nucl_type == "neutron" ){
			HCal_hit_cnt_n = nEvents - HCal_miss_cnt_p;
		}

		if( nucl_type == "proton" ){
			HCal_hit_cnt_p = nEvents- HCal_miss_cnt_p;
		}

//END of NUCLEONS loop
	}

	
	cout << "-------------------------------------------" << endl;
	cout << "                Analysis Finished     " << endl;
	cout << "Total Entries: " << nEvents << endl;
	cout << "Proton hits on HCal: " << HCal_hit_cnt_p << endl;
	cout << "Neutron hits on HCal: " << HCal_hit_cnt_n << endl;

	auto total_time_end = high_resolution_clock::now();
	auto total_time_duration = duration_cast<minutes>(total_time_end - total_time_start);
	cout << endl << "---------------------------------------------------" << endl;
	cout << "Total time for analysis: " << total_time_duration.count() << " minutes. " << endl;
	// cout << endl << "Outfile: " << outfile_name.Data() << endl;
}