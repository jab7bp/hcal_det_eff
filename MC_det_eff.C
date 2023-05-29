//Script to calculate the detection efficiency for neutrons and protons for MC output
//Has calibrate and normal run modes
//First need to calibrate so that we can pull mean values for HCalE and p_Nucleon
//Written by J. Boyd - May 18, 2023

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

double par_n[3], par_p[5];
double par_n_error[3], par_p_error[5];
double p_det_eff_err = 0.0;
double n_det_eff_err = 0.0;

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

Double_t fit_gaus(Double_t * x, Double_t *par){

	Double_t g = 0.0;

	g = par[0]*exp((-0.5)*pow(((x[0] -  par[1])/par[2]),2));
	
	return g;
}

Double_t calc_det_eff(TString nucleon, int kine, double *par ){

	double det_eff = 0.0;
	double p_central = 0.0;

	switch( kine ){
		case 4: p_central = 2.35; break;
		case 7: p_central = 6.20; break;
		case 8: p_central = 3.22; break;
		case 9: p_central = 3.21; break;
		case 14: p_central = 4.84; break;
		default: 
			p_central = 0.0;
			det_eff = 0.0; break;
	}


	if( nucleon == "neutron" ){
		det_eff = par[0] + (par[1])/( pow( par[2], p_central ));
	}
	if( nucleon == "proton" ){
		det_eff = par[0] + par[1]*p_central + par[2]*pow( p_central, 2 ) + par[3]*pow( p_central, 3 ) + par[4]*pow( p_central, 4 );
	}
	return det_eff;

}

bool calibrate = false;

int kine = 4;
int sbsfieldscale = 0;
int num_jobs;

//Select Nucleon to calculate from
//"neutron", "proton", "both"
TString Nucleon_Select = "both";

//Binning values
const Int_t nBins = 150; //Num bins for HCal E vs. Nucleon Mom (p) for 1000 evt/bin
Double_t HCalE_pN_mean[nBins] = {0.0}; //Proton Energy on HCal (obtained from mean values of HCal_E vs p_Nucleon)
Double_t HCalE_nN_mean[nBins] = {0.0}; //Neutron Energy on HCal (obtained from mean values of HCal_E vs p_Nucleon)

///////////////////////////////////////////
const Double_t p_min = 1.0; //Min. Nucleon momentum, GeV
const Double_t p_max = 9.0; //Max. Nucleon momentum, GeV
const Double_t E_min = 0.0; //Min. Nucleon Energy on HCal, GeV
const Double_t E_max = 2.0; //Max. Nucleon Energy on HCal, GeV

Double_t p_step = ( p_max - p_min )/nBins; //Number of p stored per bin

//Filenames and related
TString outfile_name, rootfile_dir, proton_infile_name, neutron_infile_name;
vector<TString> proton_infilenames_vec, neutron_infilenames_vec;
TFile *outfile, *proton_infile, *neutron_infile;
TChain *TC;

ofstream proton_E_outfile, neutron_E_outfile;
TString proton_E_outfilename = Form("datafiles/proton_E_SBS%i_mag%i.dat", kine, sbsfieldscale);
TString neutron_E_outfilename = Form("datafiles/neutron_E_SBS%i_mag%i.dat", kine, sbsfieldscale);

//----Fit values
TF1 *p_gaus_fit, *n_gaus_fit;
TF1 *tf_n_det_eff, *tf_p_det_eff;

TCanvas *c_p_Emean, *c_n_Emean;
//Proton fit parameters
Double_t p_p0 = 0.0;
Double_t p_p1 = 0.114; 

//Neutron fit parameters
Double_t n_p0 = 0.0;
Double_t n_p1 = 0.116; 

Double_t t_fac = 3.0;

//Kinematic Variables
Double_t HCal_theta = lookup_HCal_angle_by_kine(kine, "rad");
Double_t HCal_dist = lookup_HCal_dist_by_kine(kine);

TVector3 HCal_zaxis( sin(-HCal_theta ), 0, cos(-HCal_theta) );
TVector3 HCal_xaxis( 0, -1, 0 );
TVector3 HCal_yaxis = HCal_zaxis.Cross(HCal_xaxis).Unit();

TVector3 HCal_origin = HCal_dist*HCal_zaxis;

double p_central_sbs4 = 2.35;
double p_central_sbs7 = 6.20;
double p_central_sbs8 = 3.22;
double p_central_sbs9 = 3.21;
double p_central_sbs14 = 4.84;

// TVector3 HCal_origin = HCal_dist*HCal_zaxis + hcalheight*HCal_xaxis;

Double_t set_param[3];
Double_t fit_len, fit_height;

TH1D *p_bin_slices_prot[nBins], *p_bin_slices_neut[nBins];

//Histograms
TH1D *h_dx_p_tr, *h_dy_p_tr, *h_dx_n_tr, *h_dy_n_tr;
TH1D *h_dx_p_ang, *h_dy_p_ang, *h_dx_n_ang, *h_dy_n_ang;
TH1D *h_x_exp, *h_y_exp;

TH2D *h_EdepvP_p, *h_EdepvP_n, *h_EdepvP_p_Ecut, *h_EdepvP_n_Ecut;
TH2D *h_dx_vp_p, *h_dy_vp_p, *h_dx_vp_n, *h_dy_vp_n;
TH2D *h_EdepvP_nucl_n, *h_EdepvP_nucl_p;

//TGraphs
auto mg_nucl = new TMultiGraph();
auto gr_prot = new TGraph();
auto gr_neut = new TGraph();
auto *gr_n_det_eff = new TGraph();
auto *gr_p_det_eff = new TGraph();

TH1D *h_n_det_eff, *h_p_det_eff;

//Counter bins for Nucleons
Int_t p_hit[nBins] = {0};
Int_t p_total[nBins] = {0};
Int_t n_hit[nBins] = {0};
Int_t n_total[nBins] = {0};
TString nucl_type;

Double_t E_threshold_p = 0.0;
Double_t E_threshold_n = 0.0;

//Efficiency calculations:
double prot_det_eff_sbs4 = 0.0, prot_det_eff_sbs7 = 0.0, prot_det_eff_sbs8 = 0.0, prot_det_eff_sbs9 = 0.0, prot_det_eff_sbs14 = 0.0;
double neut_det_eff_sbs4 = 0.0, neut_det_eff_sbs7 = 0.0, neut_det_eff_sbs8 = 0.0, neut_det_eff_sbs9 = 0.0, neut_det_eff_sbs14 = 0.0;

//TBranch variables
//HCal
Double_t HCal_id, HCal_e, HCal_x, HCal_y, HCal_row, HCal_col, HCal_tdc, HCal_atime;
Double_t mc_omega, mc_sigma, mc_fnucl;

//MC
Double_t mc_p, mc_px, mc_py, mc_pz, mc_vx, mc_vy, mc_vz, mc_nucl, mc_posx, mc_posy;

long nEvents = 0;

//calibrate = 0 --> NORMAL: Calculate Efficiency ---- Requires the mean values of HCal_E vs p_Nucleon
//calibrate = 1 --> Mode to calculate mean values of HCalE vs p_Nucleon
void MC_det_eff(){ 

	auto total_time_start = high_resolution_clock::now();
	TStopwatch *StopWatch = new TStopwatch();

	gStyle->SetPalette(55);
	cout << "--------------------------------------" << endl;
	cout << "Analysis started for MC-based Detection Efficiencies. " << endl;
	cout << "--------------------------------------" << endl;

	cout << "Run parameters: " << endl;
	cout << "Kinematic: SBS" << kine << endl;
	cout << "SBS Field Scale: " << sbsfieldscale << endl;

	rootfile_dir = "/lustre19/expphy/volatile/halla/sbs/jboyd/simulation/out_dir/MC_REPLAY_OUT_DIR";

	num_jobs = 6;


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

	outfile_name = Form("rootfiles/mc_det_eff_SBS%i_tfac%0.0f.root", kine, t_fac);
	outfile = new TFile(outfile_name.Data(), "RECREATE");
	//Forcing TChain to reset for memory allocation
	TC = nullptr;

	//Histogram creations
	h_dx_p_tr = new TH1D("h_dx_p_tr", "proton dx (sd track); x_{HCal}-x_{Expect} (m)", 400, -2, 2);
	h_dy_p_tr = new TH1D("h_dy_p_tr","proton dy (sd track);y_{HCal}-y_{expect} (m)", 400, 3.8, 7.8);
	h_dx_n_tr = new TH1D("h_dx_n_tr","neutron dx (sd track);x_{HCal}-x_{expect} (m)", 400, -2, 2);
	h_dy_n_tr = new TH1D("h_dy_n_tr","neutron dy (sd track);y_{HCal}-y_{expect} (m)", 400, 3.8, 7.8);

	h_dx_p_ang = new TH1D("h_dx_p_ang","proton dx (angles);x_{HCal}-x_{expect} (m)", 400, -2, 2);
	h_dy_p_ang = new TH1D("h_dy_p_ang","proton dy (angles);y_{HCal}-y_{expect} (m)", 400, -2, 2);
	h_dx_n_ang = new TH1D("h_dx_n_ang","neutron dx (angles);x_{HCal}-x_{expect} (m)", 400, -2, 2);
	h_dy_n_ang = new TH1D("h_dy_n_ang","neutron dy (angles);y_{HCal}-y_{expect} (m)", 400, -2, 2);
	h_x_exp = new TH1D("h_x_exp","x exp (angles);x_{expect} (m)", 400, -2, 2);
	h_y_exp = new TH1D("h_y_exp","y exp (angles);y_{expect} (m)", 400, -2, 2);

	h_EdepvP_p = new TH2D("h_EdepvP_p","HCal E dep vs proton momentum; p_{proton} (GeV); E_{hcal} (GeV)", nBins, 1, 9, 200, E_min, E_max);
	h_EdepvP_n = new TH2D("h_EdepvP_n","HCal E dep vs neutron momentum; p_{neutron} (GeV); E_{hcal} (GeV)", nBins, 1, 9, 200, E_min, E_max);
	h_EdepvP_p_Ecut = new TH2D("h_EdepvP_p_Ecut","HCal E dep vs proton momentum, Mean E / 4 cut; p_{proton} (GeV); E_{hcal} (GeV)", nBins, 1, 9, 200, E_min, E_max);
	h_EdepvP_n_Ecut = new TH2D("h_EdepvP_n_Ecut","HCal E dep vs neutron momentum, Mean E / 4 cut; p_{neutron} (GeV); E_{hcal} (GeV)", nBins, 1, 9, 200, E_min, E_max);

	h_dx_vp_p = new TH2D("h_dx_vp_p","dx vs proton p; p_{p} (GeV); x_{HCal}-x_{expect} (m)", nBins, 1, 9, 400, -2, 2);
	h_dy_vp_p = new TH2D("h_dy_vp_p","dy vs proton p; p_{p} (GeV); y_{HCal}-y_{expect} (m)", nBins, 1, 9, 400, -2, 2);
	h_dx_vp_n = new TH2D("h_dx_vp_n","dx vs neutron p; p_{n} (GeV); x_{HCal}-x_{expect} (m)", nBins, 1, 9, 400, -2, 2);
	h_dy_vp_n = new TH2D("h_dy_vp_n","dy vs neutron p; p_{n} (GeV); y_{HCal}-y_{expect} (m)", nBins, 1, 9, 400, -2, 2);

	//We can now loop over the two nucleons
	//neutron: nucl = 0; proton: nucl = 1

	for(Int_t nucl = 0; nucl < 2; nucl++){
		TC = new TChain("T");

		if( nucl == 0 ){
			nucl_type = "neutron";

			for( size_t n_file = 0; n_file < neutron_infilenames_vec.size(); n_file++ ){
				TC->Add(neutron_infilenames_vec[n_file].Data());
			}
			// TC->Add(neutron_infile->GetName());
		}
		if( nucl == 1 ){
			nucl_type = "proton";

			for( size_t p_file = 0; p_file < proton_infilenames_vec.size(); p_file++ ){
				TC->Add(proton_infilenames_vec[p_file].Data());
			}
			// TC->Add(proton_infile->GetName());
		}
		// else{
		// 	break;
		// }
		
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
		TC->SetBranchStatus( "sbs.hcal.rowblk", 1);
		TC->SetBranchStatus( "sbs.hcal.colblk", 1);
		TC->SetBranchStatus( "sbs.hcal.tdctimeblk", 1);
		TC->SetBranchStatus( "sbs.hcal.atimeblk", 1);
		TC->SetBranchStatus( "sbs.hcal.idblk", 1);
		
		// HCal Cluster Tree Variables
		TC->SetBranchAddress( "sbs.hcal.e", &HCal_e);
		TC->SetBranchAddress( "sbs.hcal.x", &HCal_x);
		TC->SetBranchAddress( "sbs.hcal.y", &HCal_y);
		TC->SetBranchAddress( "sbs.hcal.rowblk", &HCal_row);
		TC->SetBranchAddress( "sbs.hcal.colblk", &HCal_col);
		TC->SetBranchAddress( "sbs.hcal.atimeblk", &HCal_atime);
		TC->SetBranchAddress( "sbs.hcal.tdctimeblk", &HCal_tdc);
		TC->SetBranchAddress( "sbs.hcal.idblk", &HCal_id);

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
		TC->SetBranchStatus( "MC.sdtrack_posx", 1);
		TC->SetBranchStatus( "MC.sdtrack_posy", 1);
        	
		// MC Tree Variables
		TC->SetBranchAddress( "MC.mc_ep", &mc_p);
		TC->SetBranchAddress( "MC.mc_epx", &mc_px);
		TC->SetBranchAddress( "MC.mc_epy", &mc_py);
		TC->SetBranchAddress( "MC.mc_epz", &mc_pz);
		TC->SetBranchAddress( "MC.mc_vx", &mc_vx);
		TC->SetBranchAddress( "MC.mc_vy", &mc_vy);
		TC->SetBranchAddress( "MC.mc_vz", &mc_vz);
		TC->SetBranchAddress( "MC.mc_nucl", &mc_nucl);
		TC->SetBranchAddress( "MC.sdtrack_posx", &mc_posx);
		TC->SetBranchAddress( "MC.sdtrack_posy", &mc_posy);

		nEvents = TC->GetEntries();

		cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl << endl;
		cout << "Number of events to analyze: " << nEvents << endl;
		cout << endl << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;

		if( !calibrate ){
			cout << endl << " Running in normal mode " << endl << endl;

			cout << "Pulling mean E values from datafiles" << endl;

			if( Nucleon_Select == "proton" || Nucleon_Select == "both" ){
				cout << "----- Grabbing proton values... " << endl;
				ifstream proton_E_infile( proton_E_outfilename.Data() );
				if( !proton_E_infile ){
					cerr << endl << "------------------------------" << endl;
					cerr << "ERROR: missing proton E mean infile: " << proton_E_outfilename.Data() << endl;
					cerr << "------------------------------" << endl << endl;
				}
				Int_t p_E_iter = 0;
				Double_t p_E_val = 0.0;
				string p_line;

				while( getline( proton_E_infile, p_line ) ){
					if( p_line.at(0) == '#' ){
						continue;
					}
					stringstream stream( p_line );
					stream >> p_E_val;
					HCalE_pN_mean[p_E_iter] = p_E_val;
					p_E_iter++;
				}
			}

			if( Nucleon_Select == "neutron" || Nucleon_Select == "both" ){
				cout << "----- Grabbing neutron values... " << endl;
				ifstream neutron_E_infile( neutron_E_outfilename.Data() );
				if( !neutron_E_infile ){
					cerr << endl << "------------------------------" << endl;
					cerr << "ERROR: missing neutron E mean infile: " << neutron_E_outfilename.Data() << endl;
					cerr << "------------------------------" << endl << endl;
				}
				Int_t n_E_iter = 0;
				Double_t n_E_val = 0.0;
				string n_line;

				while( getline( neutron_E_infile, n_line ) ){
					if( n_line.at(0) == '#' ){
						continue;
					}
					stringstream stream( n_line );
					stream >> n_E_val;
					HCalE_nN_mean[n_E_iter] = n_E_val;
					n_E_iter++;
				}
			}
			cout << "Finished grabbing mean E values.... " << endl << endl;

		}

	//Start of EVENTS loop
		long nEvent = 0;
		while( TC->GetEntry(nEvent++) ){

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

			if( nEvent%100 == 0 ){
				cout << "Nucleon: " << nucl_type.Data() << ", MC event: " << nEvent << endl;
				cout.flush();
			}

			Int_t bin = (mc_p - p_min)/p_step;

			if( bin > 150 || bin < 0 ){
				cout << "Warning!!!! Nucleon momentum bin out of bounds at bin: " << bin << endl;
			}

			//sd Track Calculations
			Double_t dx = -HCal_x - mc_posy; //Hall coordinates in MC
			Double_t dy = HCal_y - mc_posx;

			//MC nucleon momentum calculations
			TVector3 vertex( mc_vx, mc_vy, mc_vz );
			Double_t theta_Nucleon_exp = acos( mc_pz/mc_p );
			Double_t phi_Nucleon_exp = atan2( mc_py, mc_px );
			TVector3 pN_hat( sin(theta_Nucleon_exp)*cos(phi_Nucleon_exp), sin(theta_Nucleon_exp)*sin(phi_Nucleon_exp), cos(theta_Nucleon_exp));;

		///HCal Intersect variables
			//Intersection point for a ray landing on HCal face coordinates
			Double_t scint_intersect = ( HCal_origin - vertex ).Dot( HCal_zaxis )/( pN_hat.Dot(HCal_zaxis));

			//Ray from the origin in hall to the hit point on HCal face
			TVector3 HCal_intersect = vertex + scint_intersect*pN_hat;

			Double_t x_expected_HCal = ( HCal_intersect - HCal_origin ).Dot( HCal_xaxis );
			Double_t y_expected_HCal = ( HCal_intersect - HCal_origin ).Dot( HCal_yaxis );

		//position of x and y expected on HCal
			vector<Double_t> xy_HCal_exp;

			xy_HCal_exp.push_back(x_expected_HCal);
			xy_HCal_exp.push_back(y_expected_HCal);

			h_x_exp->Fill(x_expected_HCal);
			h_y_exp->Fill(y_expected_HCal);

			Double_t dx_v2 = HCal_x - x_expected_HCal;
			Double_t dy_v2 = HCal_y - y_expected_HCal;

			if( nucl_type == "proton" ){
				p_total[bin]++;
				E_threshold_p = HCalE_pN_mean[bin]/t_fac;

				if( HCal_e > E_threshold_p ){
					p_hit[bin]++;
					h_EdepvP_p->Fill(mc_p, HCal_e);
				}

				h_dx_p_tr->Fill(dx);
				h_dy_p_tr->Fill(dy);
				h_dx_vp_p->Fill(mc_p, dx_v2);
				h_dy_vp_p->Fill(mc_p, dy_v2);
				h_dx_p_ang->Fill(dx_v2);
				h_dy_p_ang->Fill(dy_v2);
			}

			if( nucl_type == "neutron" ){
				n_total[bin]++;
				E_threshold_n = HCalE_nN_mean[bin]/t_fac;
				if( HCal_e > E_threshold_n){
					n_hit[bin]++;
					h_EdepvP_n->Fill(mc_p, HCal_e);
				}
				h_EdepvP_n->Fill(mc_p, HCal_e);
				h_dx_n_tr->Fill(dx);
				h_dy_n_tr->Fill(dy);
				h_dx_vp_n->Fill(mc_p, dx_v2);
				h_dy_vp_n->Fill(mc_p, dy_v2);
				h_dx_n_ang->Fill(dx_v2);
				h_dy_n_ang->Fill(dy_v2);
			}

		///End of EVENTS while loop
		}
		TC->Reset();

	}

	Double_t bin_p[nBins] = {0.0};
	Double_t E_mean_p[nBins] = {0.0};
	Double_t E_sig_p[nBins] = {0.0};
	Double_t E_mean_n[nBins] = {0.0};
	Double_t E_sig_n[nBins] = {0.0};
	Double_t bin_err_p[nBins] = {0.0};
	Double_t bin_err_n[nBins] = {0.0};


	if( calibrate ){
		cout << endl << "-----------------------------" << endl;
		cout << "      Running in Calibration Mode "   << endl;
		cout << "-----------------------------" << endl;

		//Need to grab info for each/either nucleon
		// nucl = 0 for neutron, nucl = 1 for proton
		for(Int_t nucl = 0; nucl < 2; nucl++ ){
			Double_t p0_p = 0.0;
			Double_t p1_p = 0.0;
			Double_t p0_n = 0.0;
			Double_t p1_n = 0.0;

			if( nucl == 0 ){
				nucl_type = "neutron";
			}
			if( nucl == 1 ){
				nucl_type = "proton";
			}
			if( nucl_type == "neutron"){
				p0_n = n_p0;
				p1_n = n_p1;

				h_EdepvP_nucl_n = (TH2D*)(h_EdepvP_n->Clone("h_EdepvP_nucl_n"));

				c_n_Emean = new TCanvas("c_n_Emean", "c_n_Emean", 600, 500);
			}
			if( nucl_type == "proton" ){
				p0_p = p_p0;
				p1_p = p_p1;

				h_EdepvP_nucl_p = (TH2D*)(h_EdepvP_p->Clone("h_EdepvP_nucl_p"));

				c_p_Emean = new TCanvas("c_p_Emean", "c_p_Emean", 600, 500);
			}

		//Plotting HCalE and p_Nucleon to get mean
			for(Int_t bin = 0; bin < nBins; bin++ ){
				Double_t p_N = bin*p_step + p_min;
				Double_t fit_p1_exp;

				if( nucl_type == "neutron" ){
					fit_p1_exp = p0_n - 0.1 + p1_n*p_N;
				} 

				if( nucl_type == "proton" ){
					fit_p1_exp = p0_p - 0.1 + p1_p*p_N;
				} 

				Double_t fit_p2_exp = 0.002*p_N + 0.327; //Guestimate
				Int_t slice_N = 0;

				if( nucl_type == "neutron" ){
					p_bin_slices_neut[bin] = h_EdepvP_nucl_n->ProjectionY(Form("p_bin_slices_neut_%d", bin+1), bin+1, bin+1);					
					slice_N = p_bin_slices_neut[bin]->GetEntries();
				}
				
				if( nucl_type == "proton" ){
					p_bin_slices_prot[bin] = h_EdepvP_nucl_p->ProjectionY(Form("p_bin_slices_prot_%d", bin+1), bin+1, bin+1);					
					slice_N = p_bin_slices_prot[bin]->GetEntries();
				}

				set_param[0] = 800;
				set_param[1] = fit_p1_exp;
				set_param[2] = fit_p2_exp;

				if( bin < 30 ){
					fit_len = fit_p1_exp - 1.0*fit_p2_exp;
					if( fit_len < 0 ){
						fit_len	= 0;
					}
					fit_height = fit_p1_exp + 1.0*fit_p2_exp;
				}
				else{
					fit_len = fit_p1_exp - 3.0*fit_p2_exp;
					if( fit_len < 0 ){
						fit_len = 0;
					}
					fit_height = fit_p1_exp + 3.0*fit_p2_exp;
				}

				if( nucl_type == "proton" ){
					c_p_Emean->cd();
					p_gaus_fit = new TF1("p_gaus_fit", fit_gaus, fit_len, fit_height, 3);
					p_gaus_fit->SetLineWidth(4);
					p_gaus_fit->SetParameter(0, set_param[0]);
					p_gaus_fit->SetParameter(1, set_param[1]);
					p_gaus_fit->SetParLimits(1, fit_len, fit_height);
					p_gaus_fit->SetParameter(2, set_param[2]);
					p_gaus_fit->SetParLimits(2, E_min, E_max);

					p_bin_slices_prot[bin]->Fit("p_gaus_fit", "RBM");
					p_bin_slices_prot[bin]->Draw();

					bin_p[bin] = p_N;

			      	E_mean_p[bin] = p_gaus_fit->GetParameter(1);
			      	E_sig_p[bin] = p_gaus_fit->GetParameter(2);
			      	p_bin_slices_prot[bin]->SetTitle(Form("Loop:%d Np:%f Nuc:%d Mean:%f Sigma:%f MeanExp:%f SigExp:%f",bin, p_N, nucl, E_mean_p[bin],E_sig_p[bin],fit_p1_exp,fit_p2_exp));    

				}

				if( nucl_type == "neutron" ){
					c_n_Emean->cd();

					n_gaus_fit = new TF1("n_gaus_fit", fit_gaus, fit_len, fit_height, 3);
					n_gaus_fit->SetLineWidth(4);
					n_gaus_fit->SetParameter(0, set_param[0]);
					n_gaus_fit->SetParameter(1, set_param[1]);
					n_gaus_fit->SetParLimits(1, fit_len, fit_height);
					n_gaus_fit->SetParameter(2, set_param[2]);
					n_gaus_fit->SetParLimits(2, E_min, E_max);

					p_bin_slices_neut[bin]->Fit("n_gaus_fit", "RBM");
					p_bin_slices_neut[bin]->Draw();
					bin_p[bin] = p_N;

			      	E_mean_n[bin] = n_gaus_fit->GetParameter(1);
			      	E_sig_n[bin] = n_gaus_fit->GetParameter(2);
			      	p_bin_slices_neut[bin]->SetTitle(Form("Loop:%d Np:%f Nuc:%d Mean:%f Sigma:%f MeanExp:%f SigExp:%f", bin, p_N, nucl, E_mean_n[bin], E_sig_n[bin], fit_p1_exp, fit_p2_exp));    

				}

		    //END OF LOOP OVER BINS
			}
		//END OF LOOP OVER NUCLEONS
		}

	//PLOTTING
		if( Nucleon_Select == "proton" || Nucleon_Select == "both"){
			// TCanvas *c_prot = new TCanvas("c_prot", "c_prot", 600, 500);
		    gr_prot = new TGraphErrors(nBins, bin_p,E_mean_p,bin_err_p,E_sig_p);
		    gr_prot->SetTitle("Proton");
		    gr_prot->SetMarkerColor(kRed);
		    gr_prot->SetMarkerStyle(33);
		    gr_prot->SetMarkerSize(2);
		    gr_prot->SetLineColor(kRed);
		    gr_prot->SetLineWidth(2);
		    // gr_prot->Draw("AP*");
		    gr_prot->Write();
		    mg_nucl->Add(gr_prot);			
		}

		if( Nucleon_Select == "neutron" || Nucleon_Select == "both"){
			// TCanvas *c_neut = new TCanvas("c_neut", "c_neut", 600, 500);
		    gr_neut = new TGraphErrors(nBins, bin_p,E_mean_p,bin_err_p,E_sig_p);
		    gr_neut->SetTitle("Neutron");
		    gr_neut->SetMarkerColor(kBlue);
		    gr_neut->SetMarkerStyle(34);
		    gr_neut->SetMarkerSize(2);
		    gr_neut->SetLineColor(kBlue);
		    gr_neut->SetLineWidth(2);
		    // gr_neut->Draw("AP*");
		    gr_neut->Write();
		    mg_nucl->Add(gr_neut);
		}

		TCanvas *c_mg_nucl = new TCanvas("c_mg_nucl", "c_mg_nucl", 600, 500);
	    mg_nucl->SetTitle("HCal Energy vs Nucleon Momentum");
	    mg_nucl->GetXaxis()->SetTitle("Nucleon momentum, p (GeV)");
	    mg_nucl->GetYaxis()->SetTitle("E_{HCal}");
	    mg_nucl->Draw("AP*");

	    mg_nucl->Write();

		if( Nucleon_Select == "proton" || Nucleon_Select == "both" ){
			proton_E_outfile.open( proton_E_outfilename.Data() );
			proton_E_outfile << "##HCal mean E from proton data Gaussian fit" << endl;

			for( Int_t bin = 0; bin < nBins; bin++ ){
				proton_E_outfile << E_mean_p[bin] << endl;
			}
			proton_E_outfile.close();

			cout << "Proton mean Energy written to datafiles/" << proton_E_outfilename.Data() << endl;

			cout << "E_mean_p: " << E_mean_p[nBins - 1] << endl;
		}

		if( Nucleon_Select == "neutron" || Nucleon_Select == "both" ){
			neutron_E_outfile.open( neutron_E_outfilename.Data() );
			neutron_E_outfile << "##HCal mean E from neutron data Gaussian fit" << endl;

			for( Int_t bin = 0; bin < nBins; bin++ ){
				neutron_E_outfile << E_mean_n[bin] << endl;
			}
			neutron_E_outfile.close();

			cout << "Neutron mean Energy written to datafiles/" << neutron_E_outfilename.Data() << endl;

			cout << "E_mean_n: " << E_mean_n[nBins - 1] << endl;
		}	
	}

	Double_t proton_det_eff[nBins] = {0.0};
	Double_t neutron_det_eff[nBins] = {0.0};
	Double_t pn_ratio_det_eff[nBins] = {0.0};

	if( !calibrate ){
		for( Int_t bin = 0; bin < nBins; bin++ ){
			Double_t p = bin*p_step + p_min;
			bin_p[bin] = p;

			if( Nucleon_Select == "proton" || Nucleon_Select == "both" ){
				proton_det_eff[bin] = ( (double)p_hit[bin]/(double)p_total[bin] )*100.0;
				h_p_det_eff->SetBinContent(bin, proton_det_eff[bin]);
			}
			if( Nucleon_Select == "neutron" || Nucleon_Select == "both" ){
				neutron_det_eff[bin] = ( (double)n_hit[bin]/(double)n_total[bin] )*100.0;
				h_n_det_eff->SetBinContent(bin, neutron_det_eff[bin]);
			}
			if( Nucleon_Select == "both" ){
				pn_ratio_det_eff[bin] = neutron_det_eff[bin]/proton_det_eff[bin];
			}

		}
		

		//Efficiency graphs

		int det_eff_yaxis_min = 82;

		if( Nucleon_Select == "proton" || Nucleon_Select == "both"){
			TCanvas *c_p_det_eff = new TCanvas("c_p_det_eff", "c_p_det_eff", 1200, 500);

			double p_xaxis_min = 1.3;

			gr_p_det_eff = new TGraph(nBins, bin_p, proton_det_eff);
			gr_p_det_eff->SetTitle("MC HCal Proton Detection Efficiency");
			gr_p_det_eff->GetXaxis()->SetRangeUser(p_xaxis_min, 9);
			gr_p_det_eff->GetXaxis()->SetTitle("Proton momentum (GeV)");
			gr_p_det_eff->GetYaxis()->SetRangeUser(det_eff_yaxis_min, 100);
			gr_p_det_eff->GetYaxis()->SetTitle("Detection Efficieincy (%)");
			gr_p_det_eff->SetMarkerColor(51);
			gr_p_det_eff->Draw("AP*");

			tf_p_det_eff = new TF1("tf_p_det_eff", "pol4", p_xaxis_min, 9.0);
			tf_p_det_eff->SetLineColor(1);
			tf_p_det_eff->SetLineStyle(4);
			tf_p_det_eff->SetLineColor(51);
			gr_p_det_eff->Fit("tf_p_det_eff", "RMSEF+");

			tf_p_det_eff->GetParameters( par_p );

			par_p_error[0] = tf_p_det_eff->GetParError(0);
			par_p_error[1] = tf_p_det_eff->GetParError(1);
			par_p_error[2] = tf_p_det_eff->GetParError(2);
			par_p_error[3] = tf_p_det_eff->GetParError(3);
			par_p_error[4] = tf_p_det_eff->GetParError(4);

			p_det_eff_err = calc_pol_error(par_p_error, 5);

		}

		if( Nucleon_Select == "neutron" || Nucleon_Select == "both"){
			TCanvas *c_n_det_eff = new TCanvas("c_n_det_eff", "c_n_det_eff", 1200, 500);

			double n_xaxis_min = 1.4;

			gr_n_det_eff = new TGraph(nBins, bin_p, neutron_det_eff);
			gr_n_det_eff->SetTitle("MC HCal Neutron Detection Efficiency");
			gr_n_det_eff->GetXaxis()->SetRangeUser(n_xaxis_min, 9);
			gr_n_det_eff->GetXaxis()->SetTitle("Neutron momentum (GeV)");
			gr_n_det_eff->GetYaxis()->SetRangeUser(det_eff_yaxis_min, 100);
			gr_n_det_eff->GetYaxis()->SetTitle("Detection Efficieincy (%)");
			gr_n_det_eff->SetMarkerColor(8);
			gr_n_det_eff->Draw("AP*");

			tf_n_det_eff = new TF1("tf_n_det_eff", fit_det_eff, n_xaxis_min, 9.0, 3);
			tf_n_det_eff->SetParLimits(1, -1000, -0.0001);
			tf_n_det_eff->SetParLimits(2, 0.0001, 1000);
			tf_n_det_eff->SetLineColor(1);
			tf_n_det_eff->SetLineStyle(5);
			tf_n_det_eff->SetLineColor(8);
			gr_n_det_eff->Fit("tf_n_det_eff", "RMSEF+");

			tf_n_det_eff->GetParameters( par_n );

			par_n_error[0] = tf_n_det_eff->GetParError(0);
			par_n_error[1] = tf_n_det_eff->GetParError(1);
			par_n_error[2] = tf_n_det_eff->GetParError(2);

			n_det_eff_err = calc_pol_error(par_n_error, 3);
		
		}

		prot_det_eff_sbs4 = calc_det_eff( "proton", 4, par_p );
		prot_det_eff_sbs7 = calc_det_eff( "proton", 7, par_p );
		prot_det_eff_sbs8 = calc_det_eff( "proton", 8, par_p );
		prot_det_eff_sbs9 = calc_det_eff( "proton", 9, par_p );
		prot_det_eff_sbs14 = calc_det_eff( "proton", 14, par_p );

		neut_det_eff_sbs4 = calc_det_eff( "neutron", 4, par_n );
		neut_det_eff_sbs7 = calc_det_eff( "neutron", 7, par_n );
		neut_det_eff_sbs8 = calc_det_eff( "neutron", 8, par_n );
		neut_det_eff_sbs9 = calc_det_eff( "neutron", 9, par_n );
		neut_det_eff_sbs14 = calc_det_eff( "neutron", 14, par_n );

		TCanvas *c_pn_det_eff;
		if( Nucleon_Select == "both"){
			c_pn_det_eff = new TCanvas("c_pn_det_eff", "c_pn_det_eff", 1200, 500);
			h_p_det_eff->SetMarkerStyle(21);
			h_p_det_eff->SetMarkerSize(0.5);
			h_p_det_eff->SetMarkerColor(51);
			h_p_det_eff->GetYaxis()->SetRangeUser(det_eff_yaxis_min, 100);
			h_p_det_eff->GetYaxis()->SetTitle("HCal Detection Efficiency (%)");
			h_p_det_eff->GetXaxis()->SetTitle("Nucleon Momentum (GeV/c)");
			h_p_det_eff->SetStats(false);
			h_p_det_eff->Draw("P0");

			h_n_det_eff->SetMarkerStyle(20);
			h_n_det_eff->SetMarkerSize(0.5);
			h_n_det_eff->SetMarkerColor(3);
			h_n_det_eff->GetYaxis()->SetRangeUser(det_eff_yaxis_min, 100);
			h_n_det_eff->SetStats(false);
			h_n_det_eff->Draw("SAME+P0");

			h_n_det_eff->Draw("same+P0");

			tf_p_det_eff->Draw("SAME");
			tf_n_det_eff->Draw("SAME");

			TLegend *tleg_pn_det_eff = new TLegend(0.12, 0.77, 0.22, 0.87);
			tleg_pn_det_eff->AddEntry(h_p_det_eff, "Proton");
			tleg_pn_det_eff->AddEntry(h_n_det_eff, "Neutron");
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

			TPaveText *tpt_det_effs = new TPaveText(6.80, 82.5, 8.40, 92.5);
			tpt_det_effs->AddText("Detection Efficiencies");
			tpt_det_effs->AddText("");
			tpt_det_effs->AddText("SBS4:");
			tpt_det_effs->AddText(Form("proton: %0.1f%%, neutron: %0.1f%%", prot_det_eff_sbs4, neut_det_eff_sbs4));
			tpt_det_effs->AddText("");

			tpt_det_effs->AddText("SBS7:");
			tpt_det_effs->AddText(Form("proton: %0.1f%%, neutron: %0.1f%%", prot_det_eff_sbs7, neut_det_eff_sbs7));
			tpt_det_effs->AddText("");

			tpt_det_effs->AddText("SBS8:");
			tpt_det_effs->AddText(Form("proton: %0.1f%%, neutron: %0.1f%%", prot_det_eff_sbs8, neut_det_eff_sbs8));
			tpt_det_effs->AddText("");

			tpt_det_effs->AddText("SBS9:");
			tpt_det_effs->AddText(Form("proton: %0.1f%%, neutron: %0.1f%%", prot_det_eff_sbs9, neut_det_eff_sbs9));
			tpt_det_effs->AddText("");

			tpt_det_effs->AddText("SBS14:");
			tpt_det_effs->AddText(Form("proton: %0.1f%%, neutron: %0.1f%%", prot_det_eff_sbs14, neut_det_eff_sbs14));

			tpt_det_effs->Draw("same");

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

	}


	cout << "Writing to outfile...." << endl;
	outfile->Write();

	auto total_time_end = high_resolution_clock::now();
	auto total_time_duration = duration_cast<minutes>(total_time_end - total_time_start);
	cout << endl << "---------------------------------------------------" << endl;
	cout << "Total time for analysis: " << total_time_duration.count() << " minutes. " << endl;
	cout << endl << "Outfile: " << outfile_name.Data() << endl;
	// cout << endl << "Infile: " << infile_name.Data() << endl;

}