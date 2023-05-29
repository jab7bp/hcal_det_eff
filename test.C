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

TString run_target = "LH2";
int kine = 4;
int sbsfieldscale = 0;

vector<TString> filenames;
TString rootfile_dir = "/work/halla/sbs/sbs-gmn/pass0/SBS4/LH2/rootfiles";
TString infile = "e1209019_fullreplay_11616_stream0_seg";
TString infile_name;
vector<int> runnum_vec;
TFile *temp_file;
TH1D *hin_bb_gem_Ntracks;

void test(){

	hin_bb_gem_Ntracks = new TH1D("hin_bb_gem_Ntracks", "tracks", 11, -0.5, 10.5);
	for( int run = 0; run < 6; run++ ){
		runnum_vec.push_back(lookup_parsed_runnums( run_target, kine, sbsfieldscale, run ));
		infile_name = Form("e1209019_fullreplay_%i_stream0_seg*", runnum_vec[run]);
		cout << runnum_vec[run] << "  ";
		list_files(rootfile_dir, filenames, infile_name.ReplaceAll("*", ""));
	}
	for( int file = 0; file < filenames.size(); file++ ){
		cout << "File " << file << ", filename: " << filenames[file].Data() << endl;
		temp_file = new TFile(filenames[file].Data(), "READ");


				hin_bb_gem_Ntracks->Add( static_cast<TH1D*>(temp_file->Get("hbb_gem_Ntracks")) );	
		
	}	
	
}
