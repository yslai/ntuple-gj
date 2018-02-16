#include <cstdio>
#include <cmath>
#include <sys/time.h>

#include <vector>
#include <list>
#include <algorithm>

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <iostream>

//#define HI_TREE "hiEvtAnalyzer/HiTree"
#define HI_TREE "AliAnalysisTaskNTGJ/_tree_event"

namespace {

	typedef unsigned short index_t;

	size_t nevent(const char *filename)
	{
		TFile *root_file = TFile::Open(filename);

		if (root_file == NULL) {
			return 0;
		}

		TTree *hi_tree = dynamic_cast<TTree *>(root_file->Get(HI_TREE));

		if (hi_tree == NULL) {
			return 0;
		}

		const size_t ret = hi_tree->GetEntries();

		root_file->Close();

		return ret;
	}
}
  void range_extract(const char *filename, 
		     const int n_mult_bins, 
		     const int n_vert_bins, 
		     std::vector <float> &m_ranges, 
		     std::vector <float> &v_ranges) {
    
    TFile *root_file = TFile::Open(filename);
    
    if (root_file == NULL) {
    fprintf(stderr, "%s:%d:%s\n",__FILE__, __LINE__, "bin_extract: TFile Failed");
    return;
    }

    TTree *hi_tree = dynamic_cast<TTree *>(root_file->Get(HI_TREE));

    if (hi_tree == NULL) {
      fprintf(stderr, "%s:%d:%s\n",__FILE__, __LINE__, "bin_extract: TTree Failed");
      return; 
    }

    size_t nentries = hi_tree->GetEntries();

    //float centrality_v0m;

    double vtx[3];
    int multiplicity_size = 64; //64 channels to V0 detector
    float multiplicity_v0[multiplicity_size];

    hi_tree->SetBranchAddress("primary_vertex", vtx);
    hi_tree->SetBranchAddress("multiplicity_v0", multiplicity_v0);

    fprintf(stderr, "%s:%d:%s\n",__FILE__, __LINE__, "range_extract: Read TFile");

    std::vector <float> multps, vertxs;

    for (size_t i = 0; i < nentries; i++){
      hi_tree->GetEntry(i);
      
      float multp_sum = 0;
      for (int k = 0; k < multiplicity_size; k++) {
        multp_sum += multiplicity_v0[k];
      }
      
      fprintf(stderr, "%s:%d:%s %f %f\n",__FILE__, __LINE__, "range_extract:",vtx[2], multp_sum);
      multps.push_back(multp_sum);
      vertxs.push_back(vtx[2]);
    }
    
    //    std::sort(multps.begin(), multps.end(), std::greater<float>());
    std::sort(multps.begin(), multps.end());
    std::sort(vertxs.begin(), vertxs.end());
    
    fprintf(stderr, "%s:%d:%s\n",__FILE__, __LINE__, "range_extract: Sorted Vectors");
    
    for (int m = 0; m < (n_mult_bins+1); m++){
      float bin_start = multps[(m*nentries)/n_mult_bins];
      m_ranges.push_back(bin_start);
      fprintf(stderr, "%s:%d:%s: %f\n",__FILE__, __LINE__, "range_extract:",bin_start);
    }

    for(int v = 0; v<n_vert_bins+1; v++){
      float bin_start = vertxs[(v*nentries)/n_vert_bins];
      v_ranges.push_back(bin_start);
    }
  }


  std::vector<std::pair<int,int> >bin_extract(const char *filename,
					     std::vector <float> &m_ranges,
					     std::vector <float> &v_ranges)
  {
    TFile *root_file = TFile::Open(filename);
    
    if (root_file == NULL) {
      return std::vector<std::pair<int,int> >();
    }
    
    TTree *hi_tree = dynamic_cast<TTree *>(root_file->Get(HI_TREE));
    
    if (hi_tree == NULL) {
      return std::vector<std::pair<int,int> >();
    }
   
    size_t nentries = hi_tree->GetEntries(); 
    double vtx[3];
    // float centrality_v0m;
    int multiplicity_size = 64;
    float multiplicity_v0[multiplicity_size];

    hi_tree->SetBranchAddress("primary_vertex", vtx);
    hi_tree->SetBranchAddress("multiplicity_v0", multiplicity_v0);
    
    fprintf(stderr, "%s:%d:%s\n",__FILE__, __LINE__, "bin_extract: Read TFile");
    
    std::vector<std::pair<int,int> > ret;
    int mbin = 0;
    int vbin = 0;
    
    for (size_t i = 0; i < nentries; i++) {
      hi_tree->GetEntry(i);

      float multp_sum = 0;
      for (int k = 0; k < multiplicity_size; k++) {
	multp_sum += multiplicity_v0[k];
      }

      for (unsigned int m = 0; m<m_ranges.size()-1; m++){
	if (multp_sum >= m_ranges[m] && multp_sum <= m_ranges[m+1]){
	  mbin = m;
	  break;
	    }
      }
      for (unsigned int v = 0; v <v_ranges.size()-1; v++){
	if (vtx[2] >= v_ranges[v] && vtx[2] <= v_ranges[v+1]){
	  vbin = v;
	  break;
	}
      }
      
      fprintf(stderr, "%s:%d:%s:%lu %i %i\n",__FILE__, __LINE__, "bin_extract", i,mbin,vbin);
      
      ret.push_back(std::make_pair(mbin,vbin));
    }
    
    root_file->Close();
    
    return ret;
  }
  
void Write_TTree (const char *filename, TString ProximityBranch, TString LimitUseBranch, 
		  std::vector <std::vector<size_t> > Proximity_Matches,
		  std::vector <std::vector<size_t> > LimitUse_Matches){

    TFile *root_file = new TFile(filename,"update");
    TTree *hi_tree = dynamic_cast<TTree *>
      (dynamic_cast<TDirectoryFile *>
       (root_file->Get("AliAnalysisTaskNTGJ"))->Get("_tree_event"));

    TFile *newfile = new TFile("Bin_Mixed.root","recreate");
    TTree *newtree = hi_tree->CloneTree(0);
    
    unsigned int n_mix_events = 10;
    ULong64_t nentries = hi_tree->GetEntries();
    Long64_t Prox_Mix_Events[n_mix_events];
    Long64_t LimUse_Mix_Events[n_mix_events];

    fprintf(stderr, "%llu\n",nentries);
    
    //TBranch *MixE = newtree->Branch("Simple_Mix_Events", Simple_Mix_Events, "&Simple_Mix_Events[10]/L");
    //    TString Branch = Branchname + "[10]/L";
    TString ProximityBranchForm = ProximityBranch + "[10]/L";
    TBranch *PMixE = newtree->Branch(ProximityBranch, Prox_Mix_Events, ProximityBranchForm);
    TString LimitUseBranchForm = LimitUseBranch + "[10]/L";
    TBranch *LMixE = newtree->Branch(LimitUseBranch, LimUse_Mix_Events, LimitUseBranchForm);

    for (ULong64_t t = 0; t<nentries;t++){
      hi_tree->GetEntry(t);
      
      if(t < Proximity_Matches.size()){
	for (size_t s=1; s<(Proximity_Matches[t]).size();s++){
	  Prox_Mix_Events[s-1] = Proximity_Matches[t][s];
	  LimUse_Mix_Events[s-1] = LimitUse_Matches[t][s];
	  fprintf(stderr, "%llu:%lld %lld\n", t,Prox_Mix_Events[s-1],LimUse_Mix_Events[s-1]);
	}
      }
      
      else if (t >= Proximity_Matches.size()){
	for(size_t u = 0; u<n_mix_events; u++){
	  Prox_Mix_Events[u] = t; //Fill with own event number. Skip During correlation function                                                
	  LimUse_Mix_Events[u] = t;
	  fprintf(stderr, "%llu:%lld %lld\n", t,Prox_Mix_Events[u],LimUse_Mix_Events[u]);
	}
      }
      
      fprintf(stderr, "%s\n","");
      newtree->Fill();
      
    }//End loop over entries                                                                                                                                                                              
    newtree->Write();   
    delete root_file;
    delete newfile;
    
  }
 
  bool bin_compare(const std::pair<int, int> u,
		   const std::pair<int, int> v)
  {
    return u.first == v.first && u.second == v.second;
  }
  

  void Match_Events(const char *filename){
    
    int n_multp_bins = 20;
    int n_vert_bins = 7;

    std::vector <float> multp_ranges, vertx_ranges;
    range_extract(filename, n_multp_bins, n_vert_bins, multp_ranges, vertx_ranges) ;

    std::vector<std::pair<int,int> > Binned_Events = bin_extract(filename, multp_ranges, vertx_ranges);
    fprintf(stderr, "%s:%d:%s: %s\n", __FILE__, __LINE__,"Match_Events","Bins Extracted");

    size_t n_mix = 10;
    std::vector<std::vector<size_t> > Proximity_Matches, LimitUse_Matches;
    std::vector <size_t> flat;

    fprintf(stderr, "%s:%d:%s: %s\n", __FILE__, __LINE__,"Match_Events","Starting Matching");

    for (size_t i = 0; i < Binned_Events.size(); i++){
      fprintf(stderr, "%s:%d:%s: %lu\n", __FILE__, __LINE__,"Match_Events Proximity",i);
      flat.push_back(i);
      size_t j;
      //      if (i>500) j = i-500;//rough sliding window
      j = i+1;
      //      else j = 0;
      while (flat.size() < n_mix+1){
	if (j==i) {
	  j++;
	  continue;
	}
	//	if (j == Binned_Events.size()) j = i-1000;
	if (j == Binned_Events.size()) j = i-1500;
	if (bin_compare(Binned_Events[i],Binned_Events[j])){
	  flat.push_back(j);
	}
	j++;
      }

      Proximity_Matches.push_back(flat);
      flat.clear();
    }
    fprintf(stderr, "%s:%d:%s: %s\n", __FILE__, __LINE__,"Match_Events","Proximity Mixed Events Done");

    int n_events = Binned_Events.size();
    int counters [n_events];
    counters [n_events] = {0};
    std::cout<<counters[3]<<std::endl;

    for (size_t i = 0; i < Binned_Events.size(); i++){
    //for(size_t i = 0; i < 1; i++){
      fprintf(stderr, "%s:%d:%s: %lu\n", __FILE__, __LINE__,"Match_Events Limit Use",i);
      flat.push_back(i);
      size_t j=0;
      while (flat.size() < n_mix+1){
	if (j == Binned_Events.size()) j = 0;
	if (j==i){
	    j++;
	    continue;
	  }
	if (counters[j] >= 50){
	  j++;
	  continue;
	}
	if (bin_compare(Binned_Events[i],Binned_Events[j])){
	  flat.push_back(j);
	  counters[j] = counters[j]+1;
	  if(i%10000==0){
	    fprintf(stderr, "%s:%d:%s: %i\n", __FILE__, __LINE__,"Match_Events Limit Use",counters[j]); 
	  }
	}
	j++;
      }
      LimitUse_Matches.push_back(flat);
      flat.clear();
    }

    TString ProximityBranch = "Proximity_Mixed_Events";
    TString LimitUseBranch = "LimitUse_Mixed_Events";
    fprintf(stderr, "%s:%d:%s: %s\n", __FILE__, __LINE__,"Match_Events","Matches Filled. Writing to TTree");

    Write_TTree(filename,ProximityBranch,LimitUseBranch,Proximity_Matches,LimitUse_Matches);
    fprintf(stderr, "%s:%d:%s: %s\n", __FILE__, __LINE__,"Match_Events","Done");

    Proximity_Matches.clear();
    LimitUse_Matches.clear();
  }

int main(int argc, char *argv[])
{
	if (argc < 2) {
	  fprintf(stderr,"%s\n","Argument Syntax is [Command] [File]");
		return EXIT_FAILURE;
	}
	Match_Events(argv[1]);
}
