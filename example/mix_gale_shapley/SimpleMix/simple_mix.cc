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
      
      fprintf(stderr, "%s:%d:%f %f\n",__FILE__, __LINE__, vtx[2], multp_sum);
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
  
  void Write_TTree (const char *filename, std::vector <std::vector<size_t> > Matches){

    TFile *root_file = new TFile(filename,"update");
    TTree *hi_tree = dynamic_cast<TTree *>
      (dynamic_cast<TDirectoryFile *>
       (root_file->Get("AliAnalysisTaskNTGJ"))->Get("_tree_event"));

    TFile *newfile = new TFile("mixed_simple.root","recreate");
    TTree *newtree = hi_tree->CloneTree(0);
    
    unsigned int n_mix_events = 10;
    ULong64_t nentries = hi_tree->GetEntries();
    Long64_t Simple_Mix_Events[n_mix_events];
    
    fprintf(stderr, "%llu\n",nentries);
    
    TBranch *MixE = newtree->Branch("Simple_Mix_Events", Simple_Mix_Events, "&Simple_Mix_Events[10]/L");
    
    for (ULong64_t t = 0; t<nentries;t++){
      hi_tree->GetEntry(t);
      
      if(t < Matches.size()){
	for (size_t s=1; s<(Matches[t]).size();s++){
	  Simple_Mix_Events[s-1]=Matches[t][s];
	  fprintf(stderr, "%llu:%lld\n", t,Simple_Mix_Events[s-1]);
	}
      }
      
      else if (t >= Matches.size()){
	for(size_t u = 0; u<n_mix_events; u++){
	  Simple_Mix_Events[u] = t; //Fill with own event number. Skip During correlation function                                                                                                               
	  fprintf(stderr, "%llu:%lld\n",t,Simple_Mix_Events[u]);
	}
      }
      
      fprintf(stderr, "%s\n","");

      newtree->Fill();
      
    }//End loop over entries                                                                                                                                                                              
    newtree->Write();   
    delete root_file;
    delete newfile;
    
    
    gSystem->Exit(0);
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
    std::vector<std::vector<size_t> > Matches;
    std::vector <size_t> flat;

    for (size_t i = 0; i < Binned_Events.size(); i++){
      flat.push_back(i);
      size_t j;
      if (i>500) j = i-500;//rough sliding window
      else j = 0;
      while (flat.size() < n_mix+1){
	if (j==i) continue;
	if (j == Binned_Events.size()) j = 0;
	if (bin_compare(Binned_Events[i],Binned_Events[j])){
	    flat.push_back(j);
	  }
	j++;
      }
      Matches.push_back(flat);
      flat.clear();
    }
    Write_TTree(filename, Matches);    
  }




 //  void feature_normalize(std::vector<float> &u,
// 						   std::vector<float> &v, const size_t n)
// 	{
// 		std::vector<float> s(n, 0);

// 		for (size_t j = 0; j < n; j++) {
// 			float s_j = 0;
// #ifdef _OPENMP
// #pragma omp parallel for shared(u) reduction(+: s_j)
// #endif // _OPENMP
// 			for (size_t i = 0; i < u.size(); i += n) {
// 				s_j += fabsf(u[i + j]);
// 			}
// 			s[j] = s_j;
// 		}
// 		for (size_t j = 0; j < n; j++) {
// 			float s_j = 0;

// #ifdef _OPENMP
// #pragma omp parallel for shared(u) reduction(+: s_j)
// #endif // _OPENMP
// 			for (size_t i = 0; i < v.size(); i += n) {
// 				s_j += fabsf(v[i + j]);
// 			}
// 			s[j] += s_j;
// 		}
// 		for (size_t j = 0; j < n; j++) {
// 			s[j] = (u.size() + v.size()) / s[j];

// 			fprintf(stderr, "%s:%d: %lu %f\n", __FILE__, __LINE__, j, s[j]);
// 		}

// #ifdef _OPENMP
// #pragma omp parallel for shared(u, s)
// #endif // _OPENMP
// 		for (size_t i = 0; i < u.size(); i += n) {
// 			for (size_t j = 0; j < n; j++) {
// 				u[i + j] *= s[j];
// 			}
// 		}
// #ifdef _OPENMP
// #pragma omp parallel for shared(v, s)
// #endif // _OPENMP
// 		for (size_t i = 0; i < v.size(); i += n) {
// 			for (size_t j = 0; j < n; j++) {
// 				v[i + j] *= s[j];
// 			}
// 		}
// 	}
// }



// void mix_gale_shapley(const char *filename_0, const char *filename_1,
// 					  const int nfeature, const int nduplicate)
// {
// 	const size_t nevent_0 = nevent(filename_0);
// 	const size_t nevent_1 = nevent(filename_1);

// 	const size_t block_size_max = 2000;

// 	const size_t nblocks = std::min(nevent_0, nevent_1 * nduplicate) /
// 		block_size_max + 1;
// 	const size_t nblock = nblocks - 1; // FIXME:Use % and rounding to get all events 

// 	size_t lmin,lmax; 
// 	size_t width = 5; //if changed, also must change when writing to Tree

// 	std::vector<std::vector<Long64_t> > Matches;

// 	for(size_t h = 0; h < nblock+1; h++){
// 	  //for(size_t h = 0; h < 1; h++){
// 	    const size_t event_start_0 = h * nevent_0 / (nblock + 1); 
// 	    const size_t event_end_0 = (h + 1) * nevent_0 / (nblock + 1);
// 	    const size_t nevents_0 = event_end_0 - event_start_0;	  
	    
// 	    std::vector<std::vector<Long64_t> > k;		  
	  
// 	    std::vector<float> feature_0 =
// 	      feature_extract(filename_0, event_start_0, event_end_0,
// 			      nfeature);
	    
// 	    fprintf(stderr,"%s %lu %s %lu\n","Block",h,"of",nblock);
	    
// 	    if (h < width) {
// 	      lmin = 0; 
// 	      lmax = 2*width+1;
// 	    }

// 	    else if (h+width > nblock) {
// 	    lmin = nblock-2*width;
// 	    lmax = nblock+1;
// 	    }

// 	    else {
// 	      lmin = h-width;  	 
// 	      lmax = h+width+1;
// 	    }
	    
// 	    for (size_t i = lmin; i < lmax; i++) {
	      
// 	      if (i==h) continue;
	     	      
// 	      size_t event_start_1 = i * nevent_1 / (nblock+1);
// 	      size_t event_end_1 = (i + nduplicate) * nevent_1 /
// 		(nblock + nduplicate);
	      
// 	      const size_t nevents_1 = event_end_1 - event_start_1;
// 	      if(nevents_1<nevents_0) event_end_1 += (nevents_0-nevents_1);
// 	      //FIXME:small # events mix with themselves once. need conditional using last block
	      
// 		std::vector<float> feature_1 =
// 			feature_extract(filename_1, event_start_1, event_end_1,
// 							nfeature);

// 		{
// 			std::vector<float> feature_0_scaled = feature_0;
// 			std::vector<float> feature_1_scaled = feature_1;
			
// 			feature_normalize(feature_0_scaled, feature_1_scaled,
// 							  nfeature);

// 			std::vector<std::list<index_t> > preference_0;
// 			std::vector<std::list<index_t> > preference_1;

// 			order_preference(preference_0, preference_1,
// 							 feature_0_scaled, feature_1_scaled,
// 							 nfeature, nduplicate);

// 			std::vector <index_t> m;
// 			m = gale_shapley(preference_0, preference_1);

// 			const size_t feature_1_size_nfeature =
// 			  feature_1.size() / nfeature;

// 			std::vector <Long64_t>tempflat;

// 			for (size_t j = 0; j < m.size(); j++) {                                        
// 			  Long64_t q = event_start_1 + (m[j] % feature_1_size_nfeature);
// 			  tempflat.push_back(q);
// 			}
			
// 			k.push_back(tempflat);
// 		}
// 	  } //i
// 	  for (size_t j = 0; j < k[0].size(); j++){
// 	    std::vector <Long64_t> P;
// 	    P.push_back(event_start_0+j);
// 	    for (size_t l = 0; l < k.size(); l++){
// 	      P.push_back(k[l][j]);
// 	    }
// 	    Matches.push_back(P);
// 	  }

// 	}//h

// 	  // //	  write to txt
// 	  // FILE * txtfile = fopen ("pairs.txt","w");
// 	  // for (size_t t=0; t<Matches.size();t++){
// 	  //   for (size_t s=1; s<Matches[t].size();s++){
// 	  //     fprintf(txtfile, "%lld ", Matches[t][s]);
// 	  //   }
// 	  //   fprintf(txtfile, "%s\n","");
// 	  // }
// 	  // fclose (txtfile);
	  
// 	  //	  write to TTree
// 	if (strcmp(filename_0,filename_1) == 0){

// 	  TFile *root_file = new TFile(filename_0,"update");
// 	  TTree *hi_tree = dynamic_cast<TTree *>
// 	    (dynamic_cast<TDirectoryFile *>
// 	     (root_file->Get("AliAnalysisTaskNTGJ"))->Get("_tree_event"));

// 	  TFile *newfile = new TFile("mixed.root","recreate");	  
// 	  TTree *newtree = hi_tree->CloneTree(0);

// 	  unsigned int n_mix_events = 2*width;
// 	  ULong64_t nentries = hi_tree->GetEntries();    
// 	  Long64_t Mix_Events[n_mix_events];

// 	  fprintf(stderr, "%llu\n",nentries);
	  
// 	  TBranch *MixE = newtree->Branch("Mix_Events", Mix_Events, "&Mix_Events[10]/L");
	  
// 	  for (ULong64_t t = 0; t<nentries;t++){
// 	    hi_tree->GetEntry(t);
	    
// 	    if(t < Matches.size()){
// 	      for (size_t s=1; s<(Matches[t]).size();s++){
// 		Mix_Events[s-1]=Matches[t][s]; 
// 		fprintf(stderr, "%llu:%lld\n", t,Mix_Events[s-1]);
// 	      }
// 	    }
	    
// 	    else if (t >= Matches.size()){
// 	      for(size_t u = 0; u<n_mix_events; u++){
// 		Mix_Events[u] = t; //Fill with own event number. Skip During correlation function
// 		fprintf(stderr, "%llu:%lld\n",t,Mix_Events[u]);
// 	      }
// 	    }
	    
// 	    fprintf(stderr, "%s\n","");
// 	    //MixE->Fill();
// 	    newtree->Fill();  
	    
// 	  }//End loop over entries
// 	  //newtree->AutoSave();    
// 	  newtree->Write();
// 	  //newfile->ls;

// 	  delete root_file;
// 	  delete newfile;
// 	}
// 	else fprintf(stderr, "%s\n","Nothing written to root file.");
	
// 	gSystem->Exit(0);
// }

int main(int argc, char *argv[])
{
	if (argc < 2) {
	  fprintf(stderr,"%s\n","Argument Syntax is [Command] [File]");
		return EXIT_FAILURE;
	}
	Match_Events(argv[1]);
}
