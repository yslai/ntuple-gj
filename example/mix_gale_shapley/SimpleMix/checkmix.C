#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <iostream>
#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>

#define HI_TREE "AliAnalysisTaskNTGJ/_tree_event"


int checkmix(){
  TFile *root_file = new TFile("mixed_simple.root","update");
  // TTree *hi_tree = dynamic_cast<TTree *>
  // (dynamic_cast<TDirectoryFile *>
  // (root_file->Get("AliAnalysisTaskNTGJ"))->Get("_tree_event"));

  TTree *hi_tree = dynamic_cast<TTree *>(root_file->Get("_tree_event"));

 Long64_t Simple_Mix_Events[10];
 Long64_t run_number;
 hi_tree->SetBranchAddress("Simple_Mix_Events",Simple_Mix_Events);

 Long64_t nentries = hi_tree->GetEntries();
 std::cout<<nentries<<std::endl;
 for (Long64_t n = 0; n<nentries;n++){
   if (n%10000 == 0){
     hi_tree->GetEntry(n);
     for(int i = 0; i< 10;i++){
       std::cout<<n<<": "<<Simple_Mix_Events[i]<<std::endl;
     }
     std::cout<<std::endl;
   }
 }

 return 0;
}
