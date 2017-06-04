// Low memory complexity ntuple merger. Adapted from CMS Heavy Ion's
// ntuple merger (https://github.com/richard-cms/hiForestMerging),
// with performance optimization and compatibility with ROOT 5.x/CINT.

#include <TSystem.h>
#include <TChain.h>
#include <TFile.h>

TString fix_tchain_glob(const char *glob, const char *output_filename)
{
    // Note: Only one level of subdirectory will get cleaned up

    const TString temp_base = "merge_temp/merge_temp_" +
        TString(gSystem->BaseName(output_filename)).
        ReplaceAll(".root", "");

    if (temp_base.Index("/") != -1) {
        if (gSystem->AccessPathName(gSystem->DirName(temp_base))) {
            gSystem->Unlink(gSystem->DirName(temp_base));
        }
        gSystem->mkdir(gSystem->DirName(temp_base), kTRUE);
    }
    gSystem->Exec(TString("i=0; for f in ") + glob +
                  "; do ln -sf \"$(readlink -f \"$f\")\" " +
                  temp_base +
                  "_\"$(printf %03d \"$i\")\".root; i=$(($i + 1));" +
                  " done");

    return temp_base;
}

void clean_up(TString temp_base)
{
    gSystem->Exec(TString("for f in ") + temp_base +
                  "_*.root; do [ -L \"$f\" ] && rm -f \"$f\"; done");
    if (temp_base.Index("/") != -1) {
        gSystem->Unlink(gSystem->DirName(temp_base));
    }
}

void search_ntuple(TObjArray &dir_name, TObjArray &dir_tree_name,
                   const char *glob, Bool_t extended_glob,
                   TString temp_base)
{
    TChain *test_chain =
        new TChain("AliAnalysisTaskNTGJ/_tree_event");

    fprintf(stderr, "Search for ntuples... ");
    test_chain->Add(glob, 1);
    if (!(test_chain->GetEntries() > 0)) {
        fprintf(stderr, "error: no entries found, abort.\n");
        if (extended_glob) {
            clean_up(temp_base);
        }
        gSystem->Exit(1);
    }

    TFile *test_file = test_chain->GetFile();
    TList *key_list_1 = test_file->GetListOfKeys();

    for (Int_t i = 0; i < key_list_1->GetEntries(); i++) {
        TDirectoryFile *dir_file = (TDirectoryFile *)
            test_file->Get(key_list_1->At(i)->GetName());
        if (strcmp(dir_file->ClassName(), "TDirectoryFile") != 0) {
            continue;
        }
        TList *key_list_2 = dir_file->GetListOfKeys();
        for (Int_t j = 0; j < key_list_2->GetEntries(); j++) {
            TString n = dir_file->GetName();

            n += "/";
            n += key_list_2->At(j)->GetName();

            const TTree *t = (TTree *)test_file->Get(n);

            if (t != NULL &&
                (strcmp(t->ClassName(), "TTree") == 0 ||
                 strcmp(t->ClassName(), "TNtuple") == 0) &&
                ((dir_tree_name.GetEntries() == 0) ||
                 // skip duplicate tree entries
                 n != *(TString *)dir_tree_name.Last())) {
                dir_name.AddLast((TObject *)
                                 (new TString(dir_file->GetName())));
                dir_tree_name.AddLast((TObject *)(new TString(n)));
            }
        }
    }

    fprintf(stderr, "done.\n");

    delete test_chain;
}

void merge_ntuple(const char *output_filename, TObjArray &dir_name,
                  TObjArray &dir_tree_name, const char *glob)
{
    TObjArray chain;

    for (Int_t i = 0; i < dir_tree_name.GetEntries(); i++) {
        TString *dtn = (TString *)dir_tree_name.At(i);

        chain.AddLast((TObject *)(new TChain(dtn->Data())));
        ((TChain *)chain.Last())->Add(glob);
        fprintf(stderr, "Created merged ntuple %d: %s\n", i,
                dtn->Data());
    }

    TFile output_file(output_filename, "RECREATE");

    for (Int_t i = 0; i < dir_tree_name.GetEntries(); i++) {
        output_file.cd();

        TString *dn = (TString *)dir_name.At(i);
        TString *dtn = (TString *)dir_tree_name.At(i);

        fprintf(stderr, "Processing %s... ", dtn->Data());

        if (i == 0 || *dtn != *(TString *)dir_tree_name.At(i - 1)) {
            (output_file.mkdir(*dn))->cd();
        }
        else {
            output_file.cd(*dn);
        }
        ((TChain *)chain[i])->Merge(&output_file, 0, "keep");

        fprintf(stderr, "done.\n");

        delete (TChain *)chain[i];
        delete dn;
        delete dtn;
    }
    output_file.Close();

    fprintf(stderr, "Output %s produced.\n", output_filename);
}

void MergeNtuple(
    const char *glob = "save/lhc16c2-180720-1/*/AnalysisResults.root",
    const char *output_filename = "lhc16c2-180720-1.root",
    const bool extended_glob = true)
{
    TString temp_base;
    TString *new_glob = NULL;

    if (extended_glob) {
        temp_base =
            fix_tchain_glob(glob, output_filename);

        TString *new_glob = new TString(temp_base + "_*.root");

        glob = *new_glob;
        fprintf(stderr, "New glob: %s\n", glob);
    }

    TObjArray dir_name;
    TObjArray dir_tree_name;

    search_ntuple(dir_name, dir_tree_name, glob, extended_glob,
                  temp_base);
    merge_ntuple(output_filename, dir_name, dir_tree_name, glob);

    if (extended_glob) {
        clean_up(temp_base);
        if (new_glob != NULL) {
            delete new_glob;
        }
    }
    gSystem->Exit(0);
}
