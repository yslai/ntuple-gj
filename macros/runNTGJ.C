#if 0
set +o posix; function join_by { local d=$1; shift; echo -n "$1";
    shift; printf "%s" "${@/#/$d}"; }
root=root; exec $root -l -b -q "$0($(join_by \",\" \"$*\" | \
/bin/sed s/\"\\\([0-9]\\+\\\)\"/\\1/g\;s/^\"\"\$//))"; exit 0
#endif

#include <TROOT.h>
#include <TSystem.h>

void runNTGJ(const char *config_filename = "config/lhc16c2_1run.yaml",
             const char *run_mode = "test")
{
    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
    gROOT->ProcessLine(".include $ALICE_PHYSICS/include");

    FILE *fp;
    char line[4096];

    // Light-weight parsing of a YAML conform file, with what is
    // available in ROOT/CINT, pass 1
    fp = fopen(config_filename, "r");
    // Default values
    TList package_list;

    while (fgets(line, 4096, fp) != NULL) {
        // Skip comments
        if (line[0] == '#') {
            continue;
        }

        char key[4096];
        char dummy[4096];
        char value[4096];

        key[0] = '\0';
        value[0] = '\0';
        // Key is arbitrarily many characters until ':', followed by
        // arbitrarily many combinations of space/tabs (both assuming
        // user will not actually enter >= 4096 of those), then at
        // most 4096 caracters of value (which could be very long, for
        // run lists). Tailing spaces for numerical values are
        // generally passed on to CINT, which also does not become an
        // issue.
        sscanf(line, "%[^:]:%[ \t]%4096[^\n\r]", key, dummy, value);

        if (strcmp(key, "package") == 0) {
            package_list.Add((TObject *)(new TObjString(value)));
        }
    }
    package_list.SetOwner();
    fclose(fp);

    if (package_list.GetSize() == 0) {
        package_list.Add((TObject *)(new TObjString(
            "fastjet::v3.2.1_1.024-alice1-4")));
        // Compiling against CGAL requires explicit include paths to
        // Boost, MPFR, and GMP
        package_list.Add((TObject *)(new TObjString(
            "cgal::v4.6.3-18")));
        package_list.Add((TObject *)(new TObjString(
            "boost::v1.59.0-14")));
        package_list.Add((TObject *)(new TObjString(
            "MPFR::v3.1.3-4")));
        package_list.Add((TObject *)(new TObjString(
            "GMP::v6.0.0-2")));
        // The grid ROOT packages tend to lack a GCC 4.9.x =
        // CXXABI_1.3.8 dependency (see also
        // https://gcc.gnu.org/onlinedocs/libstdc++/manual/abi.html),
        // which causes jobs to fail at run time with the message:
        //
        //   root: /usr/lib64/libstdc++.so.6: version `CXXABI_1.3.8'
        //   not found (required by root)
        //   Output file AnalysisResults.root not found. Job FAILED !
        package_list.Add((TObject *)(new TObjString(
            "GCC-Toolchain::v4.9.3-alice3-1")));
    }

    for (Int_t i = 0; i < package_list.GetSize(); i++) {
        const TString &ps =
            ((TObjString *)(package_list.At(i)))->String();

        if (strncmp(ps.Data(), "GCC-Toolchain::", 15) != 0) {
            TString include = "$ALICE_ROOT/../../" +
                ps.Copy().ReplaceAll("::", "/") + "/include";

            if (gSystem->AccessPathName(include.ReplaceAll(
                "$ALICE_ROOT", gSystem->Getenv("ALICE_ROOT")))) {
                // Assuming a non-CVMFS installation, where the
                // symbolic "latest" should exist (note
                // TString::ReplaceAll already modified the content)
                include = "$ALICE_ROOT/../../" +
                    ps(0, ps.Index("/")) + "/latest/include";
            }
            gROOT->ProcessLine(".include " + include);
        }
    }

    // Load base root libraries
    gSystem->Load("libTree");
    gSystem->Load("libGeom");
    gSystem->Load("libVMC");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libPhysics");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");

    // Load analysis framework libraries
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");

    gSystem->Load("libEMCALUtils");
    gSystem->Load("libPWGPPEMCAL");

    if (gSystem->AccessPathName("libgmp.so") &&
        !gSystem->AccessPathName("/usr/lib64/libgmp.so")) {
        gSystem->Symlink("/usr/lib64/libgmp.so", "libgmp.so");
    }
    gSystem->Load("libgmp");
    if (gSystem->AccessPathName("libmpfr.so") &&
        !gSystem->AccessPathName("/usr/lib64/libmpfr.so")) {
        gSystem->Symlink("/usr/lib64/libmpfr.so", "libmpfr.so");
    }
    gSystem->Load("libmpfr");
    gSystem->Load("libCGAL");
    gSystem->Load("libfastjet");
    gSystem->Load("libsiscone");
    gSystem->Load("libsiscone_spherical");
    gSystem->Load("libfastjetplugins");
    gSystem->Load("libfastjetcontribfragile");

    gROOT->ProcessLine(".L AliAnalysisTaskNTGJ.cxx+g");

    AliAnalysisManager *mgr = new AliAnalysisManager();

    gROOT->ProcessLine(".x $ALICE_ROOT/ANALYSIS/macros/train/"
                       "AddESDHandler.C");
    gROOT->ProcessLine(".x $ALICE_ROOT/ANALYSIS/macros/train/"
                       "AddMCHandler.C");

    AliAnalysisAlien *plugin = new AliAnalysisAlien("pluginNTGJ");

    for (Int_t i = 0; i < package_list.GetSize(); i++) {
        const TString &ps =
            ((TObjString *)(package_list.At(i)))->String();

        fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, ps.Data());
        plugin->AddExternalPackage(ps.Data());
    }

    TString include_path("-I. -I$ALICE_ROOT/include "
                         "-I$ALICE_PHYSICS/include ");

    for (Int_t i = 0; i < package_list.GetSize(); i++) {
        const TString &ps =
            ((TObjString *)(package_list.At(i)))->String();

        if (strncmp(ps.Data(), "GCC-Toolchain::", 15) != 0) {
            include_path += "-I$ALICE_ROOT/../../" +
                ps.Copy().ReplaceAll("::", "/") + "/include ";
        }
    }
    package_list.Delete();

    plugin->AddIncludePath(include_path);

    // Light-weight parsing of a YAML conform file, with what is
    // available in ROOT/CINT, pass 2
    fp = fopen(config_filename, "r");

    // Default values
    TString emcal_correction_filename = "emcal_correction.yaml";
    TString emcal_geometry_filename = "";
    TString emcal_local2master_filename = "";
    bool mult_selection = true;
    bool physics_selection = false;
    bool physics_selection_mc_analysis = false;
    bool physics_selection_pileup_cut = true;
    bool force_ue_subtraction = false;
    // 1e+309 = INFINITY (for IEEE 754 double)
    TString skim_cluster_min_e = "-1e+309";
    TString skim_track_min_pt = "-1e+309";
    TString skim_muon_track_min_pt = "-1e+309";
    TString skim_jet_min_pt_1 = "-1e+309";
    TString skim_jet_min_pt_2 = "-1e+309";
    TString skim_jet_min_pt_3 = "-1e+309";
    TString skim_jet_average_pt = "-1e+309";
    // -2147483648 = INT_MIN
    TString skim_multiplicity_tracklet_min_n = "-2147483648";
    TString stored_track_min_pt = "-1e+309";
    TString stored_jet_min_pt_raw = "-1e+309";
    TString nrandom_isolation = "0";

    while (fgets(line, 4096, fp) != NULL) {
        if (line[0] == '#') {
            continue;
        }

        char key[4096];
        char dummy[4096];
        char value[4096];

        key[0] = '\0';
        value[0] = '\0';
        // See for pass 1 for explanation
        sscanf(line, "%[^:]:%[ \t]%4096[^\n\r]", key, dummy, value);

        // AliEn related options

        if (strcmp(key, "gridWorkingDir") == 0) {
            plugin->SetGridWorkingDir(value);
        }
        else if (strcmp(key, "gridOutputDir") == 0) {
            plugin->SetGridOutputDir(value);
        }
        else if (strcmp(key, "aliROOTVersion") == 0) {
            plugin->SetAliROOTVersion(value);
        }
        else if (strcmp(key, "aliPhysicsVersion") == 0) {
            plugin->SetAliPhysicsVersion(value);
        }
        else if (strcmp(key, "gridDataDir") == 0) {
            plugin->SetGridDataDir(value);
        }
        else if (strcmp(key, "dataPattern") == 0) {
            plugin->SetDataPattern(value);
        }
        else if (strcmp(key, "inputFormat") == 0) {
            plugin->SetInputFormat(value);
        }
        else if (strcmp(key, "masterResubmitThreshold") == 0) {
            plugin->SetMasterResubmitThreshold(atoi(value));
        }
        else if (strcmp(key, "maxInitFailed") == 0) {
            plugin->SetMaxInitFailed(atoi(value));
        }
        else if (strcmp(key, "nrunsPerMaster") == 0) {
            plugin->SetNrunsPerMaster(atoi(value));
        }
        else if (strcmp(key, "numberOfReplicas") == 0) {
            plugin->SetNumberOfReplicas(atoi(value));
        }
        else if (strcmp(key, "overwriteMode") == 0) {
            plugin->SetOverwriteMode(atoi(value));
        }
        else if (strcmp(key, "runNumber") == 0) {
            // Split an array of run numbers
            for (const char *v = value; *v != '\0';) {
                // Move forward until *v is a digit, but not further
                // than the end of the C string
                while (*v != '\0' && !isdigit(*v)) {
                    v++;
                }

                // Convert the digit into an integer and add to AliEn
                int n;

                sscanf(v, "%d", &n);
                plugin->AddRunNumber(n);
                // Move forward until *v is a non-digit, i.e. skipping
                // the number characters processed by sscanf(). Note
                // isdigit('\0') = 0 = false.
                while (isdigit(*v)) {
                    v++;
                }
            }
        }
        else if (strcmp(key, "runPrefix") == 0) {
            plugin->SetRunPrefix(value);
        }
        else if (strcmp(key, "splitMaxInputFileNumber") == 0) {
            plugin->SetSplitMaxInputFileNumber(atoi(value));
        }
        else if (strcmp(key, "splitMode") == 0) {
            plugin->SetSplitMode(value);
        }
        else if (strcmp(key, "ttl") == 0) {
            plugin->SetTTL(atoi(value));
        }

        else if (strcmp(key, "package") == 0) {
            // Do nothing (see pass 1)
        }

        // Task related options

        else if (strcmp(key, "emcalCorrection") == 0) {
            emcal_correction_filename = value;
        }
        else if (strcmp(key, "multSelection") == 0) {
            mult_selection = strncmp(value, "true", 4) == 0;
        }
        else if (strcmp(key, "physicsSelection") == 0) {
            physics_selection = strncmp(value, "true", 4) == 0;
        }
        else if (strcmp(key, "physicsSelectionMCAnalysis") == 0) {
            physics_selection_mc_analysis =
                strncmp(value, "true", 4) == 0;
        }
        else if (strcmp(key, "physicsSelectionPileupCut") == 0) {
            physics_selection_pileup_cut =
                strncmp(value, "true", 4) == 0;
        }
        else if (strcmp(key, "emcalGeometryFilename") == 0) {
            emcal_geometry_filename = value;
        }
        else if (strcmp(key, "emcalLocal2MasterFilename") == 0) {
            emcal_local2master_filename = value;
        }
        else if (strcmp(key, "forceUESubtraction") == 0) {
            force_ue_subtraction = strncmp(value, "true", 4) == 0;
        }
        else if (strcmp(key, "skimClusterMinE") == 0) {
            skim_cluster_min_e = value;
        }
        else if (strcmp(key, "skimTrackMinPt") == 0) {
            skim_track_min_pt = value;
        }
        else if (strcmp(key, "skimMuonTrackMinPt") == 0) {
            skim_muon_track_min_pt = value;
        }
        else if (strcmp(key, "skimJetMinPt1") == 0) {
            skim_jet_min_pt_1 = value;
        }
        else if (strcmp(key, "skimJetMinPt2") == 0) {
            skim_jet_min_pt_2 = value;
        }
        else if (strcmp(key, "skimJetMinPt3") == 0) {
            skim_jet_min_pt_3 = value;
        }
        else if (strcmp(key, "skimJetAveragePt") == 0) {
            skim_jet_average_pt = value;
        }
        else if (strcmp(key, "skimMultiplicityTrackletMinN") == 0) {
            skim_multiplicity_tracklet_min_n = value;
        }
        else if (strcmp(key, "storedTrackMinPt") == 0) {
            stored_track_min_pt = value;
        }
        else if (strcmp(key, "storedJetMinPtRaw") == 0) {
            stored_jet_min_pt_raw = value;
        }
        else if (strcmp(key, "nrandomIsolation") == 0) {
            nrandom_isolation = value;
        }
    }
    fclose(fp);

    // Intel MKL

    TString mkl_filename = "";
    TString oadb_filename = "";

    if (!gSystem->AccessPathName("libmkl_core_so")) {
        mkl_filename = "libiomp5_so libmkl_avx2_so libmkl_avx_so libmkl_core_so libmkl_intel_lp64_so libmkl_intel_thread_so";
    }
    if (strchr(emcal_geometry_filename.Data(), '/') == NULL) {
        oadb_filename += emcal_geometry_filename;
        oadb_filename += " ";
    }
    if (strchr(emcal_local2master_filename.Data(), '/') == NULL) {
        oadb_filename += emcal_local2master_filename;
        oadb_filename += " ";
    }
    plugin->SetAdditionalLibs(
        "AliAnalysisTaskNTGJ.h "
        "AliAnalysisTaskNTGJ.cxx "
        "special_function.h mc_truth.h "
        "emcal_cell.h emcal.h isolation.h jet.h "
        "bad_channel.h "
        "eLut.cpp eLut.h half.cpp halfExport.h halfFunction.h "
        "half.h halfLimits.h toFloat.h "
        "keras_model.h keras_model.cc "
        // "blasdrv.h efp7.cc einstein_sum.h "
        "photon_discr.model "
        // Not sure if this helps against the missing pyqpar_ when
        // dlopen() "libAliPythia6.so"
        "libpythia6.so libAliPythia6.so "
        "libgmp.so libmpfr.so "
        "libCGAL.so libfastjet.so libsiscone.so "
        "libsiscone_spherical.so libfastjetplugins.so "
        "libfastjetcontribfragile.so " +
        mkl_filename + " " +
        emcal_correction_filename + " " +
        oadb_filename);
    plugin->SetAnalysisSource("AliAnalysisTaskNTGJ.cxx");

    // Honor alien_CLOSE_SE for the output also, e.g. when
    // alien_CLOSE_SE=ALICE::LBL::EOS is set for NERSC.

    const char *alien_close_se = gSystem->Getenv("alien_CLOSE_SE");

    if (alien_close_se != NULL) {
        const char *file = mgr->GetCommonFileName();

        plugin->SetDefaultOutputs(kFALSE);
        plugin->SetOutputFiles(Form(
            "%s@%s", file, alien_close_se));
        plugin->SetOutputArchive(Form(
            "log_archive.zip:stdout,stderr@%s "
            "root_archive.zip:%s,*.stat@%s",
            alien_close_se, file, alien_close_se));
    }

    plugin->SetRunMode(run_mode);
    mgr->SetGridHandler(plugin);

    TString add_task_line = ".x macros/AddAliAnalysisTaskNTGJ.C("
        "\"AliAnalysisTaskNTGJ\",\"";

    add_task_line += emcal_correction_filename + "\"," +
        (mult_selection ? "true" : "false") + "," +
        (physics_selection ? "true" : "false") + "," +
        (physics_selection_mc_analysis ? "true" : "false") + "," +
        (physics_selection_pileup_cut ? "true" : "false") + ",\"" +
        emcal_geometry_filename + "\",\"" +
        emcal_local2master_filename + "\"," +
        (force_ue_subtraction ? "true" : "false") + "," +
        skim_cluster_min_e + "," +
        skim_track_min_pt + "," +
        skim_muon_track_min_pt + "," +
        skim_jet_min_pt_1 + "," + skim_jet_min_pt_2 + "," +
        skim_jet_min_pt_3 + "," +
        skim_jet_average_pt + "," +
        skim_multiplicity_tracklet_min_n + "," +
        stored_track_min_pt + "," +
        stored_jet_min_pt_raw + "," +
        nrandom_isolation + ")";
    fprintf(stderr, "%s:%d: add_task_line = `%s'\n", __FILE__,
            __LINE__, add_task_line.Data());
    gROOT->ProcessLine(add_task_line);

    if (mgr->InitAnalysis()) {
        mgr->StartAnalysis("grid");
    }
}
