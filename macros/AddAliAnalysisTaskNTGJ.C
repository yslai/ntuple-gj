AliAnalysisTaskNTGJ *
AddAliAnalysisTaskNTGJ(TString name = "AliAnalysisTaskNTGJ",
                       TString emcal_correction_filename =
                       "emcal_correction.yaml",
                       // 1e+309 = INFINITY (for IEEE 754 double)
                       double skim_cluster_min_e = -1e+309,
                       double skim_track_min_pt = -1e+309,
                       double skim_muon_track_min_pt = -1e+309,
                       double skim_jet_min_pt_1 = -1e+309,
                       double skim_jet_min_pt_2 = -1e+309,
                       double skim_jet_min_pt_3 = -1e+309,
                       // -2147483648 = INT_MIN
                       int skim_multiplicity_tracklet_min_n =
                       -2147483648,
                       double stored_track_min_pt = -1e+309,
                       double stored_jet_min_pt_raw = -1e+309)
{
    AliAnalysisManager *mgr =
        AliAnalysisManager::GetAnalysisManager();
    AliAnalysisTaskNTGJ *task =
        new AliAnalysisTaskNTGJ(name.Data());

    // gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/"
    //               "AddTaskCentrality.C");

    // AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
    // if (useMC) taskCentrality->SetMCInput();

    fprintf(stderr, "%s:%d: skim_cluster_min_e = %f\n", __FILE__,
            __LINE__, skim_cluster_min_e);

    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/"
                     "macros/AddTaskMultSelection.C");

    AliMultSelectionTask *mult_selection_task = AddTaskMultSelection(kFALSE);

    if (emcal_correction_filename != "") {
        gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/"
                         "AddTaskEmcalCorrectionTask.C");

        AliEmcalCorrectionTask *correction_task =
            AddTaskEmcalCorrectionTask();

        correction_task->SetUserConfigurationFilename(
            emcal_correction_filename.Data());
        correction_task->Initialize();
    }

    AliEMCALRecoUtils *reco_util = task->GetEMCALRecoUtils();
  
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EMCAL/macros/"
                     "ConfigureEMCALRecoUtils.C");
  
    ConfigureEMCALRecoUtils(reco_util, kFALSE, kTRUE, kTRUE,
                            kFALSE, kFALSE, kFALSE); 

    reco_util->SetNumberOfCellsFromEMCALBorder(0);
    reco_util->SwitchOnRecalibration();
    reco_util->SwitchOnRunDepCorrection();
  
    TString filename = mgr->GetCommonFileName();

    filename += ":AliAnalysisTaskNTGJ";

    mgr->AddTask(task);
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1,
                       mgr->CreateContainer("tree", TTree::Class(),
                                            AliAnalysisManager::
                                            kOutputContainer,
                                            filename.Data()));

    AliAnalysisAlien *plugin = mgr->GetGridHandler();

    if (plugin != NULL) {
        task->SetAliROOTVersion(plugin->GetAliROOTVersion());
        task->SetAliPhysicsVersion(plugin->GetAliPhysicsVersion());
        task->SetGridDataDir(plugin->GetGridDataDir());
        task->SetGridDataPattern(plugin->GetDataPattern());

        task->SetSkimClusterMinE(skim_cluster_min_e);
        task->SetSkimTrackMinPt(skim_track_min_pt);
        task->SetSkimMuonTrackMinPt(skim_muon_track_min_pt);
        task->SetSkimJetMinPt(skim_jet_min_pt_1, skim_jet_min_pt_2,
                              skim_jet_min_pt_3);
        task->SetSkimMultiplicityTrackletMinN(
            skim_multiplicity_tracklet_min_n);
        task->SetStoredTrackMinPt(stored_track_min_pt);
        task->SetStoredJetMinPtRaw(stored_jet_min_pt_raw);
    }

    return task;
}
