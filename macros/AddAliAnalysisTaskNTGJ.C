#ifdef __CLING__
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C>
#include <PWG/EMCAL/macros/AddTaskEmcalCorrectionTask.C>
#include <PWGPP/EMCAL/macros/ConfigureEMCALRecoUtils.C>
#include <OADB/macros/AddTaskPhysicsSelection.C>
#endif // __CLING__

AliAnalysisTaskNTGJ *
AddAliAnalysisTaskNTGJ(TString name,
                       TString emcal_correction_filename,
                       bool mult_selection,
                       bool physics_selection,
                       bool physics_selection_mc_analysis,
                       bool physics_selection_pileup_cut,
                       TString emcal_geometry_filename,
                       TString emcal_local2master_filename,
                       bool force_ue_subtraction,
                       double skim_cluster_min_e,
                       double skim_track_min_pt,
                       double skim_muon_track_min_pt,
                       double skim_jet_min_pt_1,
                       double skim_jet_min_pt_2,
                       double skim_jet_min_pt_3,
                       double skim_jet_average_pt,
                       int skim_multiplicity_tracklet_min_n,
                       double stored_track_min_pt,
                       double stored_jet_min_pt_raw,
                       unsigned int nrandom_isolation)
{
    AliAnalysisManager *mgr =
        AliAnalysisManager::GetAnalysisManager();

    // gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/"
    //               "AddTaskCentrality.C");

    // AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
    // if (useMC) taskCentrality->SetMCInput();

    fprintf(stderr, "%s:%d: skim_cluster_min_e = %f\n", __FILE__,
            __LINE__, skim_cluster_min_e);

    if (mult_selection) {
#ifndef __CLING__
        gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/"
                         "macros/AddTaskMultSelection.C");
#endif // __CLING__

        AliMultSelectionTask *mult_selection_task =
            AddTaskMultSelection(kFALSE);
    }

    if (emcal_correction_filename != "") {
#ifndef __CLING__
        gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/"
                         "AddTaskEmcalCorrectionTask.C");
#endif // __CLING__

        AliEmcalCorrectionTask *correction_task =
            AddTaskEmcalCorrectionTask();

        correction_task->SetUserConfigurationFilename(
            emcal_correction_filename.Data());
        correction_task->Initialize();
    }

    AliAnalysisTaskNTGJ *task =
        new AliAnalysisTaskNTGJ(name.Data());

    task->SetEMCALGeometryFilename(emcal_geometry_filename);
    task->SetEMCALLocal2MasterFilename(emcal_local2master_filename);

    if (force_ue_subtraction) {
        task->SetForceUESubtraction(force_ue_subtraction);
    }

    task->SetSkimClusterMinE(skim_cluster_min_e);
    task->SetSkimTrackMinPt(skim_track_min_pt);
    task->SetSkimMuonTrackMinPt(skim_muon_track_min_pt);
    task->SetSkimJetMinPt(skim_jet_min_pt_1, skim_jet_min_pt_2,
                          skim_jet_min_pt_3);
    task->SetSkimJetAveragePt(skim_jet_average_pt);
    task->SetSkimMultiplicityTrackletMinN(
        skim_multiplicity_tracklet_min_n);
    task->SetStoredTrackMinPt(stored_track_min_pt);
    task->SetStoredJetMinPtRaw(stored_jet_min_pt_raw);
    task->SetNRandomIsolation(nrandom_isolation);

    AliEMCALRecoUtils *reco_util = task->GetEMCALRecoUtils();
  
#ifndef __CLING__
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EMCAL/macros/"
                     "ConfigureEMCALRecoUtils.C");
#endif // __CLING__
  
    ConfigureEMCALRecoUtils(reco_util, kFALSE, kTRUE, kTRUE,
                            kFALSE, kFALSE, kFALSE); 

    reco_util->SetNumberOfCellsFromEMCALBorder(0);
    reco_util->SwitchOnRecalibration();
    reco_util->SwitchOnRunDepCorrection();

    if (physics_selection) {
#ifndef __CLING__
        gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/"
                         "AddTaskPhysicsSelection.C");
#endif // __CLING__

        AliPhysicsSelectionTask* physics_selection_task =
            AddTaskPhysicsSelection(physics_selection_mc_analysis,
                                    physics_selection_pileup_cut);
    }

    mgr->AddTask(task);
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

    TString filename = mgr->GetCommonFileName();

    filename += ":AliAnalysisTaskNTGJ";

    mgr->ConnectOutput(task, 1,
                       mgr->CreateContainer("tree", TTree::Class(),
                                            AliAnalysisManager::
                                            kOutputContainer,
                                            filename.Data()));

    AliAnalysisAlien *plugin =
        static_cast<AliAnalysisAlien *>(mgr->GetGridHandler());

    if (plugin != NULL) {
        task->SetAliROOTVersion(plugin->GetAliROOTVersion());
        task->SetAliPhysicsVersion(plugin->GetAliPhysicsVersion());
        task->SetGridDataDir(plugin->GetGridDataDir());
        task->SetGridDataPattern(plugin->GetDataPattern());
    }

    return task;
}
