AliAnalysisTaskNTGJ *
AddAliAnalysisTaskNTGJ(TString name = "AliAnalysisTaskNTGJ")
{
	AliAnalysisManager *mgr =
		AliAnalysisManager::GetAnalysisManager();
	AliAnalysisTaskNTGJ *task =
		new AliAnalysisTaskNTGJ(name.Data());

	// gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/"
	// 				 "AddTaskCentrality.C");

	// AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
	// if (useMC) taskCentrality->SetMCInput();

	gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/"
					 "macros/AddTaskMultSelection.C");

	AliMultSelectionTask *mult_selection_task = AddTaskMultSelection(kFALSE);

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

	return task;
}
