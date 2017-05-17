AliAnalysisTaskNTGJ *
AddAliAnalysisTaskNTGJ(TString name = "AliAnalysisTaskNTGJ")
{
	AliAnalysisManager *mgr =
		AliAnalysisManager::GetAnalysisManager();
	AliAnalysisTaskNTGJ *task =
		new AliAnalysisTaskNTGJ(name.Data());

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
					   mgr->CreateContainer("tree", TList::Class(),
											AliAnalysisManager::
											kOutputContainer,
											filename.Data()));

	return task;
}
