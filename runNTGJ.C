#!/usr/bin/gawk {system("root -b -q " FILENAME); exit;}
// -*- mode: c++; -*-

void runNTGJ(const char *run_mode = "full")
{
	gROOT->ProcessLine(".include $ROOTSYS/include");
	gROOT->ProcessLine(".include $ALICE_ROOT/include");
	gROOT->ProcessLine(".include $ALICE_ROOT/../../fastjet/"
					   "v3.2.1_1.024-alice1-1/include");

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

    gSystem->Load("libCGAL.so");
    gSystem->Load("libfastjet.so");
    gSystem->Load("libsiscone.so");
    gSystem->Load("libsiscone_spherical.so");
    gSystem->Load("libfastjetplugins.so");
    gSystem->Load("libfastjetcontribfragile.so");

	gROOT->ProcessLine(".L AliAnalysisTaskNTGJ.cxx+g");

	AliAnalysisManager *mgr = new AliAnalysisManager();

	gROOT->Macro("$ALICE_ROOT/ANALYSIS/macros/train/"
				 "AddESDHandler.C");
	gROOT->Macro("$ALICE_ROOT/ANALYSIS/macros/train/"
				 "AddMCHandler.C");

	AliAnalysisAlien *plugin =
		new AliAnalysisAlien("pluginNTGJ");

	plugin->SetGridWorkingDir("workdir");
	plugin->SetGridOutputDir("outputdir");

	plugin->SetAliROOTVersion("v5-09-04-1");
	plugin->SetAliPhysicsVersion("v5-09-04-01-1");
	plugin->AddExternalPackage("fastjet::v3.2.1_1.024-alice1-3");
	plugin->AddIncludePath("-I. -I$ALICE_ROOT/include "
						   "-I$ALICE_ROOT/../../fastjet/"
						   "v3.2.1_1.024-alice1-3/include");

	plugin->SetAdditionalLibs(
		"AliAnalysisTaskNTGJ.h "
		"AliAnalysisTaskNTGJ.cxx "
		"bad_channel.h "
		"eLut.cpp eLut.h half.cpp halfExport.h halfFunction.h "
		"half.h halfLimits.h toFloat.h "
		"libCGAL.so libfastjet.so libsiscone.so "
		"libsiscone_spherical.so libfastjetplugins.so "
		"libfastjetcontribfragile.so");
	plugin->SetAnalysisSource("AliAnalysisTaskNTGJ.cxx");
	plugin->SetRunPrefix("");

	const int run_number_lhc15o[] = {

		245683,
#if 0
		245683, 245700, 245702, 245705, 245829, 245831, 245833,
		245949, 245952, 245954, 246001, 246003, 246037, 246042,
		246052, 246053, 246087, 246089, 246113, 246115,
#endif

		-1
	};

	const int run_number_lhc16h3_bis[] = {
 	
		244340,

#if 0
		244343, 244351, 244355, 244359, 244364, 244377, 244411,
		244416, 244418, 244421, 244453, 244456, 244480, 244481,
		244482, 244483, 244484, 244531, 244540, 244542, 244617,
		244618, 244619, 244626, 244627, 244628,
#endif

		-1
	};

	const int *run_number;

	// plugin->SetGridDataDir("/alice/sim/2016/LHC16h3_bis/1");
	// plugin->SetDataPattern("*/*/AliESDs.root");
	// run_number = run_number_lhc16h3_bis;

	plugin->SetGridDataDir("/alice/data/2015/LHC15o");
	plugin->SetDataPattern("/pass1/*/AliESDs.root");
	plugin->SetRunPrefix("000");
	run_number = run_number_lhc15o;

	for (const int *r = run_number; *r != -1; r++) {
		plugin->AddRunNumber(*r);
	}

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
	gROOT->Macro("macros/AddAliAnalysisTaskNTGJ.C");
	if (mgr->InitAnalysis()) {
		mgr->StartAnalysis("grid");
	}
}
