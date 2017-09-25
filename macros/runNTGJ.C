#if 0
set +o posix; function join_by { local d=$1; shift; echo -n "$1";
    shift; printf "%s" "${@/#/$d}"; }
root=root; exec $root -l -b -q "$0($(join_by \",\" \"$*\" | \
/bin/sed s/\"\\\([0-9]\\+\\\)\"/\\1/g\;s/^\"\"\$//))"; exit 0
#endif

#include <TROOT.h>
#include <TSystem.h>

void runNTGJ(const char *run_mode = "full")
{
    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
    gROOT->ProcessLine(".include $ALICE_PHYSICS/include");

    const char *package[] = {
        "fastjet::v3.2.1_1.024-alice1-4",
        // Compiling against CGAL requires explicit include paths to
        // Boost, MPFR, and GMP
        "cgal::v4.6.3-18",
        "boost::v1.59.0-14",
        "MPFR::v3.1.3-4",
        "GMP::v6.0.0-2",
        // The grid ROOT package tend to lack a GCC 4.9.x =
        // CXXABI_1.3.8 dependency (see also
        // https://gcc.gnu.org/onlinedocs/libstdc++/manual/abi.html),
        // which causes jobs to fail at run time with the message:
        //
        //   root: /usr/lib64/libstdc++.so.6: version `CXXABI_1.3.8'
        //   not found (required by root)
        //   Output file AnalysisResults.root not found. Job FAILED !
        "GCC-Toolchain::v4.9.3-alice3-1",
        NULL
    };

    for (const char **p = package; *p != NULL; p++) {
        if (strncmp(*p, "GCC::", 5) != 0) {
            TString ps(*p);
            TString include = "$ALICE_ROOT/../../" +
                ps.ReplaceAll("::", "/") + "/include";

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

    plugin->SetAliROOTVersion("v5-09-11-1");
    plugin->SetAliPhysicsVersion("v5-09-11-01-1");
    for (const char **p = package; *p != NULL; p++) {
        plugin->AddExternalPackage(*p);
    }

    TString include_path("-I. -I$ALICE_ROOT/include "
                         "-I$ALICE_PHYSICS/include ");

    for (const char **p = package; *p != NULL; p++) {
        if (strncmp(*p, "GCC::", 5) != 0) {
            include_path += "-I$ALICE_ROOT/../../" +
                TString(*p).ReplaceAll("::", "/") + "/include ";
        }
    }

    plugin->AddIncludePath(include_path);

    plugin->SetAdditionalLibs(
        "AliAnalysisTaskNTGJ.h "
        "AliAnalysisTaskNTGJ.cxx "
        "special_function.h emcal.h isolation.h jet.h bad_channel.h "
        "eLut.cpp eLut.h half.cpp halfExport.h halfFunction.h "
        "half.h halfLimits.h toFloat.h "
        "keras_model.h keras_model.cc "
        "libCGAL.so libfastjet.so libsiscone.so "
        "libsiscone_spherical.so libfastjetplugins.so "
        "libfastjetcontribfragile.so");
    plugin->SetAnalysisSource("AliAnalysisTaskNTGJ.cxx");
    plugin->SetRunPrefix("");

	// 5 TeV PbPb

    const int run_number_lhc15o[] = {

        245683,
#if 0
        245683, 245700, 245702, 245705, 245829, 245831, 245833,
        245949, 245952, 245954, 246001, 246003, 246037, 246042,
        246052, 246053, 246087, 246089, 246113, 246115,
#endif

        -1
    };

	// 8 TeV PYTHIA 8 QCD

    const int run_number_lhc16c2[] = {
    
        180720,
#if 0
        182692, 184215, 185687, 187488, 189616, 190393, 192073,
        192349, 193051
#endif

        -1
    };

	const int run_number_lhc16k[] = {

		258537,

#if 0
		258499, 258477, 258456, 258454, 258426, 258393, 258387, 258359, 258336, 258299, 258278, 258274, 258273, 258271, 258270, 258258, 258257, 258256, 258204, 258203, 258202, 258198, 258197, 258178, 258117, 258114, 258113, 258109, 258108, 258107, 258063, 258062, 258059, 258049, 258048, 258045, 258042, 258019, 258017, 258014, 258012, 257963, 257960, 257958, 257957, 257939, 257937, 257936, 257893, 257892, 257855, 257850, 257803, 257800, 257799, 257798, 257797, 257773, 257765, 257754, 257737, 257735, 257734, 257733, 257724, 257697, 257694, 257692, 257691, 257689, 257687, 257682, 257642, 257606, 257605, 257594, 257590, 257587, 257566, 257562, 257561, 257560, 257541, 257540, 257539, 257537, 257531, 257530, 257492, 257491, 257490, 257487, 257474, 257457, 257320, 257260, 257224, 257209, 257206, 257204, 257145, 257144, 257142, 257141, 257140, 257139, 257138, 257137, 257136, 257100, 257092, 257084, 257083, 257082, 257080, 257077, 257026, 257021, 257012, 257011, 256944, 256942, 256941, 256697, 256695, 256694, 256692, 256691, 256684, 256681, 256677, 256676, 256658, 256620, 256619, 256592, 256591, 256589, 256567, 256565, 256564, 256562, 256561, 256560, 256556, 256554, 256552, 256514, 256512, 256510, 256506, 256504,
#endif

		-1

	};

	// 5 TeV PYTHIA + HIJING QCD

    const int run_number_lhc16h2a_bis[] = {
    
        246994,

        -1
    };

    const int *run_number;

    plugin->SetGridDataDir("/alice/sim/2016/LHC16c2/16");
    plugin->SetDataPattern("*/*/AliESDs.root");
    run_number = run_number_lhc16c2;

    // plugin->SetGridDataDir("/alice/data/2015/LHC15o");
    // plugin->SetDataPattern("/pass1/*/AliESDs.root");
    // plugin->SetRunPrefix("000");
    // run_number = run_number_lhc15o;

    // plugin->SetGridDataDir("/alice/data/2016/LHC16k");
    // plugin->SetDataPattern("/pass1/*/AliESDs.root");
    // plugin->SetRunPrefix("000");
    // run_number = run_number_lhc16k;

    for (const int *r = run_number; *r != -1; r++) {
        plugin->AddRunNumber(*r);
    }

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
    gROOT->Macro("macros/AddAliAnalysisTaskNTGJ.C");
 
    if (mgr->InitAnalysis()) {
        mgr->StartAnalysis("grid");
    }
}
