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


    FILE *fp = fopen(config_filename, "r");
    char line[4096];
    TString emcal_correction_filename;

    while (fgets(line, 4096, fp) != NULL) {
        if (line[0] == '#') {
            continue;
        }

        char key[4096];
        char dummy[4096];
        char value[4096];

        key[0] = '\0';
        value[0] = '\0';
        sscanf(line, "%[^:]:%[ \t]%4096[^\n]", key, dummy, value);
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
        else if (strcmp(key, "runNumber") == 0) {
            for (const char *v = value; *v != '\0';) {
                while (*v != '\0' && !isdigit(*v)) {
                    v++;
                }

                int n;

                sscanf(v, "%d", &n);
                plugin->AddRunNumber(n);
                while (*v != '\0' && isdigit(*v)) {
                    v++;
                }
            }
        }
        else if (strcmp(key, "runPrefix") == 0) {
            plugin->SetRunPrefix(value);
        }
        else if (strcmp(key, "emcalCorrection") == 0) {
            emcal_correction_filename = value;
        }
    }
    fclose(fp);

    TString mklml_filename;

    if (!gSystem->AccessPathName("libmklml_gnu_so")) {
        mklml_filename = "libiomp5_so libmklml_gnu_so ";
    }
    plugin->SetAdditionalLibs(
        "AliAnalysisTaskNTGJ.h "
        "AliAnalysisTaskNTGJ.cxx "
        "special_function.h mc_truth.h emcal.h isolation.h jet.h "
        "bad_channel.h "
        "eLut.cpp eLut.h half.cpp halfExport.h halfFunction.h "
        "half.h halfLimits.h toFloat.h "
        "keras_model.h keras_model.cc "
        "photon_discr.model "
        "libCGAL.so libfastjet.so libsiscone.so "
        "libsiscone_spherical.so libfastjetplugins.so "
        "libfastjetcontribfragile.so " +
        mklml_filename +
        emcal_correction_filename);
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
    gROOT->Macro("macros/AddAliAnalysisTaskNTGJ.C");

    if (mgr->InitAnalysis()) {
        mgr->StartAnalysis("grid");
    }
}
