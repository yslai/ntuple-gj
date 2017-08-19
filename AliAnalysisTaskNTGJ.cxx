#include <climits>

// This is purely to suppress warning inside ROOT once instanciating
// std::vector<bool>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#include <TCollectionProxyInfo.h>
#pragma GCC diagnostic pop

#include <Compression.h>
#include <TFile.h>
#include <TGeoManager.h>
#include <TGeoGlobalMagField.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#include <AliVEvent.h>
#include <AliAODEvent.h>
#pragma GCC diagnostic pop
#include <AliESDEvent.h>
#include <AliESDHeader.h>
#include <AliVTrack.h>
#include <AliVCluster.h>
#include <AliVCaloCells.h>
#include <AliOADBContainer.h>
#include <AliAnalysisManager.h>
#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include <AliGenEventHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliStack.h>
#include <AliAODMCParticle.h>
#include <AliESDMuonTrack.h>

#include <AliMultSelection.h>
#include <AliEventplane.h>

#include <AliMagF.h>

#include "AliAnalysisTaskNTGJ.h"

#ifndef __CINT__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "eLut.cpp"
#include "half.cpp"
#pragma GCC diagnostic pop
#include "special_function.h"
#include "emcal.h"
#include "jet.h"
#include "isolation.h"
#endif // __CINT__

// FIXME: This needs to be moved to somewhere else, later on

namespace {

    class alice_jec_t {
    protected:
        const char *_version;
    public:
        alice_jec_t(void)
            : _version("lhc16c2-2017-06-pre1")
        {
        }
        TLorentzVector jec_p(TLorentzVector p) const
        {
            const double pt = p.Pt();

            // Dummy value
            const double pt_calib = 2.0 * pt;

            const double eta = p.Eta();

            // Mean value for jet pT raw within 2.5-18 GeV
            static const double par_eta[] = { 0.015, 7.7, 0.65 };

            const double eta_calib = par_eta[0] *
                (TMath::Erf(par_eta[1] * (eta - par_eta[2])) +
                 TMath::Erf(par_eta[1] * (eta + par_eta[2])));

            TLorentzVector calib;

            calib.SetPtEtaPhiM(pt_calib, eta_calib, p.Phi(), p.M());

            return calib;
        }
        const char *version(void) const
        {
            return _version;
        }
    };

}

ClassImp(AliAnalysisTaskNTGJ);

#define EMCAL_GEOMETRY_NAME "EMCAL_COMPLETE12SMV1_DCAL_8SM"

#define ntrigger_class NTRIGGER_CLASS_MAX

#define B CHAR_MIN
#define b 0
#define S SHRT_MIN
#define s 0
#define I INT_MIN
#define i 0
#define F NAN
#define D NAN
#define L LONG_MIN
#define l 0
#define O false
#define BRANCH(b, t)                            \
    _branch_ ## b((t)),
#define BRANCH_ARRAY(b, d, t)
#define BRANCH_ARRAY2(b, d, e, t)
#define BRANCH_STR(b)
#define BRANCH_STR_ARRAY(b, d)                  \
    _branch_ ## b("TObjArray", (d)),

#define CLASS_INITIALIZATION                                \
    _emcal_geometry_name(EMCAL_GEOMETRY_NAME),              \
    _tree_event(NULL),                                      \
    MEMBER_BRANCH                                           \
    _run_number_current(INT_MIN),                           \
    _f1_ncluster_tpc_linear_pt_dep(NULL),                   \
    _track_cut(std::vector<AliESDtrackCuts>()),             \
    _reco_util(new AliEMCALRecoUtils),                      \
    _emcal_geometry(NULL),                                  \
    _muon_track_cut(new AliMuonTrackCuts),                  \
    _ncell(EMCAL_NCELL),                                    \
    _skim_cluster_min_e(-INFINITY),                         \
    _skim_track_min_pt(-INFINITY),                          \
    _skim_muon_track_min_pt(-INFINITY),                     \
    _skim_jet_min_pt(std::vector<double>(3, -INFINITY)),    \
    _skim_multiplicity_tracklet_min_n(INT_MIN),             \
    _stored_track_min_pt(-INFINITY),                        \
    _stored_jet_min_pt_raw(-INFINITY),                      \
    _emcal_mask(std::vector<bool>()),                       \
    _emcal_cell_position(NULL),                             \
    _emcal_cell_area(std::vector<double>()),                \
    _emcal_cell_incident(std::vector<std::set<size_t> >()), \
    _alien_plugin(NULL),                                    \
    _metadata_filled(false)

AliAnalysisTaskNTGJ::AliAnalysisTaskNTGJ(void)
    : AliAnalysisTaskSE(), CLASS_INITIALIZATION
{
}

AliAnalysisTaskNTGJ::AliAnalysisTaskNTGJ(
    const char *name)
    : AliAnalysisTaskSE(name), CLASS_INITIALIZATION
{
    DefineOutput(1, TTree::Class());
}

AliAnalysisTaskNTGJ::AliAnalysisTaskNTGJ(
    const AliAnalysisTaskNTGJ &x)
    : AliAnalysisTaskSE(), CLASS_INITIALIZATION
{
}

#undef ntrigger_class

#undef BRANCH
#undef BRANCH_ARRAY
#undef BRANCH_ARRAY2
#undef BRANCH_STR
#undef BRANCH_STR_ARRAY
#undef B
#undef b
#undef S
#undef s
#undef I
#undef i
#undef F
#undef D
#undef L
#undef l
#undef O

AliAnalysisTaskNTGJ &AliAnalysisTaskNTGJ::operator=(
    const AliAnalysisTaskNTGJ &x)
{
    return *this;
}

AliAnalysisTaskNTGJ::~AliAnalysisTaskNTGJ(void)
{
    if (_tree_event != NULL) {
        delete _tree_event;
    }
    if (_emcal_cell_position != NULL) {
        delete reinterpret_cast<std::vector<point_2d_t> *>
            (_emcal_cell_position);
    }

    delete _reco_util;
}

void AliAnalysisTaskNTGJ::UserCreateOutputObjects(void)
{
    TFile *file = OpenFile(1);

    if (file != NULL) {
        file->SetCompressionSettings(ROOT::CompressionSettings(
            ROOT::kZLIB, 9));
    }

    /////////////////////////////////////////////////////////////////

    _tree_event = new TTree("_tree_event", "");

#define BRANCH(b, t)                    \
    _tree_event->Branch(                \
        #b, &_branch_ ## b, #b "/" #t);
#define BRANCH_ARRAY(b, d, t)                   \
    _tree_event->Branch(                        \
        #b, _branch_ ## b, #b "[" #d "]/" #t);
#define BRANCH_ARRAY2(b, d, e, t)                       \
    _tree_event->Branch(                                \
        #b, _branch_ ## b, #b "[" #d "][" #e "]/" #t);
#define BRANCH_STR(b)                           \
    _tree_event->Branch(                        \
        #b, _branch_ ## b, #b "/C");
#define BRANCH_STR_ARRAY(b, d)                  \
    _tree_event->Branch(                        \
        #b, &_branch_ ## b, #b "[" #d "]/C");

    MEMBER_BRANCH;

#undef BRANCH
#undef BRANCH_ARRAY
#undef BRANCH_ARRAY2
#undef BRANCH_STR
#undef BRANCH_STR_ARRAY

    PostData(1, _tree_event);

    /////////////////////////////////////////////////////////////////
}

#undef MEMBER_BRANCH

void AliAnalysisTaskNTGJ::UserExec(Option_t *option)
{
    if (_emcal_mask.size() != EMCAL_NCELL) {
        _emcal_mask.resize(EMCAL_NCELL);
#if 1 // Keep = 1 to for an actual EMCAL mask (and not all channels
      // turned on)
        for (unsigned int i = 0; i < EMCAL_NCELL; i++) {
            _emcal_mask[i] = inside_edge(i, 1);
        }
#include "bad_channel.h"
        for (unsigned int i = 0; bad_channel_emcal[i] != -1; i++) {
            if (inside_edge(bad_channel_emcal[i], 1)) {
                unsigned int bad_cell_3_3[9];

                cell_3_3(bad_cell_3_3, bad_channel_emcal[i]);
                for (size_t j = 0; j < 9; j++) {
                    _emcal_mask[bad_cell_3_3[j]] = false;
                }
            }
        }
#else
        for (unsigned int i = 0; i < EMCAL_NCELL; i++) {
            _emcal_mask[i] = true;
        }
#endif
    }

    AliVEvent *event = InputEvent();

    if (event == NULL) {
        return;
    }

#if 1
    if (event->GetRunNumber() != _run_number_current) {
        _run_number_current = event->GetRunNumber();
        if (_muon_track_cut != NULL) {
            _muon_track_cut->SetAllowDefaultParams(kTRUE);
            _muon_track_cut->SetRun(
                (AliInputEventHandler *)
                ((AliAnalysisManager::GetAnalysisManager())->
                 GetInputEventHandler()));
        }
    }
#endif

    AliESDEvent *esd_event = dynamic_cast<AliESDEvent *>(event);
    AliAODEvent *aod_event = NULL;

    if (esd_event == NULL) {
        aod_event = dynamic_cast<AliAODEvent *>(event);
        if (aod_event == NULL) {
            return;
        }
    }

    alice_jec_t jec;

    if (!_metadata_filled) {
        // Use gitattributes ident mechanism to track the blob object
        // name
        strncpy(_branch_id_git, "$Id$", BUFSIZ);
        _branch_id_git[BUFSIZ - 1] = '\0';
        strncpy(_branch_version_jec, jec.version(), BUFSIZ);
        _branch_version_jec[BUFSIZ - 1] = '\0';
        if (esd_event != NULL) {
            for (size_t i = 0; i < 2; i++) {
                _branch_beam_particle[i] =
                    esd_event->GetBeamParticle(i);
            }

            const AliESDRun *esd_run = esd_event->GetESDRun();

            if (esd_run != NULL) {
                _branch_trigger_class.Delete();
                for (size_t i = 0; i < NTRIGGER_CLASS_MAX; i++) {
                    new (_branch_trigger_class[i])
                        TObjString(esd_run->GetTriggerClass(i));
                }
            }
        }
        else if (aod_event != NULL) {
            // FIXME: AOD not handled
            std::fill(_branch_beam_particle,
                      _branch_beam_particle + 2, 0);
        }

        _metadata_filled = true;
    }
    else {
        // To make sure no space is wasted, set metadata in all
        // subsequent events to empty
        _branch_id_git[0] = '\0';
        _branch_version_aliroot[0] = '\0';
        _branch_version_aliphysics[0] = '\0';
        _branch_version_jec[0] = '\0';
        _branch_grid_data_dir[0] = '\0';
        _branch_grid_data_pattern[0] = '\0';
        std::fill(_branch_beam_particle,
                  _branch_beam_particle + 2, 0);
        _branch_trigger_class.Delete();
    }

    _branch_run_number = event->GetRunNumber();
    _branch_trigger_mask[0] = esd_event != NULL ?
        esd_event->GetTriggerMask() :
        aod_event->GetTriggerMask();
    _branch_trigger_mask[1] = esd_event != NULL ?
        esd_event->GetTriggerMaskNext50() :
        aod_event->GetTriggerMaskNext50();
    _branch_has_misalignment_matrix = false;

    if (_f1_ncluster_tpc_linear_pt_dep == NULL) {
        // PWG-JE cuts, see AddTrackCutsLHC10h() in
        // AliPhysics/PWGJE/macros/AddTaskESDFilterPWGJETrain.C

        _f1_ncluster_tpc_linear_pt_dep =
            new TFormula("_f1_ncluster_tpc_linear_pt_dep",
                         "70.+30./20.*x");

        for (size_t i = 0; i < 2; i++) {
            _track_cut.push_back(AliESDtrackCuts("AliESDtrackCuts"));
            _track_cut.back().
                SetMinNClustersTPCPtDep(_f1_ncluster_tpc_linear_pt_dep,
                                        20.0);
            _track_cut.back().SetMinNClustersTPC(70);
            _track_cut.back().SetMaxChi2PerClusterTPC(4);
            _track_cut.back().SetRequireTPCStandAlone(kTRUE);
            _track_cut.back().SetAcceptKinkDaughters(kFALSE);
            _track_cut.back().SetRequireTPCRefit(kTRUE);
            _track_cut.back().SetMaxFractionSharedTPCClusters(0.4);
            _track_cut.back().SetRequireITSRefit(kTRUE);
            _track_cut.back().SetMaxDCAToVertexXY(2.4);
            _track_cut.back().SetMaxDCAToVertexZ(3.2);
            _track_cut.back().SetDCAToVertex2D(kTRUE);
            _track_cut.back().SetMaxChi2PerClusterITS(36);
            _track_cut.back().SetMaxChi2TPCConstrainedGlobal(36);
            _track_cut.back().SetRequireSigmaToVertex(kFALSE);
            _track_cut.back().SetEtaRange(-0.9, 0.9);
            _track_cut.back().SetPtRange(0.15, 1e+15);

            if (i == 0) {
                _track_cut.back().
                    SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                             AliESDtrackCuts::kAny);
            }
            else {
                _track_cut.back().SetRequireITSRefit(kFALSE);
            }
        }

        // "2015 PbPb" cuts, see GetStandardITSTPCTrackCuts2015PbPb()
        // in AliRoot/ANALYSIS/ANALYSISalice/AliESDtrackCuts.cxx .
        // Both clusterCut = 0 or 1 cases are kept, but the options
        // selPrimaries = kTRUE, cutAcceptanceEdges = kTRUE,
        // removeDistortedRegions = kFALSE are fixed

        for (size_t i = 0; i < 2; i++) {
            _track_cut.push_back(AliESDtrackCuts("AliESDtrackCuts"));

            if (i == 0) {
                _track_cut.back().SetMinNClustersTPC(50);
            }
            else {
                _track_cut.back().SetMinNCrossedRowsTPC(70);
                _track_cut.back().
                    SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
            }
            // cutAcceptanceEdges == kTRUE
            _track_cut.back().SetCutGeoNcrNcl(2., 130., 1.5, 0.0, 0.0);
            _track_cut.back().SetMaxChi2PerClusterTPC(4);
            _track_cut.back().SetAcceptKinkDaughters(kFALSE);
            _track_cut.back().SetRequireTPCRefit(kTRUE);
            _track_cut.back().SetRequireITSRefit(kTRUE);
            _track_cut.back().
                    SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                             AliESDtrackCuts::kAny);
            // selPrimaries == kTRUE
            _track_cut.back().
                SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
            // selPrimaries == kTRUE
            _track_cut.back().SetMaxChi2TPCConstrainedGlobal(36);
            _track_cut.back().SetMaxDCAToVertexZ(2);
            _track_cut.back().SetDCAToVertex2D(kFALSE);
            _track_cut.back().SetRequireSigmaToVertex(kFALSE);
            _track_cut.back().SetMaxChi2PerClusterITS(36);
        }
    }

    AliVVZERO *v0 = event->GetVZEROData();

    if (v0 != NULL) {
        for (size_t i = 0; i < 64; i++) {
            _branch_multiplicity_v0[i] = v0->GetMultiplicity(i);
        }
    }

    static const char *centrality_method[9] = {
        "V0M", "CL0", "CL1",
        "V0Mplus05", "V0Mplus10", "V0Mminus05", "V0Mminus10",
        "SPDClustersCorr", "SPDTracklets"
    };

    AliMultSelection *mult_selection = static_cast<AliMultSelection *>
        (event->FindListObject("MultSelection"));

    std::fill(_branch_centrality, _branch_centrality + 9, NAN);
    if (mult_selection != NULL) {
        for (size_t i = 0; i < 9; i++) {
            _branch_centrality[i] =
                mult_selection->GetMultiplicityPercentile(
                    centrality_method[i]);
        }
    }
    else {
        AliCentrality *centrality = event->GetCentrality();

        if (centrality != NULL) {
            for (size_t i = 0; i < 9; i++) {
                _branch_centrality[i] =
                    centrality->GetCentralityPercentile(
                        centrality_method[i]);
            }
        }
    }
    // Copy for easy cutting (where duplicated global variable are not
    // very wasteful)
    _branch_centrality_v0m = _branch_centrality[0];

    std::fill(_branch_event_plane_psi_v0,
              _branch_event_plane_psi_v0 + 3, NAN);

    AliEventplane *event_plane = event->GetEventplane();

    if (event_plane != NULL) {
        for (size_t i = 1; i < 4; i++) {
            _branch_event_plane_psi_v0[i - 1] =
                event->GetEventplane()->CalculateVZEROEventPlane(
                    event, 10, i,
                    _branch_event_plane_q_v0[i - 1][0],
                    _branch_event_plane_q_v0[i - 1][1]);
        }
    }

    if (_emcal_geometry == NULL) {
        TGeoManager::Import(
            "$ALICE_PHYSICS/OADB/EMCAL/geometry_2015.root");
        _emcal_geometry =
            AliEMCALGeometry::GetInstance(_emcal_geometry_name);

        AliOADBContainer emcal_geometry_container("emcal");
        emcal_geometry_container.InitFromFile(
            "$ALICE_PHYSICS/OADB/EMCAL/EMCALlocal2master.root",
            "AliEMCALgeo");
        TObjArray *geometry_matrix = dynamic_cast<TObjArray *>(
            emcal_geometry_container.GetObject(
                _branch_run_number, "EmcalMatrices"));
        if (geometry_matrix != NULL) {
            const Int_t nsm = _emcal_geometry->GetEMCGeometry()->
                GetNumberOfSuperModules();

            for (Int_t sm = 0; sm < nsm; sm++) {
                _emcal_geometry->SetMisalMatrix(
                    dynamic_cast<TGeoHMatrix *>(
                        geometry_matrix->At(sm)),
                    sm);
             }
             _branch_has_misalignment_matrix = true;
        }

        _emcal_cell_position = new std::vector<point_2d_t>();
        for (int cell_id = 0; cell_id < EMCAL_NCELL; cell_id++) {
            TVector3 v;

            _emcal_geometry->GetGlobal(cell_id, v);
            _branch_cell_eta[cell_id] = v.Eta();
            _branch_cell_phi[cell_id] = v.Phi();
            reinterpret_cast<std::vector<point_2d_t> *>
                (_emcal_cell_position)->push_back(
                    point_2d_t(v.Eta(), v.Phi()));
        }

        voronoi_area_incident(
            _emcal_cell_area, _emcal_cell_incident,
            *reinterpret_cast<std::vector<point_2d_t> *>
            (_emcal_cell_position));

        double sum_area_inside = 0;
        size_t count_inside = 0;

        for (int cell_id = 0; cell_id < EMCAL_NCELL; cell_id++) {
            if (inside_edge(cell_id, 1)) {
                sum_area_inside += _emcal_cell_area[cell_id];
                count_inside++;
            }
        }

        const double mean_area_inside = sum_area_inside / count_inside;

        for (int cell_id = 0; cell_id < EMCAL_NCELL; cell_id++) {
            if (!inside_edge(cell_id, 1)) {
                _emcal_cell_area[cell_id] = mean_area_inside;
            }
            _branch_cell_voronoi_area[cell_id] =
                _emcal_cell_area[cell_id];
        }
    }
    else {
        for (int cell_id = 0; cell_id < EMCAL_NCELL; cell_id++) {
            _branch_cell_eta[cell_id] = NAN;
            _branch_cell_phi[cell_id] = NAN;
            _branch_cell_voronoi_area[cell_id] = NAN;
        }
    }

    AliMCEvent *mc_truth_event = MCEvent();

    if (mc_truth_event != NULL) {
        mc_truth_event->PreReadAll();
    }

    _branch_eg_signal_process_id = INT_MIN;
    _branch_eg_mpi = INT_MIN;
    _branch_eg_pt_hat = NAN;
    _branch_eg_cross_section = NAN;
    _branch_eg_weight = NAN;
    _branch_eg_ntrial = -1;

    // Not stored by ALICE SW

    _branch_eg_scale_pdf = NAN;
    _branch_eg_alpha_qcd = NAN;
    _branch_eg_alpha_qed = NAN;
    std::fill(_branch_eg_pdf_id, _branch_eg_pdf_id + 2, INT_MIN);
    std::fill(_branch_eg_pdf_x, _branch_eg_pdf_x + 2, NAN);
    std::fill(_branch_eg_pdf_x_pdf, _branch_eg_pdf_x_pdf + 2, NAN);

    // FIXME: Weight is missing, AliGenEventHeader::EventWeight()

    AliGenEventHeader *mc_truth_header = mc_truth_event != NULL ?
        mc_truth_event->GenEventHeader() : NULL;
    AliGenPythiaEventHeader *mc_truth_pythia_header;

    if (mc_truth_header != NULL) {
        _branch_eg_weight = mc_truth_header->EventWeight();

        TArrayF eg_primary_vertex(3);

        mc_truth_header->PrimaryVertex(eg_primary_vertex);

        for (Int_t i = 0; i < 3; i++) {
            _branch_eg_primary_vertex[i] = eg_primary_vertex.At(i);
        }
        mc_truth_pythia_header =
            dynamic_cast<AliGenPythiaEventHeader *>(mc_truth_header);
        if (mc_truth_pythia_header != NULL) {
            _branch_eg_signal_process_id =
                mc_truth_pythia_header->ProcessType();
            _branch_eg_mpi = mc_truth_pythia_header->GetNMPI();
            _branch_eg_pt_hat =
                mc_truth_pythia_header->GetPtHard();
            _branch_eg_cross_section =
                mc_truth_pythia_header->GetXsection();
            _branch_eg_ntrial = mc_truth_pythia_header->Trials();
        }
    }

    AliStack *stack;

    if (mc_truth_event != NULL) {
        stack = mc_truth_event->Stack();
    }

#if 0
    AliMagF *field = (AliMagF *)TGeoGlobalMagField::Instance()->GetField();
#endif
    if (esd_event != NULL) {
        esd_event->InitMagneticField();
    }
    else if (aod_event != NULL) {
        aod_event->InitMagneticField();
    }

    std::fill(_branch_primary_vertex,
              _branch_primary_vertex + 3, NAN);
    std::fill(_branch_primary_vertex_sigma,
              _branch_primary_vertex_sigma + 3, NAN);
    _branch_primary_vertex_ncontributor = INT_MIN;

    const AliVVertex *primary_vertex = event->GetPrimaryVertex();

    if (primary_vertex != NULL) {
        primary_vertex->GetXYZ(_branch_primary_vertex);
    }
    if (esd_event != NULL) {
        const AliESDVertex *esd_primary_vertex =
            esd_event->GetPrimaryVertex();

        if (esd_primary_vertex != NULL) {
            esd_primary_vertex->GetSigmaXYZ(
                _branch_primary_vertex_sigma);
            _branch_primary_vertex_ncontributor =
                esd_primary_vertex->GetNContributors();
        }
    }

    std::fill(_branch_primary_vertex_spd,
              _branch_primary_vertex_spd + 3, NAN);
    std::fill(_branch_primary_vertex_spd_sigma,
              _branch_primary_vertex_spd_sigma + 3, NAN);
    _branch_primary_vertex_spd_ncontributor = INT_MIN;

    const AliVVertex *primary_vertex_spd =
        event->GetPrimaryVertexSPD();

    if (primary_vertex_spd != NULL) {
        primary_vertex_spd->GetXYZ(_branch_primary_vertex_spd);
    }
    if (esd_event != NULL) {
        const AliESDVertex *esd_primary_vertex_spd =
            esd_event->GetPrimaryVertexSPD();

        if (esd_primary_vertex_spd != NULL) {
            esd_primary_vertex_spd->GetSigmaXYZ(
                _branch_primary_vertex_spd_sigma);
            _branch_primary_vertex_spd_ncontributor =
                esd_primary_vertex_spd->GetNContributors();
        }
    }

    _branch_npileup_vertex_spd = INT_MIN;
    _branch_pileup_vertex_spd_ncontributor = INT_MIN;
    _branch_pileup_vertex_spd_min_z_distance = INFINITY;
    if (esd_event != NULL) {
        _branch_npileup_vertex_spd =
            esd_event->GetNumberOfPileupVerticesSPD();
        for (Int_t i = 0; i < _branch_npileup_vertex_spd; i++) {
            const AliESDVertex *pileup_vertex_spd =
                esd_event->GetPileupVertexSPD(i);

            if (pileup_vertex_spd != NULL) {
                _branch_pileup_vertex_spd_ncontributor =
                    pileup_vertex_spd->GetNContributors();
                _branch_pileup_vertex_spd_min_z_distance =
                    std::min(_branch_pileup_vertex_spd_min_z_distance,
                             fabs(pileup_vertex_spd->GetZ() -
                                  _branch_primary_vertex_spd[2]));
            }
        }
    }

    if (_skim_multiplicity_tracklet_min_n > INT_MIN) {
        AliVMultiplicity *multiplicity = event->GetMultiplicity();

        if (multiplicity == NULL ||
            !(multiplicity->GetNumberOfTracklets() >=
              _skim_multiplicity_tracklet_min_n)) {
            return;
        }
    }

    TRefArray calo_cluster;

    event->GetEMCALClusters(&calo_cluster);

    if (_skim_cluster_min_e > -INFINITY) {
        double cluster_e_max = -INFINITY;

        for (Int_t i = 0; i < calo_cluster.GetEntriesFast(); i++) {
            AliVCluster *c =
                static_cast<AliVCluster *>(calo_cluster.At(i));

            if (!(c->GetNCells() > 1) ||
                cell_masked(c, _emcal_mask)) {
                continue;
            }

            TLorentzVector p;

            c->GetMomentum(p, _branch_primary_vertex);
            cluster_e_max = std::max(cluster_e_max, p.E());
        }
        if (!(cluster_e_max >= _skim_cluster_min_e)) {
            return;
        }
    }

    std::vector<size_t> stored_mc_truth_index;
    std::vector<Int_t> reverse_stored_mc_truth_index;

    if (mc_truth_event != NULL) {
        stored_mc_truth_index.resize(
            mc_truth_event->GetNumberOfTracks(), ULONG_MAX);

        size_t nmc_truth = 0;

        for (Int_t i = 0;
             i < mc_truth_event->GetNumberOfTracks(); i++) {
            if (!mc_truth_event->IsPhysicalPrimary(i)) {
                // Keep only primary final state particles
                continue;
            }
            stored_mc_truth_index[i] = nmc_truth;
            reverse_stored_mc_truth_index.push_back(i);
            nmc_truth++;
        }
    }

    std::vector<fastjet::PseudoJet> particle_reco;

    std::fill(_branch_met_tpc, _branch_met_tpc + 2, 0);

    std::vector<size_t> reco_track_index;

    _branch_ntrack = 0;
    if (esd_event != NULL) {
        for (Int_t i = 0; i < esd_event->GetNumberOfTracks(); i++) {
            AliESDtrack *t = esd_event->GetTrack(i);

            if (t == NULL) {
                continue;
            }

            // Apply PWG-JE cuts (track cuts 0 and 1)

            if (_track_cut[0].AcceptTrack(t) ||
                _track_cut[1].AcceptTrack(t)) {
                reco_track_index.push_back(_branch_ntrack);
                particle_reco.push_back(fastjet::PseudoJet(
                    t->Px(), t->Py(), t->Pz(), t->P()));
                _branch_met_tpc[0] += t->Px();
                _branch_met_tpc[1] += t->Py();
            }

            // Store tracks passing PWG-JE *or* "2015 PbPb" cuts

            if (_track_cut[0].AcceptTrack(t) ||
                _track_cut[1].AcceptTrack(t) ||
                _track_cut[2].AcceptTrack(t) ||
                _track_cut[3].AcceptTrack(t)) {
                _branch_track_e[_branch_ntrack] = half(t->E());
                _branch_track_pt[_branch_ntrack] = half(t->Pt());
                _branch_track_eta[_branch_ntrack] = half(t->Eta());
                _branch_track_phi[_branch_ntrack] =
                    half(angular_range_reduce(t->Phi()));

                // Shortened track quality bit mask. Here bit 0 and 1
                // are the PWG-JE's bit 4 and 8. Test for
                // ((track_quality[i] & 3) != 0), i being the track
                // index, to get PWG-JE's "272" (= 1 << 4 | 1 << 8)
                // cut. Test for ((track_quality[i] & 4) == 0) and
                // ((track_quality[i] & 8) == 0) for the "2015 PbPb" cut
                // with clusterCut = 0 and 1.

                _branch_track_quality[_branch_ntrack] = 0U;
                for (size_t j = 0; j < _track_cut.size(); j++) {
                    _branch_track_quality[_branch_ntrack] |=
                        _track_cut[j].AcceptTrack(t) ? 1U << j : 0;
                }

                _branch_track_tpc_dedx[_branch_ntrack] =
                    half(t->GetTPCsignal());

                static const Int_t mode_inner_wall = 1;
                static const Double_t dead_zone_width_cm = 2;
                static const Double_t max_z_cm = 220;

                _branch_track_tpc_length_active_zone
                    [_branch_ntrack] = NAN;
                if (t->GetInnerParam() != NULL) {
                    _branch_track_tpc_length_active_zone
                        [_branch_ntrack] =
                        half(t->GetLengthInActiveZone
                             (mode_inner_wall, dead_zone_width_cm,
                              max_z_cm, event->GetMagneticField()));
                }
                _branch_track_tpc_xrow[_branch_ntrack] =
                    std::min(static_cast<Float_t>(UCHAR_MAX),
                             std::max(0.0F, t->GetTPCCrossedRows()));
                _branch_track_tpc_ncluster[_branch_ntrack] =
                    std::min(UCHAR_MAX,
                             std::max(0, t->GetNumberOfTPCClusters()));
                _branch_track_tpc_ncluster_dedx[_branch_ntrack] =
                    std::min(static_cast<UShort_t>(UCHAR_MAX),
                             std::max(static_cast<UShort_t>(0),
                                      t->GetTPCsignalN()));
                _branch_track_its_ncluster[_branch_ntrack] =
                    t->GetNumberOfITSClusters();
                _branch_track_tpc_ncluster_findable[_branch_ntrack] =
                    std::min(static_cast<UShort_t>(UCHAR_MAX),
                             std::max(static_cast<UShort_t>(0),
                                      t->GetTPCNclsF()));

                Double_t dz[2] = { NAN, NAN };
                Double_t cov[3] = { NAN, NAN, NAN };

                if (t->PropagateToDCA
                    (primary_vertex, event->GetMagneticField(),
                     kVeryBig, dz, cov) == kTRUE) {
                    _branch_track_dca_xy[_branch_ntrack] =
                        half(dz[0]);
                    _branch_track_dca_z[_branch_ntrack] =
                        half(dz[1]);
                }

                const Int_t mc_truth_index = t->GetLabel();

#define SAFE_MC_TRUTH_INDEX_TO_USHRT(mc_truth_index)            \
                !(mc_truth_index >= 0 &&                        \
                  static_cast<size_t>(mc_truth_index) <         \
                  stored_mc_truth_index.size()) ? USHRT_MAX :   \
                    stored_mc_truth_index[mc_truth_index] ==    \
                    ULONG_MAX ?                                 \
                    USHRT_MAX :                                 \
                    std::min(static_cast<size_t>(USHRT_MAX),    \
                             std::max(static_cast<size_t>(0),   \
                                      stored_mc_truth_index     \
                                      [mc_truth_index]));

                _branch_track_mc_truth_index[_branch_ntrack] =
                    SAFE_MC_TRUTH_INDEX_TO_USHRT(mc_truth_index);
                _branch_track_voronoi_area[_branch_ntrack] = 0;

                _branch_ntrack++;
                if (_branch_ntrack >= NTRACK_MAX) {
                    break;
                }
            }
        }
    }
    // FIXME: AOD not handled

    std::vector<fastjet::PseudoJet> particle_reco_ue_estimation;
    std::vector<point_2d_t> particle_reco_area_estimation;

    static const double boundary_condition_max_rapidity = 0.9;

    for (std::vector<fastjet::PseudoJet>::const_iterator iterator =
             particle_reco.begin();
         iterator != particle_reco.end(); iterator++) {
        if (fabs(iterator->rapidity()) <
            boundary_condition_max_rapidity) {
            for (int i = -1; i <= 1; i++) {
                // Knuth or euclidean modulus
                const int parity = i & 1;
                const double offset = 2 * i *
                    (boundary_condition_max_rapidity - 0.05);
                fastjet::PseudoJet replica = *iterator;

                replica.reset_momentum_PtYPhiM(
                    iterator->perp(),
                    offset + (2 * parity - 1) * iterator->rapidity(),
                    iterator->phi_std(),
                    iterator->m());
                particle_reco_ue_estimation.push_back(replica);
            }
        }
        particle_reco_area_estimation.push_back(
            point_2d_t(iterator->pseudorapidity(),
                       iterator->phi_std()));
    }

    std::vector<double> particle_reco_area;
    std::vector<std::set<size_t> > particle_reco_incident;

    voronoi_area_incident(particle_reco_area,
                          particle_reco_incident,
                          particle_reco_area_estimation);

    for (size_t i = 0; i < particle_reco_area.size(); i++) {
        if (i < reco_track_index.size() &&
            reco_track_index[i] < _branch_ntrack) {
            _branch_track_voronoi_area[reco_track_index[i]] =
                half(particle_reco_area[i]);
        }
    }

    // Shared by the isolation and jet code

    enum {
        BEAM_PARTICLE_P = 1001
    };

    const bool subtract_ue = esd_event != NULL ?
        !(esd_event->GetBeamParticle(0) == BEAM_PARTICLE_P &&
          esd_event->GetBeamParticle(1) == BEAM_PARTICLE_P) :
        false;

    static const double jet_antikt_d = 0.4;

    std::vector<fastjet::PseudoJet> particle_truth;
    std::vector<fastjet::PseudoJet> jet_truth;
    fastjet::ClusterSequenceArea *cluster_sequence_truth = NULL;

    _branch_nmc_truth = 0;

    // PDG Monte Carlo number scheme

    enum {
        PDG_CODE_ELECTRON_MINUS             = 11,
        PDG_CODE_ELECTRON_NEUTRINO          = 12,
        PDG_CODE_MUON_MINUS                 = 13,
        PDG_CODE_MUON_NEUTRINO              = 14,
        PDG_CODE_TAU_MINUS                  = 15,
        PDG_CODE_TAU_NEUTRINO               = 16,
        PDG_CODE_PHOTON                     = 22,
    };

    std::fill(_branch_met_truth, _branch_met_truth + 2, 0);

    std::vector<fastjet::PseudoJet> particle_reco_tagged =
        particle_reco;

    // A value of 2^(-30) < 10^(-9) would map a 10 TeV particle to
    // less than 10 MeV, sufficient to remove any significant momentum
    // bias while not being too greedy with the dynamic range. Ghost
    // scaling factor is chosen as power of two to maintain the
    // multiplication/division being numerically exact (provided px,
    // py, pz >= 2^(-1022+30) which os of the order 10^(-290) eV).

    static const double scale_ghost = pow(2.0, -30.0);

    if (mc_truth_event != NULL) {
        for (std::vector<Int_t>::const_iterator iterator =
                 reverse_stored_mc_truth_index.begin();
             iterator != reverse_stored_mc_truth_index.end();
             iterator++) {
            const AliMCParticle *p =
                static_cast<AliMCParticle *>(
                    mc_truth_event->GetTrack(*iterator));

            if (p == NULL) {
                continue;
            }

            if (p->GetGeneratorIndex() != 0 || !subtract_ue) {
                const unsigned int abs_pdg_code = std::abs(p->PdgCode());

                switch (abs_pdg_code) {
                case PDG_CODE_ELECTRON_NEUTRINO:
                case PDG_CODE_MUON_NEUTRINO:
                case PDG_CODE_TAU_NEUTRINO:
                    // Remove all (stable) neutrinos from the truth
                    // jet definition
                    break;
                default:
                    particle_truth.push_back(fastjet::PseudoJet(
                        p->Px(), p->Py(), p->Pz(), p->P()));
                    switch (abs_pdg_code) {
                    case PDG_CODE_ELECTRON_MINUS:
                    case PDG_CODE_PHOTON:
                        particle_truth.back().set_user_index(-2);
                        break;
                    case PDG_CODE_MUON_MINUS:
                        particle_truth.back().set_user_index(-3);
                        break;
                    }
                    _branch_met_truth[0] += p->Px();
                    _branch_met_truth[1] += p->Py();
                }
            }
            _branch_mc_truth_e[_branch_nmc_truth] = half(p->E());
            _branch_mc_truth_pt[_branch_nmc_truth] = half(p->Pt());
            _branch_mc_truth_eta[_branch_nmc_truth] = half(p->Eta());
            _branch_mc_truth_phi[_branch_nmc_truth] =
                half(angular_range_reduce(p->Phi()));
            _branch_mc_truth_pdg_code[_branch_nmc_truth] =
                std::min(SHRT_MAX,
                         std::max(SHRT_MIN, p->PdgCode()));
            _branch_mc_truth_status[_branch_nmc_truth] =
                std::min(static_cast<Int_t>(UCHAR_MAX),
                         std::max(0,
                                  mc_truth_event->
                                  Particle(*iterator)->
                                  GetStatusCode()));
            _branch_mc_truth_generator_index[_branch_nmc_truth] =
                std::min(static_cast<Short_t>(UCHAR_MAX),
                         std::max(static_cast<Short_t>(0),
                                  p->GetGeneratorIndex()));
            _branch_nmc_truth++;
            if (_branch_nmc_truth >= NMC_TRUTH_MAX) {
                break;
            }
        }

        cluster_sequence_truth = new fastjet::ClusterSequenceArea(
            particle_truth,
            fastjet::JetDefinition(fastjet::JetDefinition(
                fastjet::antikt_algorithm, jet_antikt_d)),
            fastjet::VoronoiAreaSpec());
        jet_truth = cluster_sequence_truth->inclusive_jets(0);
    }

    // The reco jets will contain EMCAL clusters as ghosts. The idea
    // is calculate a CMS L4 jet energy correction-like
    // electromagnetic fraction (EMF). Its value then can be used as a
    // parameter for the secondary correction, similar to an reversed
    // ATLAS global sequential (GS) correction (the tracking in ALICE
    // being the larger acceptance detector).

    AliVCaloCells *emcal_cell = event->GetEMCALCells();

    for (Int_t i = 0; i < calo_cluster.GetEntriesFast(); i++) {
        AliVCluster *c =
            static_cast<AliVCluster *>(calo_cluster.At(i));

        Int_t cell_id_max = -1;
        Double_t cell_energy_max = -INFINITY;
        Double_t cell_cross = NAN;

        cell_max_cross(cell_id_max, cell_energy_max, cell_cross,
                       c, emcal_cell);
        if (c->GetNCells() > 1 &&
            1 - cell_energy_max / cell_cross < 0.95 &&
            !cell_masked(c, _emcal_mask)) {
            TLorentzVector p;

            c->GetMomentum(p, _branch_primary_vertex);

            const fastjet::PseudoJet
                pj(p.Px(), p.Py(), p.Pz(), p.P());

            particle_reco_tagged.push_back(pj * scale_ghost);
            particle_reco_tagged.back().set_user_index(USER_INDEX_EM);
        }
    }

    _branch_njet_truth = 0;

    if (mc_truth_event != NULL) {
        for (std::vector<fastjet::PseudoJet>::const_iterator
                 iterator_jet = jet_truth.begin();
             iterator_jet != jet_truth.end(); iterator_jet++) {

            _branch_jet_truth_e[_branch_njet_truth] =
                half(iterator_jet->E());
            _branch_jet_truth_pt[_branch_njet_truth] =
                half(iterator_jet->perp());
            _branch_jet_truth_eta[_branch_njet_truth] =
                half(iterator_jet->pseudorapidity());
            _branch_jet_truth_phi[_branch_njet_truth] =
                half(iterator_jet->phi_std());
            _branch_jet_truth_area[_branch_njet_truth] =
                half(iterator_jet->area());

            const std::vector<fastjet::PseudoJet> constituent =
                cluster_sequence_truth->constituents(*iterator_jet);

            _branch_jet_truth_emf[_branch_njet_truth] =
                half(jet_emf(constituent));
            _branch_jet_truth_multiplicity[_branch_njet_truth] =
                jet_multiplicity(constituent);
            jet_width_sigma_h(
                _branch_jet_truth_width_sigma[_branch_njet_truth],
                *iterator_jet, constituent);
            _branch_jet_truth_ptd[_branch_njet_truth] =
                half(jet_ptd(constituent));

            // Part for reco-generator truth tagging

            for (std::vector<fastjet::PseudoJet>::const_iterator
                     iterator_constituent = constituent.begin();
                 iterator_constituent != constituent.end();
                 iterator_constituent++) {

                particle_reco_tagged.push_back(
                    *iterator_constituent * scale_ghost);
                // Positive user indices are used to tag the truth jet
                particle_reco_tagged.back().set_user_index(
                    iterator_jet - jet_truth.begin());
            }

            _branch_njet_truth++;
            if (_branch_njet_truth >= NJET_MAX) {
                break;
            }
        }
    }
    if (cluster_sequence_truth != NULL) {
        delete cluster_sequence_truth;
    }

    for (size_t i = 0; i < particle_reco.size(); i++) {
        particle_reco[i].set_user_index(static_cast<int>(i));
    }

    const fastjet::ClusterSequenceArea
        cluster_sequence_reco(
            particle_reco,
            fastjet::JetDefinition(fastjet::JetDefinition(
                fastjet::antikt_algorithm, jet_antikt_d)),
            fastjet::VoronoiAreaSpec());
    const std::vector<fastjet::PseudoJet> jet_reco =
        cluster_sequence_reco.inclusive_jets(0);

    const fastjet::ClusterSequenceArea
        cluster_sequence_reco_tagged(
            particle_reco_tagged,
            fastjet::JetDefinition(fastjet::JetDefinition(
                fastjet::antikt_algorithm, jet_antikt_d)),
            fastjet::VoronoiAreaSpec());
    const std::vector<fastjet::PseudoJet> jet_reco_tagged =
        cluster_sequence_reco_tagged.inclusive_jets(0);

    static const double jet_kt_d_ue_estimation = 0.2;
    const fastjet::ClusterSequenceArea
        cluster_sequence_ue_estimation(
            particle_reco_ue_estimation,
            fastjet::JetDefinition(fastjet::JetDefinition(
                fastjet::kt_algorithm, jet_kt_d_ue_estimation)),
            fastjet::VoronoiAreaSpec());
    const std::vector<fastjet::PseudoJet> jet_ue_estimation =
        cluster_sequence_ue_estimation.inclusive_jets(0);

    _branch_debug_njet_ue_estimation = 0;
    for (std::vector<fastjet::PseudoJet>::const_iterator
             iterator_jet = jet_ue_estimation.begin();
        iterator_jet != jet_ue_estimation.end(); iterator_jet++) {
        _branch_debug_jet_ue_estimation_pt_raw
            [_branch_debug_njet_ue_estimation] =
            iterator_jet->perp();
        _branch_debug_jet_ue_estimation_eta_raw
            [_branch_debug_njet_ue_estimation] =
            iterator_jet->pseudorapidity();
        _branch_debug_jet_ue_estimation_phi_raw
            [_branch_debug_njet_ue_estimation] =
            iterator_jet->phi_std();
        _branch_debug_jet_ue_estimation_area_raw
            [_branch_debug_njet_ue_estimation] =
            iterator_jet->area();
        _branch_debug_njet_ue_estimation++;
    }

    std::vector<double> ue_estimate =
        ue_estimation_truncated_mean(jet_ue_estimation);

    _branch_ncluster = 0;
    for (Int_t i = 0; i < calo_cluster.GetEntriesFast(); i++) {
        AliVCluster *c =
            static_cast<AliVCluster *>(calo_cluster.At(i));
        TLorentzVector p;

        c->GetMomentum(p, _branch_primary_vertex);

        _branch_cluster_e[_branch_ncluster] = p.E();
        _branch_cluster_pt[_branch_ncluster] = p.Pt();
        _branch_cluster_eta[_branch_ncluster] = p.Eta();
        _branch_cluster_phi[_branch_ncluster] =
            angular_range_reduce(p.Phi());

        _branch_cluster_m02[_branch_ncluster] = c->GetM02();
        _branch_cluster_m20[_branch_ncluster] = c->GetM20();
        _branch_cluster_tof[_branch_ncluster] = c->GetTOF();
        _branch_cluster_ncell[_branch_ncluster] = c->GetNCells();

        _branch_cluster_nmc_truth[_branch_ncluster] =
            c->GetNLabels();
        std::fill(_branch_cluster_mc_truth_index[_branch_ncluster],
                  _branch_cluster_mc_truth_index[_branch_ncluster] + 32,
                  USHRT_MAX);

        // Needed for the isolation below

        std::set<Int_t> cluster_mc_truth_index;

        for (UInt_t j = 0; j < std::min(32U, c->GetNLabels()); j++) {
            _branch_cluster_mc_truth_index[_branch_ncluster][j] =
                c->GetLabelAt(j);
            cluster_mc_truth_index.insert(c->GetLabelAt(j));
        }

        Int_t cell_id_max = -1;
        Double_t cell_energy_max = -INFINITY;
        Double_t energy_cross = NAN;

        cell_max_cross(cell_id_max, cell_energy_max, energy_cross,
                       c, emcal_cell);
        _branch_cluster_cell_id_max[_branch_ncluster] = cell_id_max;
        _branch_cluster_e_max[_branch_ncluster] = cell_energy_max;
        _branch_cluster_e_cross[_branch_ncluster] = energy_cross;

        if (esd_event != NULL) {
            _branch_cluster_iso_tpc_01[_branch_ncluster] = 0;
            _branch_cluster_iso_tpc_02[_branch_ncluster] = 0;
            _branch_cluster_iso_tpc_03[_branch_ncluster] = 0;
            _branch_cluster_iso_tpc_04[_branch_ncluster] = 0;

            std::vector<std::pair<double, double> > delta_vs_iso;

            for (Int_t j = 0; j < esd_event->GetNumberOfTracks(); j++) {
                AliESDtrack *t = esd_event->GetTrack(j);

                if (t == NULL) {
                    continue;
                }

                // Apply PWG-JE cuts (track cuts 0 and 1)

                if (_track_cut[0].AcceptTrack(t) ||
                    _track_cut[1].AcceptTrack(t)) {
                    const double dpseudorapidity = t->Eta() - p.Eta();
                    const double dazimuth = angular_range_reduce(
                        angular_range_reduce(t->Phi()) -
                        angular_range_reduce(p.Phi()));
                    const double dr_2 =
                        std::pow(dpseudorapidity, 2) +
                        std::pow(dazimuth, 2);
                    if (dr_2 < 0.1 * 0.1) {
                        _branch_cluster_iso_tpc_01
                            [_branch_ncluster] += t->Pt();
                    }
                    if (dr_2 < 0.2 * 0.2) {
                        _branch_cluster_iso_tpc_02
                            [_branch_ncluster] += t->Pt();
                    }
                    if (dr_2 < 0.3 * 0.3) {
                        _branch_cluster_iso_tpc_03
                            [_branch_ncluster] += t->Pt();
                    }
                    if (dr_2 < 0.4 * 0.4) {
                        _branch_cluster_iso_tpc_04
                            [_branch_ncluster] += t->Pt();
                        delta_vs_iso.push_back(
                            std::pair<double, double>(
                                sqrt(dr_2), t->Pt()));
                    }
                }
            }
            _branch_cluster_iso_tpc_01[_branch_ncluster] -=
                evaluate_ue(ue_estimate, p.Phi(), 0.1) *
                M_PI * std::pow(0.1, 2);
            _branch_cluster_iso_tpc_02[_branch_ncluster] -=
                evaluate_ue(ue_estimate, p.Phi(), 0.2) *
                M_PI * std::pow(0.2, 2);
            _branch_cluster_iso_tpc_03[_branch_ncluster] -=
                evaluate_ue(ue_estimate, p.Phi(), 0.3) *
                M_PI * std::pow(0.3, 2);
            _branch_cluster_iso_tpc_04[_branch_ncluster] -=
                evaluate_ue(ue_estimate, p.Phi(), 0.4) *
                M_PI * std::pow(0.4, 2);

            _branch_cluster_frixione_tpc_04_02[_branch_ncluster] =
                frixione_iso_max_x_e_eps(delta_vs_iso, 0.4, 0.2,
                                         ue_estimate, p.Phi());
            _branch_cluster_frixione_tpc_04_05[_branch_ncluster] =
                frixione_iso_max_x_e_eps(delta_vs_iso, 0.4, 0.5,
                                         ue_estimate, p.Phi());
            _branch_cluster_frixione_tpc_04_10[_branch_ncluster] =
                frixione_iso_max_x_e_eps(delta_vs_iso, 0.4, 1.0,
                                         ue_estimate, p.Phi());
        }
        else {
            _branch_cluster_iso_tpc_01[_branch_ncluster] = NAN;
            _branch_cluster_iso_tpc_02[_branch_ncluster] = NAN;
            _branch_cluster_iso_tpc_03[_branch_ncluster] = NAN;
            _branch_cluster_iso_tpc_04[_branch_ncluster] = NAN;
            _branch_cluster_frixione_tpc_04_02[_branch_ncluster] = NAN;
            _branch_cluster_frixione_tpc_04_05[_branch_ncluster] = NAN;
            _branch_cluster_frixione_tpc_04_10[_branch_ncluster] = NAN;
        }

        if (mc_truth_event != NULL) {
            _branch_cluster_iso_01_truth[_branch_ncluster] = 0;
            _branch_cluster_iso_02_truth[_branch_ncluster] = 0;
            _branch_cluster_iso_03_truth[_branch_ncluster] = 0;
            _branch_cluster_iso_04_truth[_branch_ncluster] = 0;
            _branch_cluster_frixione_04_02_truth[_branch_ncluster] = 0;
            _branch_cluster_frixione_04_05_truth[_branch_ncluster] = 0;
            _branch_cluster_frixione_04_10_truth[_branch_ncluster] = 0;

            std::vector<std::pair<double, double> > delta_vs_iso;

            for (Int_t j = 0;
                 j < mc_truth_event->GetNumberOfTracks(); j++) {
                if (mc_truth_event->IsPhysicalPrimary(j) &&
                    cluster_mc_truth_index.find(j) !=
                    cluster_mc_truth_index.end()) {
                    const AliMCParticle *t =
                        static_cast<AliMCParticle *>(
                            mc_truth_event->GetTrack(j));

                    if (t->GetGeneratorIndex() != 0 || !subtract_ue) {
                        const double dpseudorapidity = t->Eta() - p.Eta();
                        const double dazimuth = angular_range_reduce(
                            angular_range_reduce(t->Phi()) -
                            angular_range_reduce(p.Phi()));
                        const double dr_2 =
                            std::pow(dpseudorapidity, 2) +
                            std::pow(dazimuth, 2);
                        if (dr_2 < 0.1 * 0.1) {
                            _branch_cluster_iso_01_truth
                                [_branch_ncluster] += t->Pt();
                        }
                        if (dr_2 < 0.2 * 0.2) {
                            _branch_cluster_iso_02_truth
                                [_branch_ncluster] += t->Pt();
                        }
                        if (dr_2 < 0.3 * 0.3) {
                            _branch_cluster_iso_03_truth
                                [_branch_ncluster] += t->Pt();
                        }
                        if (dr_2 < 0.4 * 0.4) {
                            _branch_cluster_iso_04_truth
                                [_branch_ncluster] += t->Pt();
                            delta_vs_iso.push_back(
                                std::pair<double, double>(
                                    sqrt(dr_2), t->Pt()));
                        }
                    }
                }
            }
            _branch_cluster_frixione_04_02_truth[_branch_ncluster] =
                frixione_iso_max_x_e_eps(delta_vs_iso, 0.4, 0.2);
            _branch_cluster_frixione_04_05_truth[_branch_ncluster] =
                frixione_iso_max_x_e_eps(delta_vs_iso, 0.4, 0.5);
            _branch_cluster_frixione_04_10_truth[_branch_ncluster] =
                frixione_iso_max_x_e_eps(delta_vs_iso, 0.4, 1.0);
        }
        else {
            _branch_cluster_iso_01_truth[_branch_ncluster] = NAN;
            _branch_cluster_iso_02_truth[_branch_ncluster] = NAN;
            _branch_cluster_iso_03_truth[_branch_ncluster] = NAN;
            _branch_cluster_iso_04_truth[_branch_ncluster] = NAN;
            _branch_cluster_frixione_04_02_truth[_branch_ncluster] = NAN;
            _branch_cluster_frixione_04_05_truth[_branch_ncluster] = NAN;
            _branch_cluster_frixione_04_10_truth[_branch_ncluster] = NAN;
        }

        _branch_ncluster++;
        if (_branch_ncluster >= NCLUSTER_MAX) {
            break;
        }
    }

    calo_cluster.Delete();

    _branch_njet = 0;
    for (std::vector<fastjet::PseudoJet>::const_iterator
             iterator_jet = jet_reco.begin();
        iterator_jet != jet_reco.end(); iterator_jet++) {

        std::vector<fastjet::PseudoJet>::const_iterator
            iterator_jet_tagged = jet_reco_tagged.end();
        double dr_2_min = INFINITY;

        for (std::vector<fastjet::PseudoJet>::const_iterator
             it = jet_reco_tagged.begin();
             it != jet_reco_tagged.end(); it++) {
            const double dr_2 =
                iterator_jet->squared_distance(*it);

            if (dr_2 < dr_2_min) {
                iterator_jet_tagged = it;
                dr_2_min = dr_2;
            }
        }

        _branch_debug_jet_tag_dr_square[_branch_njet] = dr_2_min;

#if 0
        // Skip jets that only consists of tagging ghosts
        size_t ghost_only = true;

        for (std::vector<fastjet::PseudoJet>::const_iterator
                 iterator_constituent = constituent.begin();
             iterator_constituent != constituent.end();
             iterator_constituent++) {
            switch (iterator_constituent->user_index()) {
            case USER_INDEX_DEFAULT_OR_TRACK:
            case USER_INDEX_EM:
                ghost_only = false;
                break;
            }
        }
        if (ghost_only) {
            continue;
        }
#endif

        if (!(iterator_jet->perp() >= _stored_jet_min_pt_raw)) {
            continue;
        }

        // Jet quantities follow HEP convention (not ALICE so far):
        //
        // - Suffix _raw = raw, jet-uncalibrated detector quantity
        //
        // - Suffix _charged = calibrated, "charged particle-level"
        //   quantity
        //
        // - No suffix = jet-calibrated, particle-level quantity

        _branch_jet_e_raw[_branch_njet] = half(iterator_jet->E());
        _branch_jet_e[_branch_njet] = NAN;

        std::vector<fastjet::PseudoJet> constituent =
            cluster_sequence_reco.constituents(*iterator_jet);
        double area = 0;

        for (std::vector<fastjet::PseudoJet>::const_iterator
                 iterator_constituent = constituent.begin();
             iterator_constituent != constituent.end();
             iterator_constituent++) {
            const int index = iterator_constituent->user_index();

            if (index >= 0 && static_cast<size_t>(index) <
                particle_reco_area.size()) {
                area += particle_reco_area[index];
            }
        }

        const double pt_raw_ue =
            evaluate_ue(ue_estimate, iterator_jet->phi_std(),
                        jet_antikt_d) *
            area;

        _branch_jet_pt_raw_ue[_branch_njet] = half(pt_raw_ue);
        _branch_jet_pt_raw[_branch_njet] =
            half(iterator_jet->perp() - pt_raw_ue);
        _branch_jet_pt[_branch_njet] = NAN;
        _branch_jet_e_charged[_branch_njet] = NAN;
        _branch_jet_pt_charged[_branch_njet] = NAN;
        _branch_jet_eta_raw[_branch_njet] =
            half(iterator_jet->pseudorapidity());
        _branch_jet_eta[_branch_njet] =
            half(iterator_jet->pseudorapidity());
        _branch_jet_phi[_branch_njet] =
            half(iterator_jet->phi_std());
        _branch_jet_area_raw[_branch_njet] = half(area);
        _branch_jet_area[_branch_njet] = half(area);

        // Calculate the electro magnetic fraction (EMF), but without
        // a particle-flow-based removal of energy double counting.
        // Note the EM ghosts are scaled back here.

        _branch_jet_emf_raw[_branch_njet] = NAN;
        _branch_jet_emf[_branch_njet] = NAN;
        _branch_jet_multiplicity_raw[_branch_njet] = 0;
        _branch_jet_multiplicity[_branch_njet] = NAN;
        std::fill(_branch_jet_width_sigma_raw[_branch_njet],
                  _branch_jet_width_sigma_raw[_branch_njet] + 2, NAN);
        std::fill(_branch_jet_width_sigma[_branch_njet],
                  _branch_jet_width_sigma[_branch_njet] + 2, NAN);
        _branch_jet_ptd_raw[_branch_njet] = NAN;
        _branch_jet_ptd[_branch_njet] = NAN;

        if (iterator_jet_tagged != jet_reco_tagged.end()) {
            constituent = cluster_sequence_reco_tagged.
                constituents(*iterator_jet_tagged);

            _branch_jet_emf_raw[_branch_njet] =
                jet_emf(constituent, scale_ghost);
            _branch_jet_multiplicity_raw[_branch_njet] =
                jet_multiplicity(constituent);
            jet_width_sigma_h(
                _branch_jet_width_sigma_raw[_branch_njet],
                *iterator_jet, constituent, scale_ghost);
            _branch_jet_ptd_raw[_branch_njet] =
                jet_ptd(constituent, scale_ghost);
        }

        _branch_jet_e_truth[_branch_njet] = NAN;
        _branch_jet_pt_truth[_branch_njet] = NAN;
        _branch_jet_eta_truth[_branch_njet] = NAN;
        _branch_jet_phi_truth[_branch_njet] = NAN;
        _branch_jet_emf_truth[_branch_njet] = NAN;
        // Defaulting this to USHRT_MAX appears to be too unnatural
        _branch_jet_multiplicity_truth[_branch_njet] = 0;
        std::fill(_branch_jet_width_sigma_truth[_branch_njet],
                  _branch_jet_width_sigma_truth[_branch_njet] + 2, NAN);
        _branch_jet_ptd_truth[_branch_njet] = NAN;

        const size_t index_reco = iterator_jet - jet_reco.begin();

        if (mc_truth_event != NULL &&
            iterator_jet_tagged != jet_reco_tagged.end()) {
            std::map<int, double> z_ghost;

            for (std::vector<fastjet::PseudoJet>::const_iterator
                     iterator_constituent = constituent.begin();
                iterator_constituent != constituent.end();
                iterator_constituent++) {
                const int index_jet_truth =
                    iterator_constituent->user_index();

                if (index_jet_truth >= 0) {
                    if (z_ghost.find(index_jet_truth) ==
                        z_ghost.end()) {
                        z_ghost[index_jet_truth] =
                            iterator_constituent->perp() /
                            scale_ghost;
                    }
                    else {
                        z_ghost[index_jet_truth] +=
                            iterator_constituent->perp() /
                            scale_ghost;
                    }
                }
            }

            // z_truth = fraction of truth constituents inside the
            // area of the reco jet, relative to the truth jet (not
            // necessarily within the reco jet)

            std::vector<std::pair<double, int> > z_truth;

            for (std::map<int, double>::const_iterator iterator =
                     z_ghost.begin();
                 iterator != z_ghost.end(); iterator++) {
                if (jet_truth[iterator->first].perp() > 0) {
                    z_truth.push_back(std::pair<double, int>(
                        iterator->second /
                        jet_truth[iterator->first].perp(),
                        iterator->first));
                }
            }
            std::sort(z_truth.begin(), z_truth.end());

            // Note that z_truth is now in *acending* order. Any
            // information beyond 2 truth -> 1 reco jet mapping is not
            // really useful, we only need to know how fuzzy the
            // mapping was

            std::fill(
                _branch_jet_truth_index_z_truth[_branch_njet],
                _branch_jet_truth_index_z_truth[_branch_njet] + 2,
                -1);
            std::fill(
                _branch_jet_truth_z_truth[_branch_njet],
                _branch_jet_truth_z_truth[_branch_njet] + 2, NAN);
            for (size_t j = 0;
                 j < std::min(2UL, z_truth.size()); j++) {
                _branch_jet_truth_z_truth[_branch_njet][j] =
                    half(z_truth.rbegin()[j].first);
                _branch_jet_truth_index_z_truth[_branch_njet][j] =
                    z_truth.rbegin()[j].second;
            }

            // z_reco = fraction of truth constituents inside the area
            // of the reco jet, relative to the total truth particles
            // inside the reco jet

            double sum_z_ghost = 0;

            for (std::map<int, double>::const_iterator iterator =
                     z_ghost.begin();
                 iterator != z_ghost.end(); iterator++) {
                sum_z_ghost += iterator->second;
            }

            std::vector<std::pair<double, int> > z_reco;

            if (sum_z_ghost > 0) {
                for (std::map<int, double>::iterator iterator =
                         z_ghost.begin();
                     iterator != z_ghost.end(); iterator++) {
                    iterator->second /= sum_z_ghost;
                }
                for (std::map<int, double>::const_iterator iterator =
                         z_ghost.begin();
                     iterator != z_ghost.end(); iterator++) {
                    z_reco.push_back(std::pair<double, int>(
                        iterator->second, iterator->first));
                }
            }
            std::sort(z_reco.begin(), z_reco.end());

            // Note that z_truth is now in *acending* order

            std::fill(
                _branch_jet_truth_index_z_reco[_branch_njet],
                _branch_jet_truth_index_z_reco[_branch_njet] + 2,
                -1);
            std::fill(
                _branch_jet_truth_z_reco[_branch_njet],
                _branch_jet_truth_z_reco[_branch_njet] + 2, NAN);
            for (size_t j = 0;
                 j < std::min(2UL, z_reco.size()); j++) {
                _branch_jet_truth_z_reco[_branch_njet][j] =
                    half(z_reco.rbegin()[j].first);
                _branch_jet_truth_index_z_reco[_branch_njet][j] =
                    z_reco.rbegin()[j].second;
            }

            // A simplified z_reco matching, which is a more rigorous
            // version of the CMS delta R matching, for jet energy
            // correction derivation.

            _branch_jet_e_truth[_branch_njet] = NAN;
            _branch_jet_pt_truth[_branch_njet] = NAN;
            _branch_jet_eta_truth[_branch_njet] = NAN;
            _branch_jet_phi_truth[_branch_njet] = NAN;
            _branch_jet_area_truth[_branch_njet] = NAN;
            _branch_jet_emf_truth[_branch_njet] = NAN;

            if (!z_reco.empty() && z_reco.rbegin()[0].second >= 0 &&
                static_cast<size_t>(z_reco.rbegin()[0].second) <
                _branch_njet_truth) {
                const size_t k = z_reco.rbegin()[0].second;

                _branch_jet_e_truth[_branch_njet] =
                    _branch_jet_truth_e[k];
                _branch_jet_pt_truth[_branch_njet] =
                    _branch_jet_truth_pt[k];
                _branch_jet_eta_truth[_branch_njet] =
                    _branch_jet_truth_eta[k];
                _branch_jet_phi_truth[_branch_njet] =
                    _branch_jet_truth_phi[k];
                _branch_jet_area_truth[_branch_njet] =
                    _branch_jet_truth_area[k];
                _branch_jet_emf_truth[_branch_njet] =
                    _branch_jet_truth_emf[k];
                _branch_jet_multiplicity_truth[_branch_njet] =
                    _branch_jet_truth_multiplicity[k];
                for (size_t j = 0; j < 2; j++) {
                    _branch_jet_width_sigma_truth[_branch_njet][j] =
                        _branch_jet_truth_width_sigma[k][j];
                }
                _branch_jet_ptd_truth[_branch_njet] =
                    _branch_jet_truth_ptd[k];
            }
        }
        _branch_njet++;
        if (_branch_njet >= NJET_MAX) {
            break;
        }
    }

    std::fill(_branch_cell_e, _branch_cell_e + EMCAL_NCELL, NAN);
    std::fill(_branch_cell_tof, _branch_cell_tof + EMCAL_NCELL, NAN);
    std::fill(_branch_cell_mc_truth_index,
              _branch_cell_mc_truth_index + EMCAL_NCELL, USHRT_MAX);
    for (Short_t i = 0; i < emcal_cell->GetNumberOfCells(); i++) {
        Short_t cell_number = -1;
        Double_t cell_energy = NAN;
        Double_t tof = NAN;
        Int_t mc_truth_index = -1;
        Double_t efrac = NAN;

        if (emcal_cell->GetCell(
            i, cell_number, cell_energy, tof, mc_truth_index, efrac) ==
            kTRUE &&
            cell_number >= 0 && cell_number < EMCAL_NCELL) {
            _branch_cell_e[cell_number]   = half(cell_energy);
            _branch_cell_tof[cell_number] = half(tof);
            _branch_cell_mc_truth_index[cell_number] =
                SAFE_MC_TRUTH_INDEX_TO_USHRT(mc_truth_index);

#undef SAFE_MC_TRUTH_INDEX_TO_USHRT

        }
    }

    _branch_nmuon_track = 0;
    if (esd_event != NULL) {
        for (Int_t i = 0; i < esd_event->GetNumberOfMuonTracks(); i++) {
            AliESDMuonTrack *t = esd_event->GetMuonTrack(i);

            if (t == NULL) {
                continue;
            }

            _branch_muon_track_e[_branch_nmuon_track] =
                half(t->E());
            _branch_muon_track_pt[_branch_nmuon_track] =
                half(t->Pt());
            _branch_muon_track_eta[_branch_nmuon_track] =
                half(t->Eta());
            _branch_muon_track_phi[_branch_nmuon_track] =
                half(angular_range_reduce(t->Phi()));

            _branch_muon_track_r_abs[_branch_nmuon_track] =
                half(t->GetRAtAbsorberEnd());
            _branch_muon_track_p_dca[_branch_nmuon_track] = NAN;

            static const double c_505_tan_3_pi_180 =
                505 * tan(3 * M_PI / 180);

            // Default values from
            // AliPhysics/PWG/muon/buildMuonTrackCutsOADB.C

            static const double default_sigma_p_dca_23 = 80;
            static const double default_sigma_p_dca_310 = 54;

            _branch_muon_track_sigma_p_dca[_branch_nmuon_track] =
                t->GetRAtAbsorberEnd() < c_505_tan_3_pi_180 ?
                default_sigma_p_dca_23 : default_sigma_p_dca_310;
            if (_muon_track_cut != NULL) {
                _branch_muon_track_p_dca[_branch_nmuon_track] =
                    half(_muon_track_cut->GetAverageMomentum(t) *
                         _muon_track_cut->GetCorrectedDCA(t).Mag());

                // Here, muon_track_sigma_p_dca is the resoultion
                // corrected value, see
                // AliMuonTrackCuts::GetSelectionMask() in
                // AliPhysics/PWG/muon/AliMuonTrackCuts.cxx

                const AliOADBMuonTrackCutsParam param =
                    _muon_track_cut->GetMuonTrackCutsParam();
                const Double_t sigma_p_dca =
                    _muon_track_cut->IsThetaAbs23(t) ?
                    param.GetSigmaPdca23() : param.GetSigmaPdca310();

                _branch_muon_track_sigma_p_dca[_branch_nmuon_track] =
                    half(sigma_p_dca);

                // In AliMuonTrackCuts::GetSelectionMask(), it would
                // be nrp = nsigma * p * dp

                const Double_t delta_sagitta_p =
                    param.GetRelPResolution() * t->P();

                _branch_muon_track_delta_sagitta_p
                    [_branch_nmuon_track] = half(delta_sagitta_p);

                // p_resolution_effect = sigma_p_dca / (1 - nrp / (1 +
                // nrp));

                static const Double_t z_tc12_cm = 535;
                const Double_t distance_sigma_slope_p_meas =
                    z_tc12_cm * param.GetSlopeResolution() * t->P();

                _branch_muon_track_distance_sigma_slope_p
                    [_branch_nmuon_track] =
                    half(distance_sigma_slope_p_meas);
            }
            _branch_nmuon_track++;
            if (_branch_nmuon_track >= NTRACK_MAX) {
                break;
            }
        }
    }
    else if (aod_event != NULL) {
        // FIXME: Not really implemented
        TRefArray muon_track;

        aod_event->GetMuonTracks(&muon_track);
    }

    _tree_event->Fill();
}

AliEMCALRecoUtils *AliAnalysisTaskNTGJ::GetEMCALRecoUtils(void)
{
    return _reco_util;
}

void AliAnalysisTaskNTGJ::SetAliROOTVersion(const char *version)
{
    strncpy(_branch_version_aliroot, version, BUFSIZ);
    _branch_version_aliroot[BUFSIZ - 1] = '\0';
}

void AliAnalysisTaskNTGJ::SetAliPhysicsVersion(const char *version)
{
    strncpy(_branch_version_aliphysics, version, BUFSIZ);
    _branch_version_aliphysics[BUFSIZ - 1] = '\0';
}

void AliAnalysisTaskNTGJ::SetGridDataDir(const char *dir)
{
    strncpy(_branch_grid_data_dir, dir, BUFSIZ);
    _branch_grid_data_dir[BUFSIZ - 1] = '\0';
}

void AliAnalysisTaskNTGJ::SetGridDataPattern(const char *pattern)
{
    strncpy(_branch_grid_data_pattern, pattern, BUFSIZ);
    _branch_grid_data_pattern[BUFSIZ - 1] = '\0';
}

void AliAnalysisTaskNTGJ::SetSkimClusterMinE(double min_e)
{
    _skim_cluster_min_e = min_e;
}

void AliAnalysisTaskNTGJ::SetSkimTrackMinPt(double min_pt)
{
    _skim_track_min_pt = min_pt;
}

void AliAnalysisTaskNTGJ::SetSkimMuonTrackMinPt(double min_pt)
{
    _skim_muon_track_min_pt = min_pt;
}

void AliAnalysisTaskNTGJ::SetSkimJetMinPt(double min_pt_1,
                                          double min_pt_2,
                                          double min_pt_3)
{
    _skim_jet_min_pt.clear();

    if (min_pt_1 != -INFINITY) {
        _skim_jet_min_pt.push_back(min_pt_1);
        if (min_pt_2 != -INFINITY) {
            _skim_jet_min_pt.push_back(min_pt_2);
            if (min_pt_3 != -INFINITY) {
                _skim_jet_min_pt.push_back(min_pt_3);
            }
        }
    }
}

void AliAnalysisTaskNTGJ::SetSkimMultiplicityTrackletMinN(int min_n)
{
    _skim_multiplicity_tracklet_min_n = min_n;
}

void AliAnalysisTaskNTGJ::SetStoredTrackMinPt(double min_pt)
{
    _stored_track_min_pt = min_pt;
}

// This is primarily intended for JEC derivation, and therefore
// defined in pT raw

void AliAnalysisTaskNTGJ::SetStoredJetMinPtRaw(double min_pt_raw)
{
    _stored_jet_min_pt_raw = min_pt_raw;
}
