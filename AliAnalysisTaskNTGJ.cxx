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
#endif // __CINT__

namespace {

    double angular_range_reduce(const double x)
    {
        if (!std::isfinite(x)) {
            return x;
        }

        static const double cody_waite_x_max = 1608.4954386379741381;
        static const double two_pi_0 = 6.2831853071795649157;
        static const double two_pi_1 = 2.1561211432631314669e-14;
        static const double two_pi_2 = 1.1615423895917441336e-27;
        double ret;

        if(x >= -cody_waite_x_max && x <= cody_waite_x_max) {
            static const double inverse_two_pi =
                0.15915494309189534197;
            const double k = rint(x * inverse_two_pi);
            ret = ((x - (k * two_pi_0)) - k * two_pi_1) -
                k * two_pi_2;
        }
        else {
            long double sin_x;
            long double cos_x;

            sincosl(x, &sin_x, &cos_x);
            ret = (double)atan2l(sin_x, cos_x);
        }
        if(ret == -M_PI) {
            ret = M_PI;
        }

        return ret;
    }

    void to_sm_ieta_iphi(unsigned int &sm, unsigned int &ieta,
                         unsigned int &iphi, unsigned int n)
    {
        sm = n < 11520 ? n / 1152 :
            n < 12288 ? 10 + (n - 11520) / 384 :
            n < 16896 ? 12 + (n - 12288) / 768 :
            18 + (n - 16896) / 384;

        const unsigned int n0 =
            sm < 10 ? sm * 1152 :
            sm < 12 ? 11520 + (sm - 10) * 384 :
            sm < 18 ? 12288 + (sm - 12) * 768 :
            16896 + (sm - 18) * 384;
        const unsigned int n1 = n - n0;
        const unsigned int nphi =
            sm < 10 ? 24 : sm < 12 ? 8 : sm < 18 ? 24 : 8;

        ieta = 2 * (n1 / (2 * nphi)) + 1 - (n1 % 2);
        iphi = (n1 / 2) % nphi;
    }

    void neta_nphi(unsigned int &neta, unsigned int &nphi,
                    const unsigned int sm)
    {
        neta = sm < 12 ? 48 : sm < 18 ? 32 : 48;
        nphi = sm < 10 ? 24 : sm < 12 ? 8 : sm < 18 ? 24 : 8;
    }

    unsigned int flat_sm_ieta(const unsigned int sm,
                              const unsigned int ieta)
    {
        return sm < 12 ? sm * 48 + ieta :
            sm < 18 ? 576 + (sm - 12) * 32 + ieta :
            768 + (sm - 18) * 48 + ieta;
    }

    unsigned int flat_sm_iphi(const unsigned int sm,
                              const unsigned int iphi)
    {
        return sm < 10 ? sm * 24 + iphi :
            sm < 12 ? 240 + (sm - 10) * 8 + iphi :
            sm < 18 ? 256 + (sm - 12) * 24 + iphi :
            400 + (sm - 18) * 8 + iphi;
    }

    bool inside_edge(unsigned int n, unsigned int d)
    {
        unsigned int sm;
        unsigned int ieta;
        unsigned int iphi;

        to_sm_ieta_iphi(sm, ieta, iphi, n);

        unsigned int neta;
        unsigned int nphi;

        neta_nphi(neta, nphi, sm);

        return (ieta >= d && iphi >= d &&
                ieta < neta - d && iphi < nphi - d);
    }

    void cell_3_3(unsigned int n_3_3[], const unsigned int n,
                  const unsigned int ld = 3)
    {
        const unsigned int sm = n < 11520 ? n / 1152 :
            n < 12288 ? 10 + (n - 11520) / 384 :
            n < 16896 ? 12 + (n - 12288) / 768 :
            18 + (n - 16896) / 384;
        const unsigned int nphi =
            sm < 10 ? 24 : sm < 12 ? 8 : sm < 18 ? 24 : 8;

        if (n % 2 == 0) {
            n_3_3[0 * ld + 0] = n - 1;
            n_3_3[0 * ld + 1] = n + 1;
            n_3_3[0 * ld + 2] = n + 3;
        }
        else {
            n_3_3[0 * ld + 0] = n - 2 * nphi - 3;
            n_3_3[0 * ld + 1] = n - 2 * nphi - 1;
            n_3_3[0 * ld + 2] = n - 2 * nphi + 1;
        }
        n_3_3[1 * ld + 0] = n - 2;
        n_3_3[1 * ld + 1] = n;
        n_3_3[1 * ld + 2] = n + 2;
        if (n % 2 == 0) {
            n_3_3[2 * ld + 0] = n + 2 * nphi - 1;
            n_3_3[2 * ld + 1] = n + 2 * nphi + 1;
            n_3_3[2 * ld + 2] = n + 2 * nphi + 3;
        }
        else {
            n_3_3[2 * ld + 0] = n - 3;
            n_3_3[2 * ld + 1] = n - 1;
            n_3_3[2 * ld + 2] = n + 1;
        }
    }

    bool cell_masked(AliVCluster *c, std::vector<bool> emcal_mask)
    {
        for (Int_t j = 0; j < c->GetNCells(); j++) {
            const Int_t cell_id = c->GetCellsAbsId()[j];

            if (cell_id >= 0 &&
                cell_id < static_cast<Int_t>(emcal_mask.size()) &&
                !emcal_mask[cell_id]) {
                return true;
            }
        }

        return false;
    }

}

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
#define BRANCH_STR_ARRAY(b, d)					\
	_branch_ ## b("TObjArray", (d)),

// FIXME: I need to create an interface for:
// _cluster_trigger_min_e(6),
// _jet_min_pt_raw(6),

#define CLASS_INITIALIZATION                    \
    _emcal_geometry_name(EMCAL_GEOMETRY_NAME),  \
    _tree_event(NULL),                          \
    MEMBER_BRANCH                               \
    _f1_ncluster_tpc_linear_pt_dep(NULL),       \
    _track_cut(std::vector<AliESDtrackCuts>()), \
    _reco_util(new AliEMCALRecoUtils),          \
    _emcal_geometry(NULL),                      \
    _ncell(EMCAL_NCELL),                        \
    _cluster_trigger_min_e(-INFINITY),          \
    _jet_min_pt_raw(-INFINITY),                 \
    _emcal_mask(std::vector<bool>()),           \
    _alien_plugin(NULL),                        \
    _metadata_filled(false),                    \
    _prng()

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
#define BRANCH_STR_ARRAY(b, d)					\
    _tree_event->Branch(                        \
        #b, &_branch_ ## b, #b "[" #d "]/C");

    MEMBER_BRANCH;

#undef BRANCH
#undef BRANCH_ARRAY
#undef BRANCH_ARRAY2
#undef BRANCH_STR
#undef BRANCH_STR_ARRAY

    /////////////////////////////////////////////////////////////////

    PostData(1, _tree_event);
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

    AliESDEvent *esd_event = dynamic_cast<AliESDEvent *>(event);

	if (esd_event == NULL) {
		return;
	}

	alice_jec_t jec;

	if (!_metadata_filled) {
		strncpy(_branch_id_git, "$Id$", BUFSIZ);
		_branch_id_git[BUFSIZ - 1] = '\0';
		strncpy(_branch_version_jec, jec.version(), BUFSIZ);
		_branch_version_jec[BUFSIZ - 1] = '\0';
		for (size_t i = 0; i < 2; i++) {
			_branch_beam_particle[i] = esd_event->GetBeamParticle(i);
		}

		const AliESDRun *esd_run = esd_event->GetESDRun();

		if (esd_run != NULL) {
			_branch_trigger_class.Delete();
			for (size_t i = 0; i < NTRIGGER_CLASS_MAX; i++) {
				new (_branch_trigger_class[i]) TObjString(esd_run->GetTriggerClass(i));
			}
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
		for (size_t i = 0; i < 2; i++) {
			_branch_beam_particle[i] = 0;
		}
		_branch_trigger_class.Delete();
	}

    _branch_run_number = event->GetRunNumber();
	_branch_trigger_mask[0] = esd_event->GetTriggerMask();
	_branch_trigger_mask[1] = esd_event->GetTriggerMaskNext50();
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

    for (size_t i = 0; i < 9; i++) {
        _branch_centrality[i] = NAN;
    }
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
    }

    AliMCEvent *mc_truth_event = MCEvent();

    if (mc_truth_event != NULL) {
        mc_truth_event->PreReadAll();
    }

    _branch_eg_ntrial = -1;
    _branch_eg_perp_hat = NAN;
    _branch_eg_cross_section = NAN;

    AliGenEventHeader *mc_truth_header = mc_truth_event != NULL ?
        mc_truth_event->GenEventHeader() : NULL;
    AliGenPythiaEventHeader *mc_truth_pythia_header;

    if (mc_truth_header != NULL) {
        mc_truth_pythia_header =
            dynamic_cast<AliGenPythiaEventHeader *>(mc_truth_header);
        if (mc_truth_pythia_header != NULL) {
            _branch_eg_ntrial = mc_truth_pythia_header->Trials();
            _branch_eg_perp_hat =
                mc_truth_pythia_header->GetPtHard();
            _branch_eg_cross_section =
                mc_truth_pythia_header->GetXsection();
        }
    }

    AliStack *stack;

    if (mc_truth_event != NULL) {
        stack = mc_truth_event->Stack();
    }

#if 0
    AliMagF *field = (AliMagF *)TGeoGlobalMagField::Instance()->GetField();
#endif
    esd_event->InitMagneticField();

    const AliVVertex *primary_vertex = event->GetPrimaryVertex();

    if (primary_vertex != NULL) {
        primary_vertex->GetXYZ(_branch_primary_vertex);
    }

    TRefArray calo_cluster;

    event->GetEMCALClusters(&calo_cluster);

    double cluster_e_max = -INFINITY;

    for (Int_t i = 0; i < calo_cluster.GetEntriesFast(); i++) {
        AliVCluster *c =
            static_cast<AliVCluster *>(calo_cluster.At(i));

        if (!(c->GetNCells() > 1) || cell_masked(c, _emcal_mask)) {
            continue;
        }

        TLorentzVector p;

        c->GetMomentum(p, _branch_primary_vertex);
        cluster_e_max = std::max(cluster_e_max, p.E());
    }

    if (!(cluster_e_max >= _cluster_trigger_min_e)) {
        return;
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

    _branch_ntrack = 0;
    for (Int_t i = 0; i < esd_event->GetNumberOfTracks(); i++) {
        AliESDtrack *t = esd_event->GetTrack(i);

        if (t == NULL) {
            continue;
        }

        // Apply PWG-JE cuts (track cuts 0 and 1)

        if (_track_cut[0].AcceptTrack(t) ||
            _track_cut[1].AcceptTrack(t)) {
            particle_reco.push_back(fastjet::PseudoJet(
                t->Px(), t->Py(), t->Pz(), t->P()));
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

            // Shortened track quality bit mask. Here bit 0 and 1 are
            // the PWG-JE's bit 4 and 8. Test for (track_quality[i] &
            // 3 != 0), i being the track index, to get PWG-JE's "272"
            // (= 1 << 4 | 1 << 8) cut. Test for (track_quality[i] & 4
            // == 0) and (track_quality[i] & 8 == 0) for the "2015
            // PbPb" cut with clusterCut = 0 and 1.

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

            _branch_track_tpc_length_active_zone[_branch_ntrack] = NAN;
            if (t->GetInnerParam() != NULL) {
                _branch_track_tpc_length_active_zone[_branch_ntrack] =
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

            if (t->PropagateToDCA(primary_vertex,
                                  event->GetMagneticField(),
                                  kVeryBig, dz, cov) == kTRUE) {
                _branch_track_dca_xy[_branch_ntrack] = half(dz[0]);
                _branch_track_dca_z[_branch_ntrack] = half(dz[1]);
            }

			const Int_t mc_label = t->GetLabel();

#define SAFE_MC_LABEL_TO_USHRT(mc_label)						\
			!(mc_label >= 0 &&									\
			  static_cast<size_t>(mc_label) <					\
			  stored_mc_truth_index.size()) ? USHRT_MAX :		\
				stored_mc_truth_index[mc_label] == ULONG_MAX ?	\
				USHRT_MAX :										\
				std::min(static_cast<size_t>(USHRT_MAX),		\
						 std::max(static_cast<size_t>(0),		\
								  stored_mc_truth_index			\
								  [mc_label]));

			_branch_track_mc_truth_index[_branch_ntrack] =
				SAFE_MC_LABEL_TO_USHRT(mc_label);

            _branch_ntrack++;
        }
    }

    std::vector<fastjet::PseudoJet> particle_truth;
    std::vector<fastjet::PseudoJet> jet_truth;
    fastjet::ClusterSequenceArea *cluster_sequence_truth = NULL;

    _branch_nmc_truth = 0;

    // fastjet::PseudoJet user indices -2 and -3 are used to tag the
    // EM particles/EMCAL clusters and muons. The index -1 is already
    // taken, being the fastjet::PseudoJet default initializer. After
    // the removal of EM and muons, -1 then implicitly means hadronic

    enum {
        USER_INDEX_DEFAULT_OR_TRACK = -1,
        USER_INDEX_EM               = -2,
        USER_INDEX_MUON             = -3
    };

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

            const unsigned int abs_pdg_code = std::abs(p->PdgCode());

            switch (abs_pdg_code) {
            case PDG_CODE_ELECTRON_NEUTRINO:
            case PDG_CODE_MUON_NEUTRINO:
            case PDG_CODE_TAU_NEUTRINO:
                // Remove all (stable) neutrinos from the truth jet
                // definition
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
            _branch_nmc_truth++;
        }

        cluster_sequence_truth = new fastjet::ClusterSequenceArea(
                particle_truth,
                fastjet::JetDefinition(fastjet::JetDefinition(
                    fastjet::antikt_algorithm, 0.4)),
                fastjet::VoronoiAreaSpec());
        jet_truth = cluster_sequence_truth->inclusive_jets(0);
    }

    // A value of 2^(-30) < 10^(-9) would map a 10 TeV particle to
    // less than 10 MeV, sufficient to remove any significant momentum
    // bias while not being too greedy with the dynamic range. Ghost
    // scaling factor is chosen as power of two to maintain the
    // multiplication/division being numerically exact (provided px,
    // py, pz >= 2^(-1022+30) which os of the order 10^(-290) eV).

    static const double scale_ghost = pow(2.0, -30.0);

    // The reco jets will contain EMCAL clusters as ghosts. The idea
    // is calculate a CMS L4 jet energy correction-like
    // electromagnetic fraction (EMF). Its value then can be used as a
    // parameter for the secondary correction, similar to an reversed
    // ATLAS global sequential (GS) correction (the tracking in ALICE
    // being the larger acceptance detector).

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

        _branch_ncluster++;

        // Reuse the loop and also fill the reco pseudojets with EMCAL
        // clusters here (as ghosts)

        if (c->GetNCells() > 1 && !cell_masked(c, _emcal_mask)) {
            const fastjet::PseudoJet
                pj(p.Px(), p.Py(), p.Pz(), p.P());

            particle_reco.push_back(pj * scale_ghost);
            particle_reco.back().set_user_index(USER_INDEX_EM);
        }
    }

	calo_cluster.Delete();

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
            double sum_hadronic = 0;
            double sum_em = 0;

            for (std::vector<fastjet::PseudoJet>::const_iterator
                     iterator_constituent = constituent.begin();
                 iterator_constituent != constituent.end();
                 iterator_constituent++) {
                // Part to calculate the generator truth
                // electromagnetic fraction

                if (iterator_constituent->user_index() ==
                    USER_INDEX_EM) {
                    sum_em += iterator_constituent->perp();
                }
                else if (iterator_constituent->user_index() !=
                         USER_INDEX_MUON) {
                    sum_hadronic += iterator_constituent->perp();
                }

                // Part for reco-generator truth tagging

                particle_reco.push_back(
                    *iterator_constituent * scale_ghost);
                // Positive user indices are used to tag the truth jet
                particle_reco.back().set_user_index(
                    iterator_jet - jet_truth.begin());
            }
            _branch_jet_truth_emf[_branch_njet_truth] =
                sum_em / (sum_hadronic + sum_em);

            _branch_njet_truth++;
        }
    }
    if (cluster_sequence_truth != NULL) {
        delete cluster_sequence_truth;
    }

	enum {
		BEAM_PARTICLE_P = 1001
	};

	const bool jet_subtract_ue = esd_event != NULL &&
		!(esd_event->GetBeamParticle(0) == BEAM_PARTICLE_P &&
		  esd_event->GetBeamParticle(1) == BEAM_PARTICLE_P);

    const fastjet::ClusterSequenceArea
        cluster_sequence_reco(
            particle_reco,
            fastjet::JetDefinition(fastjet::JetDefinition(
                fastjet::antikt_algorithm, 0.4)),
            fastjet::VoronoiAreaSpec());
    const std::vector<fastjet::PseudoJet> jet_reco =
        cluster_sequence_reco.inclusive_jets(0);

    _branch_njet = 0;
    for (std::vector<fastjet::PseudoJet>::const_iterator
             iterator_jet = jet_reco.begin();
        iterator_jet != jet_reco.end(); iterator_jet++) {

        const std::vector<fastjet::PseudoJet> constituent =
            cluster_sequence_reco.constituents(*iterator_jet);

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

        if (!(iterator_jet->perp() >= _jet_min_pt_raw)) {
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
        //
        // FIXME: eta, phi, area currently violates this logic, just
        // like in CMS (ATLAS at least calibrates eta)

        _branch_jet_e_raw[_branch_njet] = half(iterator_jet->E());
        _branch_jet_pt_raw[_branch_njet] = half(iterator_jet->perp());
        _branch_jet_e[_branch_njet] = NAN;
        _branch_jet_pt[_branch_njet] = NAN;
        _branch_jet_e_charged[_branch_njet] = NAN;
        _branch_jet_pt_charged[_branch_njet] = NAN;
        _branch_jet_eta_raw[_branch_njet] =
            half(iterator_jet->pseudorapidity());
        _branch_jet_eta[_branch_njet] =
            half(iterator_jet->pseudorapidity());
        _branch_jet_phi[_branch_njet] =
            half(iterator_jet->phi_std());
        _branch_jet_area[_branch_njet] =
            half(iterator_jet->area());

        // Calculate the electro magnetic fraction (EMF), but without
        // a particle-flow-based removal of energy double counting.

        double sum_track = 0;
        double sum_em = 0;

        for (std::vector<fastjet::PseudoJet>::const_iterator
                 iterator_constituent = constituent.begin();
             iterator_constituent != constituent.end();
             iterator_constituent++) {
            if (iterator_constituent->user_index() ==
                USER_INDEX_EM) {
                sum_em += iterator_constituent->perp();
            }
            else {
                sum_track += iterator_constituent->perp();
            }
        }
        _branch_jet_emf_raw[_branch_njet] =
            sum_em / (sum_track * scale_ghost + sum_em);
        // Future extension for a calibrated EMF
        _branch_jet_emf[_branch_njet] = NAN;

        _branch_jet_e_truth[_branch_njet] = NAN;
        _branch_jet_pt_truth[_branch_njet] = NAN;
        _branch_jet_eta_truth[_branch_njet] = NAN;
        _branch_jet_phi_truth[_branch_njet] = NAN;
        _branch_jet_emf_truth[_branch_njet] = NAN;

        const size_t index_reco = iterator_jet - jet_reco.begin();

        if (mc_truth_event != NULL) {
            const std::vector<fastjet::PseudoJet> constituent =
                cluster_sequence_reco.constituents(*iterator_jet);

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

            for (size_t j = 0; j < 2; j++) {
                _branch_jet_truth_index_z_truth[_branch_njet][j] = -1;
                _branch_jet_truth_z_truth[_branch_njet][j] = NAN;
            }
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

            for (size_t j = 0; j < 2; j++) {
                _branch_jet_truth_index_z_reco[_branch_njet][j] = -1;
                _branch_jet_truth_z_reco[_branch_njet][j] = NAN;
            }
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
            }
        }
        _branch_njet++;
    }

	AliVCaloCells *emcal_cell = event->GetEMCALCells();

	for (size_t i = 0; i < EMCAL_NCELL; i++) {
		_branch_cell_amplitude[i] = NAN;
		_branch_cell_time[i] = NAN;
		_branch_cell_mc_truth_index[i] = USHRT_MAX;
		_branch_cell_efrac[i] = 0;
	}
	for (Short_t i = 0; i < emcal_cell->GetNumberOfCells(); i++) {
		Short_t cell_number = -1;
		Double_t amplitude = NAN;
		Double_t time = NAN;
		Int_t mc_label = -1;
		Double_t efrac = NAN;

		if (emcal_cell->GetCell(
			i, cell_number, amplitude, time, mc_label, efrac) ==
			kTRUE &&
            cell_number >= 0 && cell_number < EMCAL_NCELL) {
            _branch_cell_amplitude[cell_number] = half(amplitude);
            _branch_cell_time[cell_number]      = half(time);


            _branch_cell_mc_truth_index[cell_number] =
				SAFE_MC_LABEL_TO_USHRT(mc_label);

#undef SAFE_MC_LABEL_TO_USHRT

			// efrac is stored as INT8
            _branch_cell_efrac[cell_number] = 
				std::min(static_cast<double>(UCHAR_MAX),
						 std::max(0.0, efrac * UCHAR_MAX));
        }
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
