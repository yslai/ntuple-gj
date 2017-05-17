// This is purely to suppress warning inside ROOT once instanciating
// std::vector<bool>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#include <TCollectionProxyInfo.h>
#pragma GCC diagnostic pop

#include <TGeoManager.h>
#include <TFile.h>
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

#include <TGeoGlobalMagField.h>
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
		if(!std::isfinite(x)) {
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

}

ClassImp(AliAnalysisTaskNTGJ);

#define EMCAL_GEOMETRY_NAME	"EMCAL_COMPLETE12SMV1_DCAL_8SM"

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
#define BRANCH(b, t)							\
	_branch_ ## b((t)),
#define BRANCH_ARRAY(b, d, t)
#define BRANCH_ARRAY2(b, d, e, t)

#define CLASS_INITIALIZATION						\
	_emcal_geometry_name(EMCAL_GEOMETRY_NAME),		\
		_list(NULL),								\
		_tree_event(NULL),							\
		MEMBER_BRANCH								\
		_f1_ncluster_tpc_linear_pt_dep(NULL),		\
		_track_cut(std::vector<AliESDtrackCuts>()),	\
		_reco_util(new AliEMCALRecoUtils),			\
		_emcal_geometry(NULL),						\
		_ncell(EMCAL_NCELL),						\
		_cluster_trigger_min_e(6),					\
		_jet_min_pt_raw(6),							\
		_emcal_mask(std::vector<bool>()),			\
		_prng()

AliAnalysisTaskNTGJ::AliAnalysisTaskNTGJ(void)
	: AliAnalysisTaskSE(), CLASS_INITIALIZATION
{
}

AliAnalysisTaskNTGJ::AliAnalysisTaskNTGJ(
	const char *name)
	: AliAnalysisTaskSE(name), CLASS_INITIALIZATION
{
	DefineOutput(1, TList::Class());
}

AliAnalysisTaskNTGJ::AliAnalysisTaskNTGJ(
	const AliAnalysisTaskNTGJ &x)
	: AliAnalysisTaskSE(), CLASS_INITIALIZATION
{
}

#undef BRANCH
#undef BRANCH_ARRAY
#undef BRANCH_ARRAY2
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
	if (_list != NULL) {
		delete _list;
	}

	delete _reco_util;
}

void AliAnalysisTaskNTGJ::UserCreateOutputObjects(void)
{
	_list = new TList();
	_list->SetOwner(kTRUE);

	/////////////////////////////////////////////////////////////////

	_tree_event =
		new TTree("_tree_event", "");

#define BRANCH(b, t)					\
	_tree_event->Branch(				\
		#b, &_branch_ ## b, #b "/" #t);
#define BRANCH_ARRAY(b, d, t)					\
	_tree_event->Branch(						\
		#b, _branch_ ## b, #b "[" #d "]/" #t);
#define BRANCH_ARRAY2(b, d, e, t)						\
	_tree_event->Branch(								\
		#b, _branch_ ## b, #b "[" #d "][" #e "]/" #t);

	MEMBER_BRANCH;

#undef BRANCH
#undef BRANCH_ARRAY
#undef BRANCH_ARRAY2

	_list->Add(_tree_event);

	/////////////////////////////////////////////////////////////////

	PostData(1, _list);
}

#undef MEMBER_BRANCH

void AliAnalysisTaskNTGJ::UserExec(Option_t *option)
{
	AliVEvent *event = InputEvent();

	if (event == NULL) {
		return;
	}

    AliESDEvent *esd_event = dynamic_cast<AliESDEvent *>(event);
#if 0
    AliAODEvent *aod_event = dynamic_cast<AliAODEvent *>(event);
#endif

	if (_emcal_mask.size() != _ncell) {
		_emcal_mask.resize(_ncell);
#if 0
		for (unsigned int i = 0; i < _ncell; i++) {
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
		for (unsigned int i = 0; i < _ncell; i++) {
			_emcal_mask[i] = true;
		}
#endif
	}

	_branch_run_number = event->GetRunNumber();
	_branch_has_misalignment_matrix = false;

	if (_f1_ncluster_tpc_linear_pt_dep == NULL) {
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
	}

	static const char *centrality_method[41] = {
		"V0M", "V0A", "V0A0", "V0A123", "V0C", "V0A23", "V0C01",
		"V0S", "V0MEq", "V0AEq", "V0CEq", "FMD", "TRK", "TKL", "CL0",
		"CL1", "CND", "ZNA", "ZNC", "ZPA", "ZPC", "NPA", "V0MvsFMD",
		"TKLvsV0M", "ZEMvsZDC", "V0Mtrue", "V0Atrue", "V0Ctrue",
		"V0MEqtrue", "V0AEqtrue", "V0CEqtrue", "FMDtrue", "TRKtrue",
		"TKLtrue", "CL0true", "CL1true", "CNDtrue", "ZNAtrue",
		"ZNCtrue", "ZPAtrue", "ZPCtrue"
	};

	AliVVZERO *v0 = event->GetVZEROData();

	if (v0 != NULL) {
		for (size_t i = 0; i < 4; i++) {
			_branch_multiplicity_v0a[i] = v0->GetMRingV0A(i);
			_branch_multiplicity_v0c[i] = v0->GetMRingV0C(i);
		}
	}

	AliCentrality *centrality = event->GetCentrality();

	if (centrality != NULL) {
		for (size_t i = 0; i < 41; i++) {
			_branch_centrality[i] =
				centrality->GetCentralityPercentile(
					centrality_method[i]);
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
	event->GetPrimaryVertex()->GetXYZ(_branch_primary_vertex);

	TRefArray *calo_cluster = new TRefArray();

	event->GetEMCALClusters(calo_cluster);

	double cluster_e_max = -INFINITY;

	for (Int_t i = 0; i < calo_cluster->GetEntriesFast(); i++) {
		AliVCluster *c =
			static_cast<AliVCluster *>(calo_cluster->At(i));
		Int_t cell_id_max = -1;
		Double_t amplitude_fraction_max = -INFINITY;

		for (Int_t j = 0; j < c->GetNCells(); j++) {
			const Int_t cell_id = c->GetCellsAbsId()[j];
			const Double_t amplitude_fraction =
				 c->GetCellsAmplitudeFraction()[j];

			if (amplitude_fraction > amplitude_fraction_max) {
				cell_id_max = cell_id;
				amplitude_fraction_max = amplitude_fraction;
			}
		}

		TLorentzVector p;

		c->GetMomentum(p, _branch_primary_vertex);

		if (c->GetNCells() > 1 &&
			cell_id_max >= 0 &&
			cell_id_max < static_cast<Int_t>(_ncell) &&
			_emcal_mask[cell_id_max]) {
			cluster_e_max = std::max(cluster_e_max, p.E());
		}
	}

#if 1
	if (!(cluster_e_max >= _cluster_trigger_min_e)) {
		return;
	}
#endif

	_branch_ncluster = 0;
	for (Int_t i = 0; i < calo_cluster->GetEntriesFast(); i++) {
		AliVCluster *c =
			static_cast<AliVCluster *>(calo_cluster->At(i));
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
	}

	std::vector<fastjet::PseudoJet> particle_reco;

	_branch_ntrack = 0;
	for (Int_t i = 0; i < esd_event->GetNumberOfTracks(); i++) {
        AliESDtrack *t = esd_event->GetTrack(i);

        if (t == NULL) {
			continue;
		}

		if (_track_cut[0].AcceptTrack(t) ||
			_track_cut[1].AcceptTrack(t)) {
			particle_reco.push_back(fastjet::PseudoJet(
				t->Px(), t->Py(), t->Pz(), t->P()));
			_branch_track_e[_branch_ntrack] = half(t->E());
			_branch_track_pt[_branch_ntrack] = half(t->Pt());
			_branch_track_eta[_branch_ntrack] = half(t->Eta());
			_branch_track_phi[_branch_ntrack] =
				half(angular_range_reduce(t->Phi()));
			_branch_ntrack++;
		}
	}

	// fprintf(stdout, "%s:%d: %f %d %d %lu\n", __FILE__, __LINE__, cluster_e_max, event->GetNumberOfTracks(), esd_event->GetNumberOfESDTracks(), _branch_ntrack);

	std::vector<fastjet::PseudoJet> particle_truth;
	std::vector<fastjet::PseudoJet> jet_truth;
	fastjet::ClusterSequenceArea *cluster_sequence_truth = NULL;

	_branch_nmc_truth = 0;

	if (mc_truth_event != NULL) {
		for (Int_t i = 0;
			 i < mc_truth_event->GetNumberOfTracks(); i++) {
			const AliMCParticle *p =
				static_cast<AliMCParticle *>(
					mc_truth_event->GetTrack(i));

			if (p == NULL) {
				continue;
			}

			switch (std::abs(p->PdgCode())) {
			case 12:
			case 14:
			case 16:
			case 18:
				// Remove neutrinos from the truth jet definition
				break;
			default:
				particle_truth.push_back(fastjet::PseudoJet(
					p->Px(), p->Py(), p->Pz(), p->P()));
			}
			_branch_mc_truth_e[_branch_nmc_truth] = half(p->E());
			_branch_mc_truth_pt[_branch_nmc_truth] = half(p->Pt());
			_branch_mc_truth_eta[_branch_nmc_truth] = half(p->Eta());
			_branch_mc_truth_phi[_branch_nmc_truth] =
				half(angular_range_reduce(p->Phi()));
			_branch_nmc_truth++;
		}

		cluster_sequence_truth = new fastjet::ClusterSequenceArea(
				particle_truth,
				fastjet::JetDefinition(fastjet::JetDefinition(
					fastjet::antikt_algorithm, 0.4)),
				fastjet::VoronoiAreaSpec());
		jet_truth = cluster_sequence_truth->inclusive_jets(0);
	}

	// 10^(-8) would map a 10 TeV particle to 100 MeV, sufficient to
	// remove any significant momentum bias

	static const double scale_ghost = 1e-8;

	if (mc_truth_event != NULL) {
		for(std::vector<fastjet::PseudoJet>::const_iterator
				iterator_jet = jet_truth.begin();
			iterator_jet != jet_truth.end(); iterator_jet++) {
			const std::vector<fastjet::PseudoJet> constituent =
				cluster_sequence_truth->constituents(*iterator_jet);
			for(std::vector<fastjet::PseudoJet>::const_iterator
					iterator_constituent = constituent.begin();
				iterator_constituent != constituent.end();
				iterator_constituent++) {
				particle_reco.push_back(
					*iterator_constituent * scale_ghost);
				particle_reco.back().set_user_index(
					iterator_jet - jet_truth.begin());
			}
		}
	}
	if (cluster_sequence_truth != NULL) {
		delete cluster_sequence_truth;
	}

	const fastjet::ClusterSequenceArea
		cluster_sequence_reco(
			particle_reco,
			fastjet::JetDefinition(fastjet::JetDefinition(
				fastjet::antikt_algorithm, 0.4)),
			fastjet::VoronoiAreaSpec());
	const std::vector<fastjet::PseudoJet> jet_reco =
		cluster_sequence_reco.inclusive_jets(0);

#if 1
	_branch_njet = 0;
	for(std::vector<fastjet::PseudoJet>::const_iterator
			iterator_jet = jet_reco.begin();
		iterator_jet != jet_reco.end(); iterator_jet++) {
		if (!(iterator_jet->perp() >= _jet_min_pt_raw)) {
			continue;
		}

		_branch_jet_e_raw[_branch_njet] = half(iterator_jet->E());
		_branch_jet_pt_raw[_branch_njet] = half(iterator_jet->perp());
		_branch_jet_e[_branch_njet] = NAN;
		_branch_jet_pt[_branch_njet] = NAN;
		_branch_jet_e_charged[_branch_njet] = NAN;
		_branch_jet_pt_charged[_branch_njet] = NAN;
		_branch_jet_eta[_branch_njet] =
			half(iterator_jet->pseudorapidity());
		_branch_jet_phi[_branch_njet] =
			half(iterator_jet->phi_std());

		const size_t index_reco = iterator_jet - jet_reco.begin();

#if 1
		if (mc_truth_event != NULL) {
			const std::vector<fastjet::PseudoJet> constituent =
				cluster_sequence_reco.constituents(*iterator_jet);

			std::map<int, double> z_ghost;

			for(std::vector<fastjet::PseudoJet>::const_iterator
					iterator_constituent = constituent.begin();
				iterator_constituent != constituent.end();
				iterator_constituent++) {
				if (iterator_constituent->user_index() >= 0) {
					const int index_jet_truth =
						iterator_constituent->user_index();

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

				// z_truth = fraction of truth constituents inside the
				// area of the reco jet, relative to the truth jet
				// (not necessarily within the reco jet)

				std::vector<std::pair<double, int> > z_truth;

				double sum = 0;
				for (std::map<int, double>::iterator iterator =
						 z_ghost.begin();
					 iterator != z_ghost.end(); iterator++) {
					sum += iterator->second;
				}

				for (std::map<int, double>::iterator iterator =
						 z_ghost.begin();
					 iterator != z_ghost.end(); iterator++) {
					z_truth.push_back(std::pair<double, int>(
						iterator->second /
						jet_truth[iterator->first].perp(),
						iterator->first));
				}
				std::sort(z_truth.begin(), z_truth.end());

				// Note that z_truth is now in *acending* order. Any
				// information beyond 2 truth -> 1 reco jet mapping is
				// not really useful, we only need to know how fuzzy
				// the mapping was

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

				// z_reco = fraction of truth constituents inside the
				// area of the reco jet, relative to the total truth
				// particles inside the reco jet

				std::vector<std::pair<double, int> > z_reco;

				if (sum > 0) {
					for (std::map<int, double>::iterator iterator =
							 z_ghost.begin();
						 iterator != z_ghost.end(); iterator++) {
						iterator->second /= sum;
					}
					for (std::map<int, double>::iterator iterator =
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
			}
		}
#endif
		_branch_njet++;
	}
#endif

	_tree_event->Fill();
}

AliEMCALRecoUtils *AliAnalysisTaskNTGJ::GetEMCALRecoUtils(void)
{
	return _reco_util;
}
