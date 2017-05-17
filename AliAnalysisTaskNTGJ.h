// -*- mode: c++; -*-

#ifndef ALIANALYSISTASKPHOTONDISC_H_
#define ALIANALYSISTASKPHOTONDISC_H_

#include <limits.h>
#include <vector>
#include <TList.h>
#include <TH2D.h>
#include <TFormula.h>
#include <AliESDtrackCuts.h>
#include <AliEMCALGeometry.h>
#include <AliEMCALRecoUtils.h>
#include <AliAnalysisTaskSE.h>

#define EMCAL_NCELL 17664

#define NCLUSTER_MAX 131072
#define NTRACK_MAX 131072
#define NMC_TRUTH_MAX 131072
#define NJET_MAX 131072

class AliAnalysisTaskNTGJ : public AliAnalysisTaskSE {
private:
	TString _emcal_geometry_name; //!
	TList *_list; //!

	TTree *_tree_event; //!

#define MEMBER_BRANCH												\
	BRANCH(run_number, I)											\
	BRANCH(mixed_event, B)											\
	BRANCH_ARRAY(multiplicity_v0a, 4, F)							\
	BRANCH_ARRAY(multiplicity_v0c, 4, F)							\
	BRANCH_ARRAY(centrality, 41, F)									\
	BRANCH(has_misalignment_matrix, O)								\
	BRANCH(eg_ntrial, I)											\
	BRANCH(eg_perp_hat, F)											\
	BRANCH(eg_cross_section, F)										\
	BRANCH_ARRAY(primary_vertex, 3, D)								\
	BRANCH(ncluster, l)												\
	BRANCH_ARRAY(cluster_e, ncluster, F)							\
	BRANCH_ARRAY(cluster_pt, ncluster, F)							\
	BRANCH_ARRAY(cluster_eta, ncluster, F)							\
	BRANCH_ARRAY(cluster_phi, ncluster, F)							\
	BRANCH_ARRAY(cluster_m02, ncluster, F)							\
	BRANCH_ARRAY(cluster_m20, ncluster, F)							\
	BRANCH_ARRAY(cluster_tof, ncluster, F)							\
	BRANCH_ARRAY(cluster_ncell, ncluster, I)						\
	BRANCH(ntrack, l)												\
	BRANCH_ARRAY(track_e, ntrack, F)								\
	BRANCH_ARRAY(track_pt, ntrack, F)								\
	BRANCH_ARRAY(track_eta, ntrack, F)								\
	BRANCH_ARRAY(track_phi, ntrack, F)								\
	BRANCH(nmc_truth, l)											\
	BRANCH_ARRAY(mc_truth_e, nmc_truth, F)							\
	BRANCH_ARRAY(mc_truth_pt, nmc_truth, F)							\
	BRANCH_ARRAY(mc_truth_eta, nmc_truth, F)						\
	BRANCH_ARRAY(mc_truth_phi, nmc_truth, F)						\
	BRANCH_ARRAY(mc_truth_label, nmc_truth, I)						\
	BRANCH_ARRAY(mc_truth_pdg_id, nmc_truth, I)						\
	BRANCH_ARRAY(mc_truth_status, nmc_truth, I)						\
	BRANCH(njet, l)													\
	BRANCH_ARRAY(jet_e_raw, njet, F)								\
	BRANCH_ARRAY(jet_e, njet, F)									\
	BRANCH_ARRAY(jet_e_charged, njet, F)							\
	BRANCH_ARRAY(jet_pt_raw, njet, F)								\
	BRANCH_ARRAY(jet_pt, njet, F)									\
	BRANCH_ARRAY(jet_pt_charged, njet, F)							\
	BRANCH_ARRAY(jet_eta, njet, F)									\
	BRANCH_ARRAY(jet_phi, njet, F)									\
	BRANCH_ARRAY2(jet_truth_index_z_truth, njet, 2, I)				\
	BRANCH_ARRAY2(jet_truth_z_truth, njet, 2, F)					\
	BRANCH_ARRAY2(jet_truth_index_z_reco, njet, 2, I)				\
	BRANCH_ARRAY2(jet_truth_z_reco, njet, 2, F)						\


#define B Char_t
#define b UChar_t
#define S Short_t
#define s UShort_t
#define I Int_t
#define i UInt_t
#define F Float_t
#define D Double_t
#define L Long_t
#define l ULong_t
#define O Bool_t
#define ncluster NCLUSTER_MAX
#define ntrack NTRACK_MAX
#define nmc_truth NMC_TRUTH_MAX
#define njet NJET_MAX
#define full_emcal_ncell EMCAL_NCELL
#define BRANCH(b, t)							\
	t _branch_ ## b;
#define BRANCH_ARRAY(b, d, t)					\
	t _branch_ ## b [(d)];
#define BRANCH_ARRAY2(b, d, e, t)				\
	t _branch_ ## b [(d)][(e)];

	MEMBER_BRANCH;

#undef BRANCH
#undef BRANCH_ARRAY
#undef BRANCH_ARRAY2
#undef nmc_truth
#undef njet
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

	TFormula *_f1_ncluster_tpc_linear_pt_dep; //!
	std::vector<AliESDtrackCuts> _track_cut; //!

	AliEMCALRecoUtils *_reco_util; //!
	AliEMCALGeometry *_emcal_geometry; //!

	size_t _ncell; //!
	double _cluster_trigger_min_e; //!
	double _jet_min_pt_raw; //!

	std::vector<bool> _emcal_mask;

	TRandom3 _prng; //!
public:
	AliAnalysisTaskNTGJ(void);
	AliAnalysisTaskNTGJ(const char *name);
	AliAnalysisTaskNTGJ(const AliAnalysisTaskNTGJ &);
	AliAnalysisTaskNTGJ &operator=
		(const AliAnalysisTaskNTGJ &);
	~AliAnalysisTaskNTGJ(void);
	virtual void UserCreateOutputObjects(void);
	virtual void UserExec(Option_t *);
	AliEMCALRecoUtils *GetEMCALRecoUtils(void);
	ClassDef(AliAnalysisTaskNTGJ, 1);
};

#endif // ALIANALYSISTASKPHOTONDISC_H_
