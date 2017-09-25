// -*- mode: c++; -*-

#ifndef ALIANALYSISTASKPHOTONDISC_H_
#define ALIANALYSISTASKPHOTONDISC_H_

#include <limits.h>
#include <vector>
#include <set>
#include <TList.h>
#include <TH2D.h>
#include <TFormula.h>
#include <TClonesArray.h>
#include <AliESDtrackCuts.h>
#include <AliEMCALGeometry.h>
#include <AliEMCALRecoUtils.h>
#include <AliMuonTrackCuts.h>
#include <AliAnalysisTaskSE.h>
#include <AliAnalysisAlien.h>

#define EMCAL_NCELL         17664
#define NTRIGGER_CLASS_MAX  100 // = AliESDRun::kNTriggerClasses

#define NCLUSTER_MAX        (1U << 17)
#define NTRACK_MAX          (1U << 17)
#define NMC_TRUTH_MAX       (1U << 17)
#define NJET_MAX            (1U << 17)

class AliAnalysisTaskNTGJ : public AliAnalysisTaskSE {
private:
    TString _emcal_geometry_name; //!
    TTree *_tree_event; //!

#define MEMBER_BRANCH                                               \
    /* */                                                           \
    BRANCH_STR(id_git)                                              \
    BRANCH_STR(version_aliroot)                                     \
    BRANCH_STR(version_aliphysics)                                  \
    BRANCH_STR(version_jec)                                         \
    BRANCH_STR(grid_data_dir)                                       \
    BRANCH_STR(grid_data_pattern)                                   \
    BRANCH_ARRAY(beam_particle, 2, I)                               \
    BRANCH(ntrigger_class, b)                                       \
    BRANCH_STR_ARRAY(trigger_class, ntrigger_class)                 \
    BRANCH(run_number, I)                                           \
    /* */                                                           \
    BRANCH_ARRAY(trigger_mask, 2, l)                                \
    BRANCH_ARRAY(multiplicity_v0, 64, F)                            \
    BRANCH(centrality_v0m, F)                                       \
    BRANCH_ARRAY(centrality, 9, F)                                  \
    BRANCH_ARRAY(event_plane_psi_v0, 3, F)                          \
    BRANCH_ARRAY2(event_plane_q_v0, 3, 2, D)                        \
    BRANCH(has_misalignment_matrix, O)                              \
    BRANCH_ARRAY(cell_eta, 17664, F)                                \
    BRANCH_ARRAY(cell_phi, 17664, F)                                \
    BRANCH_ARRAY(cell_voronoi_area, 17664, F)                       \
    BRANCH_ARRAY(primary_vertex, 3, D)                              \
    BRANCH_ARRAY(primary_vertex_sigma, 3, D)                        \
    BRANCH(primary_vertex_ncontributor, I)                          \
    BRANCH_ARRAY(primary_vertex_spd, 3, D)                          \
    BRANCH_ARRAY(primary_vertex_spd_sigma, 3, D)                    \
    BRANCH(primary_vertex_spd_ncontributor, I)                      \
    BRANCH(npileup_vertex_spd, I)                                   \
    BRANCH(pileup_vertex_spd_ncontributor, I)                       \
    BRANCH(pileup_vertex_spd_min_z_distance, D)                     \
    BRANCH(eg_signal_process_id, I)                                 \
    BRANCH(eg_mpi, I)                                               \
    BRANCH(eg_pt_hat, F)                                            \
    BRANCH(eg_cross_section, F)                                     \
    BRANCH(eg_weight, F)                                            \
    BRANCH_ARRAY(eg_primary_vertex, 3, F)                           \
    BRANCH(eg_ntrial, I)                                            \
    BRANCH(eg_scale_pdf, F)                                         \
    BRANCH(eg_alpha_qcd, F)                                         \
    BRANCH(eg_alpha_qed, F)                                         \
    BRANCH_ARRAY(eg_pdf_id, 2, I)                                   \
    BRANCH_ARRAY(eg_pdf_x, 2, F)                                    \
    BRANCH_ARRAY(eg_pdf_x_pdf, 2, F)                                \
    /* */                                                           \
    BRANCH(ncluster, i)                                             \
    BRANCH_ARRAY(cluster_e, ncluster, F)                            \
    BRANCH_ARRAY(cluster_pt, ncluster, F)                           \
    BRANCH_ARRAY(cluster_eta, ncluster, F)                          \
    BRANCH_ARRAY(cluster_phi, ncluster, F)                          \
    BRANCH_ARRAY2(cluster_lambda_square, ncluster, 2, F)            \
    BRANCH_ARRAY(cluster_tof, ncluster, F)                          \
    BRANCH_ARRAY(cluster_ncell, ncluster, I)                        \
    BRANCH_ARRAY(cluster_cell_id_max, ncluster, s)                  \
    BRANCH_ARRAY(cluster_e_max, ncluster, F)                        \
    BRANCH_ARRAY(cluster_e_cross, ncluster, F)                      \
    BRANCH_ARRAY(cluster_nmc_truth, ncluster, i)                    \
    BRANCH_ARRAY2(cluster_mc_truth_index, ncluster, 32, s)          \
    BRANCH_ARRAY(cluster_iso_tpc_01, ncluster, F)                   \
    BRANCH_ARRAY(cluster_iso_tpc_02, ncluster, F)                   \
    BRANCH_ARRAY(cluster_iso_tpc_03, ncluster, F)                   \
    BRANCH_ARRAY(cluster_iso_tpc_04, ncluster, F)                   \
    BRANCH_ARRAY(cluster_frixione_tpc_04_02, ncluster, F)           \
    BRANCH_ARRAY(cluster_frixione_tpc_04_05, ncluster, F)           \
    BRANCH_ARRAY(cluster_frixione_tpc_04_10, ncluster, F)           \
    BRANCH_ARRAY(cluster_iso_01_truth, ncluster, F)                 \
    BRANCH_ARRAY(cluster_iso_02_truth, ncluster, F)                 \
    BRANCH_ARRAY(cluster_iso_03_truth, ncluster, F)                 \
    BRANCH_ARRAY(cluster_iso_04_truth, ncluster, F)                 \
    BRANCH_ARRAY(cluster_frixione_04_02_truth, ncluster, F)         \
    BRANCH_ARRAY(cluster_frixione_04_05_truth, ncluster, F)         \
    BRANCH_ARRAY(cluster_frixione_04_10_truth, ncluster, F)         \
    BRANCH_ARRAY(cell_e, 17664, F)                                  \
    BRANCH_ARRAY(cell_tof, 17664, F)                                \
    BRANCH_ARRAY(cell_mc_truth_index, 17664, s)                     \
    /* */                                                           \
    BRANCH(ntrack, i)                                               \
    BRANCH_ARRAY(track_e, ntrack, F)                                \
    BRANCH_ARRAY(track_pt, ntrack, F)                               \
    BRANCH_ARRAY(track_eta, ntrack, F)                              \
    BRANCH_ARRAY(track_phi, ntrack, F)                              \
    BRANCH_ARRAY(track_quality, ntrack, b)                          \
    BRANCH_ARRAY(track_tpc_dedx, ntrack, F)                         \
    BRANCH_ARRAY(track_tpc_length_active_zone, ntrack, F)           \
    BRANCH_ARRAY(track_tpc_xrow, ntrack, b)                         \
    BRANCH_ARRAY(track_tpc_ncluster, ntrack, b)                     \
    BRANCH_ARRAY(track_tpc_ncluster_dedx, ntrack, b)                \
    BRANCH_ARRAY(track_tpc_ncluster_findable, ntrack, b)            \
    BRANCH_ARRAY(track_its_ncluster, ntrack, b)                     \
    BRANCH_ARRAY(track_its_chi_square, ntrack, F)                   \
    BRANCH_ARRAY(track_dca_xy, ntrack, F)                           \
    BRANCH_ARRAY(track_dca_z, ntrack, F)                            \
    BRANCH_ARRAY(track_mc_truth_index, ntrack, s)                   \
    BRANCH_ARRAY(track_voronoi_area, ntrack, F)                     \
    /* */                                                           \
    BRANCH(nmuon_track, i)                                          \
    BRANCH_ARRAY(muon_track_e, nmuon_track, F)                      \
    BRANCH_ARRAY(muon_track_pt, nmuon_track, F)                     \
    BRANCH_ARRAY(muon_track_eta, nmuon_track, F)                    \
    BRANCH_ARRAY(muon_track_phi, nmuon_track, F)                    \
    BRANCH_ARRAY(muon_track_r_abs, nmuon_track, F)                  \
    BRANCH_ARRAY(muon_track_p_dca, nmuon_track, F)                  \
    BRANCH_ARRAY(muon_track_sigma_p_dca, nmuon_track, F)            \
    BRANCH_ARRAY(muon_track_delta_sagitta_p, nmuon_track, F)        \
    BRANCH_ARRAY(muon_track_distance_sigma_slope_p, nmuon_track, F) \
    BRANCH_ARRAY(muon_track_mc_truth_index, nmuon_track, s)         \
    /* */                                                           \
    BRANCH(nmc_truth, i)                                            \
    BRANCH_ARRAY(mc_truth_e, nmc_truth, F)                          \
    BRANCH_ARRAY(mc_truth_pt, nmc_truth, F)                         \
    BRANCH_ARRAY(mc_truth_eta, nmc_truth, F)                        \
    BRANCH_ARRAY(mc_truth_phi, nmc_truth, F)                        \
    BRANCH_ARRAY(mc_truth_pdg_code, nmc_truth, S)                   \
    BRANCH_ARRAY(mc_truth_status, nmc_truth, b)                     \
    BRANCH_ARRAY(mc_truth_generator_index, nmc_truth, b)            \
    BRANCH(debug_njet_ue_estimation, i)                             \
    BRANCH_ARRAY(debug_jet_ue_estimation_pt_raw,                    \
                 debug_njet_ue_estimation, F)                       \
    BRANCH_ARRAY(debug_jet_ue_estimation_eta_raw,                   \
                 debug_njet_ue_estimation, F)                       \
    BRANCH_ARRAY(debug_jet_ue_estimation_phi_raw,                   \
                 debug_njet_ue_estimation, F)                       \
    BRANCH_ARRAY(debug_jet_ue_estimation_area_raw,                  \
                 debug_njet_ue_estimation, F)                       \
    /* */                                                           \
    BRANCH(njet, i)                                                 \
    BRANCH_ARRAY(debug_jet_tag_dr_square, njet, F)                  \
    BRANCH_ARRAY(jet_e_raw, njet, F)                                \
    BRANCH_ARRAY(jet_e, njet, F)                                    \
    BRANCH_ARRAY(jet_e_charged, njet, F)                            \
    BRANCH_ARRAY(jet_pt_raw_ue, njet, F)                            \
    BRANCH_ARRAY(jet_pt_raw, njet, F)                               \
    BRANCH_ARRAY(jet_pt, njet, F)                                   \
    BRANCH_ARRAY(jet_pt_charged, njet, F)                           \
    BRANCH_ARRAY(jet_eta_raw, njet, F)                              \
    BRANCH_ARRAY(jet_eta, njet, F)                                  \
    BRANCH_ARRAY(jet_phi, njet, F)                                  \
    BRANCH_ARRAY(jet_area_raw, njet, F)                             \
    BRANCH_ARRAY(jet_area, njet, F)                                 \
    BRANCH_ARRAY(jet_emf_raw, njet, F)                              \
    BRANCH_ARRAY(jet_emf, njet, F)                                  \
    BRANCH_ARRAY(jet_multiplicity_raw, njet, s)                     \
    BRANCH_ARRAY(jet_multiplicity, njet, F)                         \
    BRANCH_ARRAY2(jet_width_sigma_raw, njet, 2, F)                  \
    BRANCH_ARRAY2(jet_width_sigma, njet, 2, F)                      \
    BRANCH_ARRAY(jet_ptd_raw, njet, F)                              \
    BRANCH_ARRAY(jet_ptd, njet, F)                                  \
    BRANCH_ARRAY2(jet_truth_index_z_truth, njet, 2, I)              \
    BRANCH_ARRAY2(jet_truth_z_truth, njet, 2, F)                    \
    BRANCH_ARRAY2(jet_truth_index_z_reco, njet, 2, I)               \
    BRANCH_ARRAY2(jet_truth_z_reco, njet, 2, F)                     \
    BRANCH_ARRAY(jet_e_truth, njet, F)                              \
    BRANCH_ARRAY(jet_pt_truth, njet, F)                             \
    BRANCH_ARRAY(jet_eta_truth, njet, F)                            \
    BRANCH_ARRAY(jet_phi_truth, njet, F)                            \
    BRANCH_ARRAY(jet_area_truth, njet, F)                           \
    BRANCH_ARRAY(jet_emf_truth, njet, F)                            \
    BRANCH_ARRAY(jet_multiplicity_truth, njet, s)                   \
    BRANCH_ARRAY2(jet_width_sigma_truth, njet, 2, F)                \
    BRANCH_ARRAY(jet_ptd_truth, njet, F)                            \
    BRANCH(njet_truth, i)                                           \
    BRANCH_ARRAY(jet_truth_e, njet_truth, F)                        \
    BRANCH_ARRAY(jet_truth_pt, njet_truth, F)                       \
    BRANCH_ARRAY(jet_truth_eta, njet_truth, F)                      \
    BRANCH_ARRAY(jet_truth_phi, njet_truth, F)                      \
    BRANCH_ARRAY(jet_truth_area, njet_truth, F)                     \
    BRANCH_ARRAY(jet_truth_emf, njet_truth, F)                      \
    BRANCH_ARRAY(jet_truth_multiplicity, njet, s)                   \
    BRANCH_ARRAY2(jet_truth_width_sigma, njet, 2, F)                \
    BRANCH_ARRAY(jet_truth_ptd, njet, F)                            \
    BRANCH_ARRAY(met_tpc, 2, D)                                     \
    BRANCH_ARRAY(met_truth, 2, D)                                   \


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

#define ntrigger_class NTRIGGER_CLASS_MAX
#define ncluster NCLUSTER_MAX
#define ntrack NTRACK_MAX
#define nmuon_track NTRACK_MAX
#define nmc_truth NMC_TRUTH_MAX
#define debug_njet_ue_estimation NJET_MAX
#define njet NJET_MAX
#define njet_truth NJET_MAX
#define full_emcal_ncell EMCAL_NCELL

#define BRANCH(b, t)                            \
    t _branch_ ## b;
#define BRANCH_ARRAY(b, d, t)                   \
    t _branch_ ## b [(d)];
#define BRANCH_ARRAY2(b, d, e, t)               \
    t _branch_ ## b [(d)][(e)];
#define BRANCH_STR(b)                           \
    char _branch_ ## b[BUFSIZ];
#define BRANCH_STR_ARRAY(b, d)                  \
    TClonesArray _branch_ ## b;

    MEMBER_BRANCH;

#undef BRANCH
#undef BRANCH_ARRAY
#undef BRANCH_ARRAY2
#undef BRANCH_STR
#undef BRANCH_STR_ARRAY

#undef ntrigger_class
#undef ncluster
#undef ntrack
#undef nmc_truth
#undef debug_njet_ue_estimation
#undef njet
#undef njet_truth
#undef full_emcal_ncell

#undef C
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

    Int_t _run_number_current; //!

    TFormula *_f1_ncluster_tpc_linear_pt_dep; //!
    std::vector<AliESDtrackCuts> _track_cut; //!

    AliEMCALRecoUtils *_reco_util; //!
    AliEMCALGeometry *_emcal_geometry; //!

    AliMuonTrackCuts *_muon_track_cut; //!

    size_t _ncell; //!
    double _skim_cluster_min_e; //!
    double _skim_track_min_pt; //!
    double _skim_muon_track_min_pt; //!
    std::vector<double> _skim_jet_min_pt; //!
    int _skim_multiplicity_tracklet_min_n; //!
    double _stored_track_min_pt; //!
    double _stored_jet_min_pt_raw; //!

    std::vector<bool> _emcal_mask; //!

    void *_emcal_cell_position; //!
    std::vector<double> _emcal_cell_area; //!
    std::vector<std::set<size_t> > _emcal_cell_incident; //!

	void *_keras_model_photon_discrimination;

    AliAnalysisAlien *_alien_plugin; //!
    bool _metadata_filled; //!
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
    void SetAliROOTVersion(const char *version);
    void SetAliPhysicsVersion(const char *version);
    void SetGridDataDir(const char *dir);
    void SetGridDataPattern(const char *pattern);
    //
    void SetSkimClusterMinE(double min_e = -INFINITY);
    void SetSkimTrackMinPt(double min_pt = -INFINITY);
    void SetSkimMuonTrackMinPt(double min_pt = -INFINITY);
    void SetSkimJetMinPt(double min_pt_1 = -INFINITY,
                         double min_pt_2 = -INFINITY,
                         double min_pt_3 = -INFINITY);
    void SetSkimMultiplicityTrackletMinN(int min_n = INT_MIN);
    //
    void SetStoredTrackMinPt(double min_pt = -INFINITY);
    void SetStoredJetMinPtRaw(double min_pt_raw = -INFINITY);
    ClassDef(AliAnalysisTaskNTGJ, 1);
};

#endif // ALIANALYSISTASKPHOTONDISC_H_
