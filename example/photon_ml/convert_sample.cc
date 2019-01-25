#include <numeric>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TProfile.h>

#include <H5Cpp.h>

// This is chosen to be the CPU L2 cache size, which should exceed 512
// kB for many years now
#ifndef HDF5_DEFAULT_CACHE
#define HDF5_DEFAULT_CACHE (512 * 1024)
#endif // HDF5_DEFAULT_CACHE

#ifndef HDF5_USE_DEFLATE
#define HDF5_USE_DEFLATE
#endif // HDF5_USE_DEFLATE

#define RANK 2

#include <special_function.h>

namespace {

    void to_sm_nphi(unsigned int &sm, unsigned int &nphi,
                    unsigned int n)
    {
        sm = n < 11520 ? n / 1152 :
            n < 12288 ? 10 + (n - 11520) / 384 :
            n < 16896 ? 12 + (n - 12288) / 768 :
            18 + (n - 16896) / 384;
        nphi = sm < 10 ? 24 : sm < 12 ? 8 : sm < 18 ? 24 : 8;
    }

    void to_sm_ieta_iphi(unsigned int &sm, unsigned int &ieta,
                         unsigned int &iphi, unsigned int n)
    {
        unsigned int nphi;

        to_sm_nphi(sm, nphi, n);

        const unsigned int n0 =
            sm < 10 ? sm * 1152 :
            sm < 12 ? 11520 + (sm - 10) * 384 :
            sm < 18 ? 12288 + (sm - 12) * 768 :
            16896 + (sm - 18) * 384;
        const unsigned int n1 = n - n0;

        ieta = 2 * (n1 / (2 * nphi)) + 1 - (n1 % 2);
        iphi = (n1 / 2) % nphi;
    }

    void neta_nphi(unsigned int &neta, unsigned int &nphi,
                    const unsigned int sm)
    {
        neta = sm < 12 ? 48 : sm < 18 ? 32 : 48;
        nphi = sm < 10 ? 24 : sm < 12 ? 8 : sm < 18 ? 24 : 8;
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

    void cell_5_5(unsigned int n_5_5[], const unsigned int n,
                  const unsigned int ld = 5)
    {
        const unsigned int sm = n < 11520 ? n / 1152 :
            n < 12288 ? 10 + (n - 11520) / 384 :
            n < 16896 ? 12 + (n - 12288) / 768 :
            18 + (n - 16896) / 384;
        const unsigned int nphi =
            sm < 10 ? 24 : sm < 12 ? 8 : sm < 18 ? 24 : 8;

        n_5_5[0 * ld + 0] = n - 2 * nphi - 4;
        n_5_5[0 * ld + 1] = n - 2 * nphi - 2;
        n_5_5[0 * ld + 2] = n - 2 * nphi;
        n_5_5[0 * ld + 3] = n - 2 * nphi + 2;
        n_5_5[0 * ld + 4] = n - 2 * nphi + 4;
        if (n % 2 == 0) {
            n_5_5[1 * ld + 0] = n - 3;
            n_5_5[1 * ld + 1] = n - 1;
            n_5_5[1 * ld + 2] = n + 1;
            n_5_5[1 * ld + 3] = n + 3;
            n_5_5[1 * ld + 4] = n + 5;
        }
        else {
            n_5_5[1 * ld + 0] = n - 2 * nphi - 5;
            n_5_5[1 * ld + 1] = n - 2 * nphi - 3;
            n_5_5[1 * ld + 2] = n - 2 * nphi - 1;
            n_5_5[1 * ld + 3] = n - 2 * nphi + 1;
            n_5_5[1 * ld + 4] = n - 2 * nphi + 3;
        }
        n_5_5[2 * ld + 0] = n - 4;
        n_5_5[2 * ld + 1] = n - 2;
        n_5_5[2 * ld + 2] = n;
        n_5_5[2 * ld + 3] = n + 2;
        n_5_5[2 * ld + 4] = n + 4;
        if (n % 2 == 0) {
            n_5_5[3 * ld + 0] = n + 2 * nphi - 3;
            n_5_5[3 * ld + 1] = n + 2 * nphi - 1;
            n_5_5[3 * ld + 2] = n + 2 * nphi + 1;
            n_5_5[3 * ld + 3] = n + 2 * nphi + 3;
            n_5_5[3 * ld + 4] = n + 2 * nphi + 5;
        }
        else {
            n_5_5[3 * ld + 0] = n - 5;
            n_5_5[3 * ld + 1] = n - 3;
            n_5_5[3 * ld + 2] = n - 1;
            n_5_5[3 * ld + 3] = n + 1;
            n_5_5[3 * ld + 4] = n + 3;
        }
        n_5_5[4 * ld + 0] = n + 2 * nphi - 4;
        n_5_5[4 * ld + 1] = n + 2 * nphi - 2;
        n_5_5[4 * ld + 2] = n + 2 * nphi;
        n_5_5[4 * ld + 3] = n + 2 * nphi + 2;
        n_5_5[4 * ld + 4] = n + 2 * nphi + 4;
    }

    void cell_neighbor(unsigned int n_m_m[], const unsigned int n,
                       const unsigned int m = 5, unsigned int ld = 0)
    {
        if (ld == 0) {
            ld = m;
        }

        const unsigned int sm = n < 11520 ? n / 1152 :
            n < 12288 ? 10 + (n - 11520) / 384 :
            n < 16896 ? 12 + (n - 12288) / 768 :
            18 + (n - 16896) / 384;
        const unsigned int nphi =
            sm < 10 ? 24 : sm < 12 ? 8 : sm < 18 ? 24 : 8;

        for (unsigned int i = 0; i < m; i++) {
            const unsigned int i_centered = i - (m - 1) / 2;
            const unsigned int offset_i = (i_centered & 1) == 0 ?
                i_centered * nphi :
                (n & 1) == 0 ?
                ((i_centered + 2) & ~1) * nphi + 1 :
                (i_centered & ~1) * nphi - 1;
            for (unsigned int j = 0; j < m; j++) {
                const unsigned int j_times_2_centered =
                    2 * j - (m - 1);

                n_m_m[i * ld + j] =
                    n + offset_i + j_times_2_centered;
            }
        }
    }

    float unit_step(float x)
    {
        return std::isfinite(x) ? x >= 0.0F ? 1.0F : 0.0F : x;
    }

    void cluster_center_recal(float &pseudorapidity, float &azimuth)
    {
        static const float alice_emcal_granularity = 0.0143F;
        const float p = pseudorapidity - alice_emcal_granularity * (
	        6.26120e-03F -
            7.93806e-02F * unit_step(azimuth - 1.745F) +
            7.01157e-02F * unit_step(azimuth - 2.094F) -
            5.67665e-02F * unit_step(azimuth - 2.443F) -
            3.12072e-02F * unit_step(azimuth - 2.792F) -
            0.0681282 + 0.30579 * pseudorapidity +
            0.174691 * std::pow(pseudorapidity, 2));
        const float a = azimuth - alice_emcal_granularity * (
	        (1.80892e-01F + 3.88693e-02F + 2.17416e-02F + 9.48413e-03F) -
            (2.88726e-01F + 7.01681e-02F + 1.31708e-02F - 7.20175e-04F) *
            unit_step(azimuth - 1.745F) -
            (1.77386e-01F + 8.19726e-02F + 4.11356e-02F + 1.04680e-03F) *
            unit_step(azimuth - 2.094F) -
            (1.38388e-01F + 3.72398e-02F + 2.81211e-03F + 7.69352e-03F) *
            unit_step(azimuth - 2.443F) -
            (1.02816e-01F + 3.30964e-02F - 8.73742e-03F + 6.21222e-03F) *
            unit_step(azimuth - 2.792F) -
            (0.00410491F + 0.000346243F) -
            (0.0682701F + 0.0270903F) * pseudorapidity -
            (0.060021F - 0.0362166F) * std::pow(pseudorapidity, 2));

        pseudorapidity = p;
        azimuth = a;
    }

    float min_sm_dazimuth(float pseudorapidity, float azimuth)
    {
        static const float sm_center[20] = {
            1.57079583F, 1.57079670F, 1.91986285F, 1.91986059F,
            2.26892674F, 2.26892899F, 2.61799219F, 2.61799514F,
            2.96706059F, 2.96705799F, 3.20794364F, 3.20839677F,
            -1.57079635F, -1.57079688F, -1.22173021F, -1.22173229F,
            -0.87266589F, -0.87266172F, -0.63178698F, -0.63132969F
        };
        float min_dazimuth = M_PI;

        for (size_t i = 0; i < 20; i++) {
            if ((i % 2 == 0 ? 1.0F : -1.0F) * pseudorapidity >= 0) {
                const float dazimuth =
                    angular_range_reduce(azimuth - sm_center[i]);

                if (fabsf(dazimuth) < min_dazimuth) {
                    min_dazimuth = dazimuth;
                }
            }
        }

        return min_dazimuth;
    }

}

#define NTRACK_MAX (1U << 15)
#define NCLUSTER_MAX (1U << 15)
#define NMC_TRUTH_MAX (1U << 15)
#define CLUSTER_NMC_TRUTH_MAX 32

int main(int argc, char *argv[])
{
    unsigned int id[48][24];

#if 1
    for (unsigned int n = 0; ; n++) {
        unsigned int sm;
        unsigned int ieta;
        unsigned int iphi;

        to_sm_ieta_iphi(sm, ieta, iphi, n);

        if (sm >= 1) {
            break;
        }
        id[ieta][iphi] = n;
    }
    for (size_t i = 0; i < 10; i++) {
        for (size_t j = 6; j < 9 + 6; j++) {
            fprintf(stderr, "%3u ", id[i][j]);
        }
        fprintf(stderr, "\n");
    }

    for (unsigned int n = 117; n >= 116; n--) {
        unsigned int n_9_9[7 * 7];

        cell_neighbor(n_9_9, n, 7);
        for (size_t i = 0; i < 7; i++) {
            for (size_t j = 0; j < 7; j++) {
                fprintf(stderr, "%3u ", n_9_9[i * 7 + j]);
            }
            fprintf(stderr, "< %3u\n", n);
        }
    }
#endif

    if (argc < 3) {
        exit(EXIT_FAILURE);
    }

    // Access mode H5F_ACC_TRUNC truncates any existing file, while
    // not throwing any exception (unlike H5F_ACC_RDWR)
    fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, argv[argc - 1]);
    H5::H5File hdf5_file(argv[argc - 1], H5F_ACC_TRUNC);
    // How many properties per photon is written
#define NCELL 5
    static const size_t row_size_X = NCELL * NCELL + 4;
    static const size_t row_size_y = 1;
    static const size_t row_size_lam = 1;
    // The tensor dimension increment for each new event
    hsize_t dim_extend_X[RANK] = { 1, row_size_X };
    hsize_t dim_extend_y[RANK] = { 1, row_size_y };
    hsize_t dim_extend_lam[RANK] = { 1, row_size_lam };
    // The maximum tensor dimension, for unlimited number of events
    hsize_t dim_max_X[RANK] = { H5S_UNLIMITED, row_size_X };
    hsize_t dim_max_y[RANK] = { H5S_UNLIMITED, row_size_y };
    hsize_t dim_max_lam[RANK] = { H5S_UNLIMITED, row_size_lam };
    // The extensible HDF5 data space
	H5::DataSpace data_space_X(RANK, dim_extend_X, dim_max_X);
	H5::DataSpace data_space_y(RANK, dim_extend_y, dim_max_y);
	H5::DataSpace data_space_lam(RANK, dim_extend_lam, dim_max_lam);

    // To enable zlib compression (there will be many NANs) and
    // efficient chunking (splitting of the tensor into contingous
    // hyperslabs), a HDF5 property list is needed
    H5::DSetCreatPropList property_X = H5::DSetCreatPropList();
    H5::DSetCreatPropList property_y = H5::DSetCreatPropList();
    H5::DSetCreatPropList property_lam = H5::DSetCreatPropList();

#ifdef HDF5_USE_DEFLATE
    // Check for zlib (deflate) availability and enable only if
    // present
    if (!H5Zfilter_avail(H5Z_FILTER_DEFLATE)) {
        fprintf(stderr, "%s:%d: warning: deflate filter not "
                "available\n", __FILE__, __LINE__);
    }
    else {
        unsigned int filter_info;

        H5Zget_filter_info(H5Z_FILTER_DEFLATE, &filter_info);
        if (!(filter_info & H5Z_FILTER_CONFIG_ENCODE_ENABLED)) {
            fprintf(stderr, "%s:%d: warning: deflate filter not "
                    "available for encoding\n", __FILE__, __LINE__);
        }
        else {
            property_X.setDeflate(1);
            property_y.setDeflate(1);
            property_lam.setDeflate(1);
        }
    }
#endif // HDF5_USE_DEFLATE

    // Activate chunking, while observing the HDF5_DEFAULT_CACHE being
    // the CPU L2 cache size
    hsize_t dim_chunk_X[RANK] = {
        std::max(static_cast<unsigned long long>(1),
                 HDF5_DEFAULT_CACHE /
                 std::max(static_cast<unsigned long long>(1),
                          dim_extend_X[1] * sizeof(float))),
        dim_extend_X[1]
    };
    hsize_t dim_chunk_y[RANK] = {
        std::max(static_cast<unsigned long long>(1),
                 HDF5_DEFAULT_CACHE /
                 std::max(static_cast<unsigned long long>(1),
                          dim_extend_y[1] * sizeof(float))),
        dim_extend_y[1]
    };
    hsize_t dim_chunk_lam[RANK] = {
        std::max(static_cast<unsigned long long>(1),
                 HDF5_DEFAULT_CACHE /
                 std::max(static_cast<unsigned long long>(1),
                          dim_extend_lam[1] * sizeof(float))),
        dim_extend_lam[1]
    };

    property_X.setChunk(RANK, dim_chunk_X);
    property_y.setChunk(RANK, dim_chunk_y);
    property_lam.setChunk(RANK, dim_chunk_lam);

    // Create the data set, which will have space for the first event
    H5::DataSet data_set_X =
        hdf5_file.createDataSet("X", H5::PredType::NATIVE_FLOAT,
                                data_space_X, property_X);
    H5::DataSet data_set_y =
        hdf5_file.createDataSet("y", H5::PredType::NATIVE_FLOAT,
                                data_space_y, property_y);
    H5::DataSet data_set_lam =
        hdf5_file.createDataSet("lam", H5::PredType::NATIVE_FLOAT,
                                data_space_lam, property_lam);
    hsize_t offset[RANK] = {0, 0};

    size_t count_prompt = 0;
    size_t count_nonprompt = 0;

    for (int iarg = 1; iarg < argc - 1; iarg++) {
        TFile *root_file = TFile::Open(argv[iarg]);

        if (root_file == NULL) {
            continue;
        }

        TDirectoryFile *df = dynamic_cast<TDirectoryFile *>
            (root_file->Get("AliAnalysisTaskNTGJ"));

        if (df == NULL) {
            continue;
        }

        TTree *_tree_event = dynamic_cast<TTree *>
            (df->Get("_tree_event"));

        if (_tree_event == NULL) {
            continue;
        }

        fprintf(stderr, "%s:%d: %s (%d / %d)\n", __FILE__, __LINE__,
                argv[iarg], iarg - 1, argc - 2);

        Float_t multiplicity_v0[64];

        _tree_event->SetBranchAddress("multiplicity_v0", multiplicity_v0);

        UInt_t ncluster;
        Float_t cluster_e[NCLUSTER_MAX];
        Float_t cluster_pt[NCLUSTER_MAX];
        Float_t cluster_eta[NCLUSTER_MAX];
        Float_t cluster_phi[NCLUSTER_MAX];
        Float_t cluster_lambda_square[NCLUSTER_MAX][2];
        Float_t cluster_tof[NCLUSTER_MAX];
        Int_t cluster_ncell[NCLUSTER_MAX];
        UShort_t cluster_cell_id_max[NCLUSTER_MAX];
        Float_t cluster_e_max[NCLUSTER_MAX];
        Float_t cluster_e_cross[NCLUSTER_MAX];
        Float_t cluster_s_nphoton[NCLUSTER_MAX][4];
        UInt_t cluster_nmc_truth[NCLUSTER_MAX];
        UShort_t cluster_mc_truth_index[NCLUSTER_MAX][CLUSTER_NMC_TRUTH_MAX];

        _tree_event->SetBranchAddress("ncluster", &ncluster);
        _tree_event->SetBranchAddress("cluster_e", cluster_e);
        _tree_event->SetBranchAddress("cluster_pt", cluster_pt);
        _tree_event->SetBranchAddress("cluster_eta", cluster_eta);
        _tree_event->SetBranchAddress("cluster_phi", cluster_phi);
        _tree_event->SetBranchAddress("cluster_lambda_square",
                                      cluster_lambda_square);
        _tree_event->SetBranchAddress("cluster_tof", cluster_tof);
        _tree_event->SetBranchAddress("cluster_ncell",
                                      cluster_ncell);
        _tree_event->SetBranchAddress("cluster_cell_id_max",
                                      cluster_cell_id_max);
        _tree_event->SetBranchAddress("cluster_e_max",
                                      cluster_e_max);
        _tree_event->SetBranchAddress("cluster_e_cross",
                                      cluster_e_cross);
        _tree_event->SetBranchAddress("cluster_nmc_truth",
                                      cluster_nmc_truth);
        _tree_event->SetBranchAddress("cluster_mc_truth_index",
                                      cluster_mc_truth_index);

        UInt_t ntrack;
        Float_t track_e[NTRACK_MAX];
        Float_t track_pt[NTRACK_MAX];
        Float_t track_eta[NTRACK_MAX];
        Float_t track_phi[NTRACK_MAX];
        UChar_t track_quality[NTRACK_MAX];
        Float_t track_eta_emcal[NTRACK_MAX];
        Float_t track_phi_emcal[NTRACK_MAX];

        _tree_event->SetBranchAddress("ntrack", &ntrack);
        _tree_event->SetBranchAddress("track_e", track_e);
        _tree_event->SetBranchAddress("track_pt", track_pt);
        _tree_event->SetBranchAddress("track_eta", track_eta);
        _tree_event->SetBranchAddress("track_phi", track_phi);
        _tree_event->SetBranchAddress("track_quality",
                                      track_quality);
        _tree_event->SetBranchAddress("track_eta_emcal",
                                      track_eta_emcal);
        _tree_event->SetBranchAddress("track_phi_emcal",
                                      track_phi_emcal);

        Float_t cell_e[17664];
        Float_t cell_tof[17664];

        _tree_event->SetBranchAddress("cell_e", cell_e);
        _tree_event->SetBranchAddress("cell_tof", cell_tof);

        UInt_t nmc_truth;
        Float_t mc_truth_e[NMC_TRUTH_MAX];
        Short_t mc_truth_first_parent_pdg_code[NMC_TRUTH_MAX];

        _tree_event->SetBranchAddress("nmc_truth", &nmc_truth);
        _tree_event->SetBranchAddress("mc_truth_e", mc_truth_e);
        _tree_event->SetBranchAddress("mc_truth_first_parent_pdg_code",
                                      mc_truth_first_parent_pdg_code);

        size_t cluster_count = 0;

        for (Long64_t i = 0; i < _tree_event->GetEntries(); i++) {
            _tree_event->GetEntry(i);
            for (UInt_t j = 0; j < ncluster; j++) {
                if (cluster_e[j] >= 8.0F &&
                    inside_edge(cluster_cell_id_max[j],
                                (NCELL - 1) / 2) &&
                    // EMCAL noise suppression cut (with no particular
                    // name by the ALICE experiment)
                    cluster_ncell[j] >= 2 &&
                    // EMCAL noise suppression cut, which the ALICE
                    // collaboration calls the "exoticity cut"
                    cluster_e_cross[j] > 0.05 * cluster_e_max[j]) {
                    unsigned int neighbor[NCELL * NCELL];

                    cell_neighbor(neighbor, cluster_cell_id_max[j],
                                  NCELL);

                    // Loop over all MC (ground) truth particles and
                    // decide if the particle is prompt by the highest
                    // energy one hitting the cluster
                    float mc_truth_e_max = -INFINITY;
                    bool prompt = false;

                    for (size_t k = 0; k < cluster_nmc_truth[j]; k++) {
                        if (mc_truth_e[cluster_mc_truth_index[j][k]] >
                            mc_truth_e_max) {
#if 0
                            fprintf(stderr, "%s:%d: %d\n",
                                    __FILE__, __LINE__,
                                    mc_truth_first_parent_pdg_code
                                    [cluster_mc_truth_index[j][k]]);
#endif
                            // Prompt clusters are from partons, gauge
                            // bosons and leptons (electrons). The PDG
                            // numbering scheme of Monte Carlo
                            // generator particles
                            // (http://pdg.lbl.gov/mc_particle_id_contents.html)
                            // designates partons, gauge bosons, and
                            // leptons to have absolute values < 100
                            // and mesons to have absolute values >=
                            // 100
                            prompt = std::abs(
                                mc_truth_first_parent_pdg_code
                                [cluster_mc_truth_index[j][k]])
                                < 100;
                            mc_truth_e_max = mc_truth_e
                                  [cluster_mc_truth_index[j][k]];
                        }
                    }

                    // Guard against clusters with empty MC (ground)
                    // truth particle association
                    if (mc_truth_e_max >= 0) {
                        std::vector<float> X;
                        std::vector<unsigned int> y;

                        for (size_t k = 0; k < NCELL * NCELL; k++) {
                            X.push_back(
                                std::isfinite(cell_e[neighbor[k]]) ?
                                cell_e[neighbor[k]] / cluster_e[j] : 0);
                        }
                        X.push_back(1.0F / sqrtf(cluster_e[j]));
                        X.push_back(cluster_eta[j]);
                        X.push_back(min_sm_dazimuth(
                             cluster_eta[j], cluster_phi[j]));
                        X.push_back(std::accumulate(
                            multiplicity_v0,
                            multiplicity_v0 + 64, 0));
                        y.push_back(prompt ? 1 : 0);

                        std::vector<float> lam(1, cluster_lambda_square[j][0]);

                        if (offset[0] == 0) {
                            data_set_X.write(&X[0], H5::PredType::
                                             NATIVE_FLOAT);
                            data_set_y.write(&y[0], H5::PredType::
                                             NATIVE_UINT);
                            data_set_lam.write(&lam[0], H5::PredType::
                                               NATIVE_FLOAT);
                        }
                        else {
                            const hsize_t dim_extended_X[RANK] = {
                                offset[0] + dim_extend_X[0],
                                dim_extend_X[1]
                            };
                            const hsize_t dim_extended_y[RANK] = {
                                offset[0] + dim_extend_y[0],
                                dim_extend_y[1]
                            };
                            const hsize_t dim_extended_lam[RANK] = {
                                offset[0] + dim_extend_lam[0],
                                dim_extend_lam[1]
                            };

                            data_set_X.extend(dim_extended_X);
                            data_set_y.extend(dim_extended_y);
                            data_set_lam.extend(dim_extended_lam);

                            H5::DataSpace file_space;
                            H5::DataSpace memory_space;

                            file_space = data_set_X.getSpace();
                            file_space.selectHyperslab(
                                H5S_SELECT_SET, dim_extend_X,
                                offset);
                            memory_space =
                                H5::DataSpace(RANK, dim_extend_X,
                                              NULL);
                            data_set_X.write(&X[0], H5::PredType::
                                             NATIVE_FLOAT,
                                             memory_space,
                                             file_space);

                            file_space = data_set_y.getSpace();
                            file_space.selectHyperslab(
                                H5S_SELECT_SET, dim_extend_y,
                                offset);
                            memory_space =
                                H5::DataSpace(RANK, dim_extend_y,
                                              NULL);
                            data_set_y.write(&y[0], H5::PredType::
                                             NATIVE_UINT,
                                             memory_space,
                                             file_space);

                            file_space = data_set_lam.getSpace();
                            file_space.selectHyperslab(
                                H5S_SELECT_SET, dim_extend_lam,
                                offset);
                            memory_space =
                                H5::DataSpace(RANK, dim_extend_lam,
                                              NULL);
                            data_set_lam.write(&lam[0], H5::PredType::
                                               NATIVE_UINT,
                                               memory_space,
                                               file_space);
                        }
                        offset[0]++;
                        if (prompt) {
                            count_prompt++;
                        }
                        else {
                            count_nonprompt++;
                        }
                    }
                }
            }
            if (i % 1000 == 0) {
                fprintf(stderr, "%s:%d: %lld / %lld "
                        "(total cluster count = %lld, "
                        "prompt = %lu, nonprompt = %lu)\n",
                        __FILE__, __LINE__, i,
                        _tree_event->GetEntries(), offset[0],
                        count_prompt, count_nonprompt);
            }
#if 0
            // Abort after certain number of events
            if (i == 10000) {
                break;
            }
#endif
        }
        _tree_event->Delete();

        delete df;

        root_file->Close();

        delete root_file;
    }

    hdf5_file.close();

    return EXIT_SUCCESS;
}
