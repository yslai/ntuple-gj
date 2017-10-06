// -*- mode: c++; -*-

#ifndef EMCAL_H_
#define EMCAL_H_

#include <vector>
#include <AliVCluster.h>
#include "keras_model.h"

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

    void cell_3_3(unsigned int n_3_3[], const unsigned int n,
                  const unsigned int ld = 3)
    {
        unsigned int sm;
        unsigned int nphi;

        to_sm_nphi(sm, nphi, n);

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

    void cell_cross(unsigned int n_cross[], const unsigned int n)
    {
        // Note that the "cross" happens to correspond to the odd
        // entries in the 3x3 array

        unsigned int sm;
        unsigned int nphi;

        to_sm_nphi(sm, nphi, n);

        n_cross[0] = n % 2 == 0 ? n + 1 : n - 2 * nphi - 1;
        n_cross[1] = n - 2;
        n_cross[2] = n + 2;
        n_cross[3] = n % 2 == 0 ? n + 2 * nphi + 1 : n - 1;
    }

    void cell_max_cross(Int_t &cell_id_max,
                        Double_t &cell_energy_max,
                        Double_t &energy_cross,
                        AliVCluster *c, AliVCaloCells *emcal_cell)
    {
        std::set<Int_t> cluster_cell_id;

        cell_id_max = -1;
        cell_energy_max = -INFINITY;

        for (Int_t j = 0; j < c->GetNCells(); j++) {
            const Int_t cell_id = c->GetCellsAbsId()[j];
            const Double_t cell_energy =
                emcal_cell->GetCellAmplitude(cell_id);

            cluster_cell_id.insert(cell_id);
            if (cell_energy > cell_energy_max) {
                cell_energy_max = cell_energy;
                cell_id_max = cell_id;
            }
        }

        energy_cross = NAN;

        if (cell_id_max != -1) {
            unsigned int cell_id_max_cross[4];

            cell_cross(cell_id_max_cross, cell_id_max);
            energy_cross = 0;
            for (size_t j = 0; j < 4; j++) {
                if (cluster_cell_id.find(cell_id_max_cross[j]) !=
                    cluster_cell_id.end()) {
                    energy_cross += emcal_cell->
                        GetCellAmplitude(cell_id_max_cross[j]);
                }
            }
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

    std::vector<float>
    cluster_cell_keras_inference(AliVCluster *cluster,
                                 AliVCaloCells *cell,
                                 const double *vertex,
                                 const AliVVZERO *v0,
                                 KerasModel &model)
    {
        Int_t cell_id_max = -1;
        Double_t cell_energy_max = -INFINITY;
        Double_t cell_cross = NAN;

        cell_max_cross(cell_id_max, cell_energy_max, cell_cross,
                       cluster, cell);

        unsigned int sm_max;
        unsigned int nphi_max;

        to_sm_nphi(sm_max, nphi_max, cell_id_max);

        if (!inside_edge(cell_id_max, 2)) {
            return std::vector<float>();
        }

        TLorentzVector p;

        cluster->GetMomentum(p, vertex);

        unsigned int cell_id_5_5[25];

        cell_5_5(cell_id_5_5, cell_id_max);

        Tensor feature(29);

        for (size_t i = 0; i < 25; i++) {
            unsigned int sm;
            unsigned int nphi;

            to_sm_nphi(sm, nphi, cell_id_5_5[i]);
            feature(i) = sm == sm_max ?
                cell->GetCellAmplitude(cell_id_5_5[i]) / p.E() :
                0;
        }
        feature(25) = 1 / sqrt(p.E());
        feature(26) = p.Eta();

        static const double sm_center[20][2] = {
            {  1,  1.5707958333333356 },
            { -1,  1.570796701388894  },
            {  1,  1.9198628472222254 },
            { -1,  1.9198605902777761 },
            {  1,  2.2689267361111165 },
            { -1,  2.2689289930555647 },
            {  1,  2.6179921875000085 },
            { -1,  2.617995138888896  },
            {  1,  2.9670605902777836 },
            { -1,  2.9670579861111235 },
            {  1,  3.2079436405129194 },
            { -1,  3.2083967655129197 },
            {  1, -1.5707963541666699 },
            { -1, -1.5707968750000028 },
            {  1, -1.2217302083333352 },
            { -1, -1.2217322916666704 },
            {  1, -0.8726658854166689 },
            { -1, -0.8726617187500013 },
            {  1, -0.6317869791666666 },
            { -1, -0.6313296875       }
        };
        double min_dazimuth = M_PI;

        for (size_t i = 0; i < 20; i++) {
            if (sm_center[i][0] * p.Eta() >= 0) {
                const double dazimuth = angular_range_reduce(
                    p.Phi() - sm_center[i][1]);
                if (fabs(dazimuth) < min_dazimuth) {
                    min_dazimuth = dazimuth;
                }
            }
        }
        feature(27) = min_dazimuth;

        double sum_v0 = 0;

        for (size_t i = 0; i < 64; i++) {
            sum_v0 += v0->GetMultiplicity(i);
        }
        feature(28) = sum_v0;

        static const float mean[29] = {
0.0069592884397400656, 0.016739319288611606, 0.017413634470062438, 0.016836823271191431, 0.0062490905766574402, 0.011579136673917802, 0.02108443031761912, 0.036355012949252585, 0.019960526175031221, 0.0096304383968745313, 0.012721679201750422, 0.048735977286154052, 0.55877568122273502, 0.049128016675045558, 0.011427085419368051, 0.010574512708817985, 0.019122197882992662, 0.033512532793995589, 0.018437010084581661, 0.0091730716482323824, 0.0059668110213861852, 0.014275452826075318, 0.014208260230330783, 0.014007044965148305, 0.0057210942874272057, 0.31024828824296874, -0.0047336070862199816, 0.00037758976077924402, 5359.7114849607315
        };
        static const float std[29] = {
0.03220330143028375, 0.058527186459927376, 0.058980584281578763, 0.058814659977434826, 0.030467882695618582, 0.045897783119451321, 0.057881013595359083, 0.060718691525402518, 0.056007461153791016, 0.040591070168018495, 0.046286403913191096, 0.075668568566528357, 0.17227619105109979, 0.075823573806928821, 0.043645304093659625, 0.043893754912312719, 0.054684955826642324, 0.057795085444955568, 0.053928497822212026, 0.039200660225381219, 0.02949374176044248, 0.054212090352210654, 0.053193423162352239, 0.053465619132931054, 0.02941922773989222, 0.028035783584958934, 0.39410486198367362, 0.077305135107501324, 6986.2733548371962
        };

        for (size_t i = 0; i < 29; i++) {
            feature(i) = (feature(i) - mean[i]) / std[i];
        }

        Tensor output_tensor;

        model.Apply(&feature, &output_tensor);

        return output_tensor.data_;
    }

}

#endif // EMCAL_H_
