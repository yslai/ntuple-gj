// -*- mode: c++; -*-

#ifndef EMCAL_CELL_H_
#define EMCAL_CELL_H_

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

}

#endif // EMCAL_CELL_H_
