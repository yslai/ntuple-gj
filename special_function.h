// -*- mode: c++; -*-

#ifndef SPECIAL_FUNCTION_H_
#define SPECIAL_FUNCTION_H_

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

}

#endif // SPECIAL_FUNCTION_H_
