// -*- mode: c++; -*-

#ifndef EMCAL_H_
#define EMCAL_H_

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

}

#endif // EMCAL_H_
