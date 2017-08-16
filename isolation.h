// -*- mode: c++; -*-

#ifndef ISOLATION_H_
#define ISOLATION_H_

#include <algorithm>

namespace {

	// ROOT ACLiC weirdness is preventing delta_vs_iso from being
	// declared const

	double frixione_iso_max_x_e_eps(std::vector<std::pair<
									double, double> >
									delta_vs_iso,
									double delta_0, double n,
									std::vector<double> ue_estimate =
									std::vector<double>(),
									double phi = NAN)
	{
		std::sort(delta_vs_iso.begin(), delta_vs_iso.end());

		double sum_iso = 0;
		double max_x_e_eps = -INFINITY;
		const double one_cos_delta_0 = 1 - cos(delta_0);

		for (std::vector<std::pair<double, double> >::
				 const_iterator iterator = delta_vs_iso.begin();
			 iterator != delta_vs_iso.end(); iterator++) {
			const double delta = iterator->first;
			const double iso_ue =
				evaluate_ue(ue_estimate, phi, delta) *
				M_PI * std::pow(delta, 2);

			sum_iso += iterator->second;

			const double x_e_eps = (sum_iso - iso_ue) *
				pow(one_cos_delta_0 / (1 - cos(delta)), n);

			max_x_e_eps = std::max(max_x_e_eps, x_e_eps);
		}

		return max_x_e_eps;
	}

}

#endif // ISOLATION_H_
