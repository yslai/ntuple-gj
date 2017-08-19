// -*- mode: c++; -*-

#ifndef JET_H_
#define JET_H_

#include <vector>
#include <set>
#include <map>

#include <TDecompSVD.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>
#include <CGAL/Polygon_2.h>

#include <fastjet/PseudoJet.hh>
#pragma GCC diagnostic pop

namespace {

    typedef CGAL::Delaunay_triangulation_2<
        CGAL::Exact_predicates_inexact_constructions_kernel>
    delaunay_triangulation_t;
    typedef delaunay_triangulation_t::Point point_2d_t;
    typedef CGAL::Voronoi_diagram_2<
        delaunay_triangulation_t,
        CGAL::Delaunay_triangulation_adaptation_traits_2<
            delaunay_triangulation_t>,
        CGAL::
Delaunay_triangulation_caching_degeneracy_removal_policy_2<
            delaunay_triangulation_t> > voronoi_diagram_t;
    typedef CGAL::Polygon_2<
        CGAL::Exact_predicates_inexact_constructions_kernel>
        polygon_t;

    void voronoi_area_incident(
        std::vector<double> &particle_area,
        std::vector<std::set<size_t> > &particle_incident,
        const std::vector<point_2d_t> &
        particle_pseudorapidity_azimuth)
    {
        // Make the Voronoi diagram

        voronoi_diagram_t diagram;

        // Reverse Voronoi face lookup
        std::map<voronoi_diagram_t::Face_handle, size_t> face_index;

        for (std::vector<point_2d_t>::const_iterator iterator =
                 particle_pseudorapidity_azimuth.begin();
             iterator != particle_pseudorapidity_azimuth.end();
             iterator++) {
            // Reflect at ALICE TPC boundary of |eta| = 0.9, to
            // cut-off the tesselation at the boundary condition via
            // "mirror tracks"
            for (int j = -1; j <= 1; j++) {
                // Make two additional replicas with azimuth +/- 2 pi
                // (and use only the middle) to mimick the cyclical
                // boundary condition
                for (int k = -1; k <= 1; k++) {
                    const point_2d_t
                        p((2 * (j & 1) - 1) * iterator->x() +
                          j * (2 * 0.9),
                          iterator->y() + k * (2 * M_PI));
                    const voronoi_diagram_t::Face_handle
                        handle = diagram.insert(p);

                    face_index[handle] = iterator -
                        particle_pseudorapidity_azimuth.begin();
                }
            }
        }

        particle_area.clear();
        particle_incident = std::vector<std::set<size_t> >(
            particle_pseudorapidity_azimuth.size(),
            std::set<size_t>());

        // Extract the Voronoi cells as polygon and calculate the
        // area associated with individual particles

        for (std::vector<point_2d_t>::const_iterator iterator =
                 particle_pseudorapidity_azimuth.begin();
             iterator != particle_pseudorapidity_azimuth.end();
             iterator++) {
            const voronoi_diagram_t::Locate_result result =
                diagram.locate(*iterator);
            const voronoi_diagram_t::Face_handle *face =
                boost::get<voronoi_diagram_t::Face_handle>(&result);
            double polygon_area;

            if (face != NULL) {
                voronoi_diagram_t::Ccb_halfedge_circulator
                    circulator_start = (*face)->outer_ccb();
                bool unbounded = false;
                polygon_t polygon;

                voronoi_diagram_t::Ccb_halfedge_circulator
                    circulator = circulator_start;

                // Circle around the edges and extract the polygon
                // vertices
                do {
                    if (circulator->has_target()) {
                        polygon.push_back(
                            circulator->target()->point());
                        particle_incident[face_index[*face]].insert(
                            face_index[circulator->twin()->face()]);
                    }
                    else {
                        unbounded = true;
                        break;
                    }
                }
                while (++circulator != circulator_start);
                polygon_area = unbounded ?
                    INFINITY : polygon.area();
            }
            else {
                polygon_area = NAN;
            }
            particle_area.push_back(fabs(polygon_area));
        }
    }

    void
    append_quantile(std::set<std::vector<fastjet::PseudoJet>::
                    const_iterator> &jet_truncated,
                    std::vector<std::pair<
                    double, std::vector<fastjet::PseudoJet>::
                    const_iterator> > rho_vs_jet,
                    double quantile)
    {
        if (!rho_vs_jet.empty()) {
            std::sort(rho_vs_jet.begin(), rho_vs_jet.end());

            const size_t iterator_margin =
                floor(0.5 * (1 - quantile) * rho_vs_jet.size());

            for (std::vector<std::pair<
                     double, std::vector<fastjet::PseudoJet>::
                     const_iterator> >::const_iterator
                     iterator =
                     rho_vs_jet.begin() + iterator_margin;
                 iterator != rho_vs_jet.end() - iterator_margin;
                 iterator++) {
                jet_truncated.insert(iterator->second);
            }
        }
    }

    std::vector<double> ue_estimation_truncated_mean(
        std::vector<fastjet::PseudoJet> jet,
        size_t order_fourier = 3, double quantile = 0.5)
    {
        // Since the windows are staggered by 2x, there are 2 *
        // quantile * jet.size() / nwindow per window;

        const unsigned int nwindow_phi =
            std::max(1U, std::min(48U, static_cast<unsigned int>(
                floor(quantile * jet.size()))));
        const double window_width = 4 * M_PI / nwindow_phi;
        std::set<std::vector<fastjet::PseudoJet>::const_iterator>
            jet_truncated;

        for (size_t i = 0; i < nwindow_phi; i++) {
            const double phi_0 = i * (2 * M_PI / nwindow_phi) - M_PI;
            std::vector<std::pair<
                double, std::vector<fastjet::PseudoJet>::
                const_iterator> > rho_vs_jet;

            for (std::vector<fastjet::PseudoJet>::const_iterator
                     iterator = jet.begin();
                 iterator != jet.end(); iterator++) {
                if (angular_range_reduce(
                        iterator->phi_std() - phi_0) <
                    window_width) {
                    rho_vs_jet.push_back(
                        std::pair<double, std::vector<
                        fastjet::PseudoJet>::const_iterator>(
                            iterator->perp() / iterator->area(),
                            iterator));
                }
            }
            append_quantile(jet_truncated, rho_vs_jet, quantile);
        }

        std::vector<double> ret;

        if (!jet_truncated.empty()) {
            order_fourier =
                std::min(order_fourier,
                         (jet_truncated.size() - 1) / 2);

            TMatrixD a(jet_truncated.size(), 2 * order_fourier + 1);
            TVectorD b(jet_truncated.size());

            TMatrixDColumn(a, 0) = 1;

            size_t row = 0;

            for (std::set<std::vector<fastjet::PseudoJet>::
                     const_iterator>::const_iterator iterator =
                     jet_truncated.begin();
                 iterator != jet_truncated.end(); iterator++) {
                const double phi = (*iterator)->phi_std();

                for (size_t j = 0; j < order_fourier; j++) {
                    a(row, 2 * j + 1) = cos((j + 1) * phi);
                    a(row, 2 * j + 2) = sin((j + 1) * phi);
                }
                b(row) = (*iterator)->perp() / (*iterator)->area();
                row++;
            }

            TDecompSVD a_svd(a);
            Bool_t status;
            TVectorD x = a_svd.Solve(b, status);

            if (status != kFALSE) {
                for (size_t i = 0; i < 2 * order_fourier + 1; i++) {
                    ret.push_back(x(i));
                }
            }
        }

        return ret;
    }

    double evaluate_ue(std::vector<double> ue_estimate,
                       double azimuth, double r = 0)
    {
        if (ue_estimate.empty()) {
            return 0;
        }

        double s = ue_estimate[0];

        for (size_t i = 0; i < (ue_estimate.size() - 1) / 2; i++) {
            const double v =
                sqrt(std::pow(ue_estimate[2 * i + 2], 2) +
                     std::pow(ue_estimate[2 * i + 1], 2));
            const double psi = atan2(ue_estimate[2 * i + 2],
                                     ue_estimate[2 * i + 1]);
            const double k = i + 1;

            s += r == 0 ? 1.0 : (2 * j1(k * r) / (k * r)) *
                v * cos(k * azimuth - psi);
        }

        return s;
    }

    // fastjet::PseudoJet user indices -2 and -3 are used to tag the
    // EM particles/EMCAL clusters and muons. The index -1 is already
    // taken, being the fastjet::PseudoJet default initializer. After
    // the removal of EM and muons, -1 then implicitly means hadronic

    enum {
        USER_INDEX_DEFAULT_OR_TRACK = -1,
        USER_INDEX_EM               = -2,
        USER_INDEX_MUON             = -3
    };

    double jet_emf(const std::vector<fastjet::PseudoJet> constituent,
                   double scale_em_ghost = 1)
    {
        double sum_hadronic = 0;
        double sum_em = 0;

        for (std::vector<fastjet::PseudoJet>::const_iterator
                 iterator_constituent = constituent.begin();
             iterator_constituent != constituent.end();
             iterator_constituent++) {
            switch (iterator_constituent->user_index()) {
            case USER_INDEX_DEFAULT_OR_TRACK:
                sum_hadronic += iterator_constituent->perp();
                break;
            case USER_INDEX_EM:
                sum_em += iterator_constituent->perp();
                break;
            }
        }

        return sum_hadronic + sum_em > 0 ?
            sum_em / (sum_hadronic * scale_em_ghost + sum_em) :
            NAN;
    }

    size_t jet_multiplicity(const std::vector<fastjet::PseudoJet>
                            constituent)
    {
        size_t multiplicity = 0;

        for (std::vector<fastjet::PseudoJet>::const_iterator
                 iterator_constituent = constituent.begin();
             iterator_constituent != constituent.end();
             iterator_constituent++) {
            switch (iterator_constituent->user_index()) {
            case USER_INDEX_DEFAULT_OR_TRACK:
            case USER_INDEX_EM:
            case USER_INDEX_MUON:
                multiplicity++;
                break;
            }
        }

        return multiplicity;
    }

    double constituent_perp(const fastjet::PseudoJet constituent,
                            double scale_em_ghost = 1)
    {
        switch (constituent.user_index()) {
        case USER_INDEX_DEFAULT_OR_TRACK:
        case USER_INDEX_MUON:
            return constituent.perp();
        case USER_INDEX_EM:
            return constituent.perp() / scale_em_ghost;
        default:
            return 0;
        }
    }

    double jet_ptd(const std::vector<fastjet::PseudoJet> constituent,
                   double scale_em_ghost = 1)
    {
        double sum_1 = 0;
        double sum_2 = 0;

        for (std::vector<fastjet::PseudoJet>::const_iterator
                 iterator_constituent = constituent.begin();
             iterator_constituent != constituent.end();
             iterator_constituent++) {
            const double perp =
                constituent_perp(*iterator_constituent);

            sum_1 += perp;
            sum_2 += std::pow(perp, 2);
        }

        // The default value is the one particle limit

        return sum_1 > 0 ? sqrt(sum_2) / sum_1 : 1;
    }

    void jet_width_sigma(double sigma[],
                         const fastjet::PseudoJet jet,
                         const std::vector<fastjet::PseudoJet>
                         constituent,
                         double scale_em_ghost = 1)
    {
        double m11 = 0;
        double m22 = 0;
        double m12_m21 = 0;
        double sum_2 = 0;

        for (std::vector<fastjet::PseudoJet>::const_iterator
                 iterator_constituent = constituent.begin();
             iterator_constituent != constituent.end();
             iterator_constituent++) {
            const double perp_2 =
                std::pow(constituent_perp(*iterator_constituent), 2);
            const double dpseudorapidity =
                iterator_constituent->pseudorapidity() -
                jet.pseudorapidity();
            const double dazimuth = angular_range_reduce(
                iterator_constituent->phi_std() -
                jet.phi_std());

            m11 += perp_2 * std::pow(dpseudorapidity, 2);
            m22 += perp_2 * std::pow(dazimuth, 2);
            m12_m21 -= perp_2 * dpseudorapidity * dazimuth;
            sum_2 += perp_2;
        }

        // a = -1 in the equation for eigenvalues
        const double b = m11 + m22;
        const double c = std::pow(m12_m21, 2) - m11 * m22;
        const double q = -0.5 *
            (b + copysign(1, b) * sqrt(std::pow(b, 2) + 4 * c));

        // Major axis, x2
        sigma[0] = sqrt(c / (q * sum_2));
        // Minor axis, x1
        sigma[1] = sqrt(q / (-sum_2));
    }

    void jet_width_sigma_h(float sigma[],
                           const fastjet::PseudoJet jet,
                           const std::vector<fastjet::PseudoJet>
                           constituent,
                           double scale_em_ghost = 1)
    {
        double sigma_d[2];

        jet_width_sigma(sigma_d, jet, constituent, scale_em_ghost);
        for (size_t i = 0; i < 2; i++) {
            sigma[i] = half(sigma_d[i]);
        }
    }

}

#endif // JET_H_
