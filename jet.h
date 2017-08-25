// -*- mode: c++; -*-

#ifndef JET_H_
#define JET_H_

#include <vector>
#include <set>
#include <map>

#include <TDecompSVD.h>
#include <TPolyLine.h>

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

#include <special_function.h>

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

    void voronoi_insert_alice_tpc(
        voronoi_diagram_t &diagram,
        std::map<voronoi_diagram_t::Face_handle, size_t> &face_index,
        const std::vector<point_2d_t> particle_pseudorapidity_azimuth)
    {
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
                        p(iterator->x() * (1 - 2 * (j & 1)) + j * (2 * 0.9),
                          iterator->y() + k * (2 * M_PI));
                    const voronoi_diagram_t::Face_handle
                        handle = diagram.insert(p);

                    face_index[handle] = iterator -
                        particle_pseudorapidity_azimuth.begin();
                }
            }
        }
    }

    void voronoi_area_incident(
        std::vector<double> &particle_area,
        std::vector<std::set<size_t> > &particle_incident,
        const std::vector<point_2d_t> particle_pseudorapidity_azimuth)
    {
        voronoi_diagram_t diagram;
        std::map<voronoi_diagram_t::Face_handle, size_t> face_index;

        voronoi_insert_alice_tpc(diagram, face_index,
                                 particle_pseudorapidity_azimuth);

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

    void voronoi_polygon(
        std::vector<TPolyLine> &polyline,
        const std::vector<point_2d_t> &
        particle_pseudorapidity_azimuth)
    {
        voronoi_diagram_t diagram;
        std::map<voronoi_diagram_t::Face_handle, size_t> face_index;

        voronoi_insert_alice_tpc(diagram, face_index,
                                 particle_pseudorapidity_azimuth);

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
            std::vector<double> x;
            std::vector<double> y;

            if (face != NULL) {
                voronoi_diagram_t::Ccb_halfedge_circulator
                    circulator_start = (*face)->outer_ccb();

                voronoi_diagram_t::Ccb_halfedge_circulator
                    circulator = circulator_start;

                // Circle around the edges and extract the polygon
                // vertices
                do {
                    if (circulator->has_target()) {
                        x.push_back(circulator->target()->point().x());
                        y.push_back(circulator->target()->point().y());
                    }
                }
                while (++circulator != circulator_start);
            }
            if (!x.empty()) {
                x.push_back(x.front());
                y.push_back(y.front());
            }
            polyline.push_back(TPolyLine(x.size(), &x[0], &y[0]));
        }
    }

    void
    append_quantile(std::vector<fastjet::PseudoJet> &
                    constituent_truncated,
                    std::set<int> &constituent_truncated_user_index,
                    const std::vector<std::pair<
                    double, std::vector<fastjet::PseudoJet>::
                    const_iterator> > &rho_vs_jet_unsorted,
                    const fastjet::ClusterSequenceArea
                    cluster_sequence,
                    const std::vector<double> &particle_area,
                    double quantile)
    {
        if (rho_vs_jet_unsorted.empty()) {
            return;
        }

        std::vector<std::pair<
            double, std::vector<fastjet::PseudoJet>::
            const_iterator> > rho_vs_jet = rho_vs_jet_unsorted;

        std::sort(rho_vs_jet.begin(), rho_vs_jet.end());

        const size_t iterator_margin =
            floor(0.5 * (1 - quantile) * rho_vs_jet.size());

        for (std::vector<std::pair<
                 double, std::vector<fastjet::PseudoJet>::
                 const_iterator> >::const_iterator
                 iterator_rho_vs_jet =
                 rho_vs_jet.begin() + iterator_margin;
             iterator_rho_vs_jet !=
                 rho_vs_jet.end() - iterator_margin;
             iterator_rho_vs_jet++) {
            std::vector<fastjet::PseudoJet> constituent =
                cluster_sequence.
                constituents(*iterator_rho_vs_jet->second);

            for (std::vector<fastjet::PseudoJet>::const_iterator
                     iterator_constituent = constituent.begin();
                 iterator_constituent != constituent.end();
                 iterator_constituent++) {
                const int index = iterator_constituent->user_index();

                if (index >= 0 && static_cast<size_t>(index) <
                    particle_area.size() &&
                    std::isfinite(particle_area[index]) &&
                    constituent_truncated_user_index.find(index) ==
                    constituent_truncated_user_index.end()) {
                    constituent_truncated.push_back(
                        *iterator_constituent);
                    constituent_truncated_user_index.insert(index);
                }
            }
        }
    }

    void constituent_quantile(
        std::vector<fastjet::PseudoJet> &constituent_truncated,
        std::set<int> &constituent_truncated_user_index,
        fastjet::ClusterSequenceArea cluster_sequence,
        std::vector<double> particle_area, double quantile)
    {
        const std::vector<fastjet::PseudoJet> jet =
            cluster_sequence.inclusive_jets(0);

        static const unsigned int nwindow_azimuth_max = 24U;

        // Since the windows are staggered by 2x, there are 2 *
        // quantile * jet.size() / nwindow per window;

        const unsigned int nwindow_azimuth =
            std::max(1U, std::min(
                nwindow_azimuth_max, static_cast<unsigned int>(
                    floor(quantile * jet.size()))));
        const double azimuth_window_width =
            2 * M_PI / nwindow_azimuth;

        for (size_t i = 0; i < nwindow_azimuth; i++) {
            const double azimuth_window_center =
                i * (2 * M_PI / nwindow_azimuth) - M_PI;
            std::vector<std::pair<
                double, std::vector<fastjet::PseudoJet>::
                const_iterator> > rho_vs_jet;

            for (std::vector<fastjet::PseudoJet>::const_iterator
                     iterator = jet.begin();
                 iterator != jet.end(); iterator++) {
                if (angular_range_reduce(
                        iterator->phi_std() - azimuth_window_center) <
                    azimuth_window_width) {
                    rho_vs_jet.push_back(
                        std::pair<double, std::vector<
                        fastjet::PseudoJet>::const_iterator>(
                            iterator->perp() / iterator->area(),
                            iterator));
                }
            }
            append_quantile(constituent_truncated,
                            constituent_truncated_user_index,
                            rho_vs_jet, cluster_sequence,
                            particle_area, quantile);
        }
    }

    std::pair<std::vector<double>, std::vector<double> >
    ue_estimation_truncated_mean(
        fastjet::ClusterSequenceArea cluster_sequence,
        std::vector<double> particle_area,
        size_t order_pseudorapidity_chebyshev = 4,
        size_t order_azimuth_fourier = 3,
        double quantile = 0.5)
    {
        std::vector<fastjet::PseudoJet> constituent_truncated;
        std::set<int> constituent_truncated_user_index;

        constituent_quantile(constituent_truncated,
                             constituent_truncated_user_index,
                             cluster_sequence, particle_area,
                             quantile);

        std::vector<double> pseudorapidity_dependence;

        if (!constituent_truncated.empty()) {
            order_pseudorapidity_chebyshev =
                std::min(order_pseudorapidity_chebyshev,
                         constituent_truncated.size() - 1);

            TMatrixD a(constituent_truncated.size(),
                       order_pseudorapidity_chebyshev + 1);
            TVectorD b(constituent_truncated.size());
            size_t row = 0;

            for (std::vector<fastjet::PseudoJet>::const_iterator
                     iterator = constituent_truncated.begin();
                 iterator != constituent_truncated.end();
                 iterator++) {
                const double area =
                    particle_area[iterator->user_index()];
                // The convenience of ALICE central tracks being from
                // pseudorapidity -0.9 to 0.9 (close to -1 to 1) is
                // taken advantage to avoid a linear transform for the
                // Chebyshev polynomials
                const double x = iterator->pseudorapidity();

                a(row, 0) = area;
                if (order_pseudorapidity_chebyshev >= 1) {
                    a(row, 1) = x * area;
                }

                // t[0] is T_n(x), t[1] is T_{n - 1}(x)
                double t[2] = { x, 1 };

                for (size_t j = 2;
                     j < order_pseudorapidity_chebyshev + 1; j++) {
                    const double tn1 = 2 * x * t[0] - t[1];

                    a(row, j) = tn1 * area;
                    t[0] = tn1;
                    t[1] = t[0];
                }
                b(row) = iterator->perp();
                row++;
            }

            TDecompSVD a_svd(a);
            Bool_t status;
            TVectorD x = a_svd.Solve(b, status);

            if (status != kFALSE) {
                for (size_t i = 0;
                     i < order_pseudorapidity_chebyshev + 1; i++) {
                    pseudorapidity_dependence.push_back(x(i));
                }
            }
        }

        std::vector<double> azimuth_dependence;

        if (!constituent_truncated.empty()) {
            order_azimuth_fourier =
                std::min(order_azimuth_fourier,
                         (constituent_truncated.size() - 1) / 2);

            TMatrixD a(constituent_truncated.size(),
                       2 * order_azimuth_fourier + 1);
            TVectorD b(constituent_truncated.size());
            size_t row = 0;

            for (std::vector<fastjet::PseudoJet>::const_iterator
                     iterator = constituent_truncated.begin();
                 iterator != constituent_truncated.end();
                 iterator++) {
                const double azimuth = iterator->phi_std();
                const double area =
                    particle_area[iterator->user_index()];

                a(row, 0) = area;
                for (size_t j = 0; j < order_azimuth_fourier; j++) {
                    a(row, 2 * j + 1) = cos((j + 1) * azimuth) * area;
                    a(row, 2 * j + 2) = sin((j + 1) * azimuth) * area;
                }
                b(row) = iterator->perp();
                row++;
            }

            TDecompSVD a_svd(a);
            Bool_t status;
            TVectorD x = a_svd.Solve(b, status);

            if (status != kFALSE) {
                for (size_t i = 0; i < 2 * order_azimuth_fourier + 1;
                     i++) {
                    azimuth_dependence.push_back(x(i));
                }
            }
        }

        return std::pair<
            std::vector<double>, std::vector<double> >(
                pseudorapidity_dependence, azimuth_dependence);
    }

    std::set<int> ue_user_index_truncated_mean(
        fastjet::ClusterSequenceArea cluster_sequence,
        std::vector<double> particle_area,
        size_t order_pseudorapidity_chebyshev = 4,
        size_t order_azimuth_fourier = 3,
        double quantile = 0.5)
    {
        std::vector<fastjet::PseudoJet> constituent_truncated;
        std::set<int> constituent_truncated_user_index;

        constituent_quantile(constituent_truncated,
                             constituent_truncated_user_index,
                             cluster_sequence, particle_area,
                             quantile);

        return constituent_truncated_user_index;
    }

    double evaluate_ue(std::pair<std::vector<double>,
                       std::vector<double> > ue_estimate,
                       double pseudorapidity, double azimuth)
    {
        if (ue_estimate.first.empty() ||
            ue_estimate.second.empty()) {
            return 0;
        }

        double p = ue_estimate.first[0];

        if (ue_estimate.first.size() >= 1) {
            p += ue_estimate.first[1] * pseudorapidity;
        }

        const double x = pseudorapidity;
        // t[0] is T_n(x), t[1] is T_{n - 1}(x)
        double t[2] = { x, 1 };

        for (size_t i = 2; i < ue_estimate.first.size(); i++) {
            const double tn1 = 2 * x * t[0] - t[1];

            p += ue_estimate.first[i] * tn1;
            t[0] = tn1;
            t[1] = t[0];
        }

        double a = ue_estimate.second[0];

        for (size_t i = 0; i < (ue_estimate.second.size() - 1) / 2; i++) {
            const double v =
                sqrt(std::pow(ue_estimate.second[2 * i + 2], 2) +
                     std::pow(ue_estimate.second[2 * i + 1], 2));
            const double psi = atan2(ue_estimate.second[2 * i + 2],
                                     ue_estimate.second[2 * i + 1]);
            const double k = i + 1;

            a += v * cos(k * azimuth - psi);
        }

        return p * a / ue_estimate.second[0];
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
                constituent_perp(*iterator_constituent,
                                 scale_em_ghost);

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
                std::pow(constituent_perp(*iterator_constituent,
                                          scale_em_ghost), 2);
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

}

#endif // JET_H_
