#ifndef UTILS_HPP
#define UTILS_HPP

#include <math.h>
#include "diy/link.hpp"
#include <cstdio>
#include "misc.h"

// This utility is the same as diy's pick.hpp, but ensures that distance computation is
// done in double precision even though the bounds are integer
//
// This is a safety measure to reduce the chance of numerical errors and should be used
namespace utl{

    // Find the distance between point `p` and box `bounds`
    template<class Point, class Bounds>
        double
        distance(const Bounds& bounds, const Point& p)
        {
            // fmt::print(stderr,"bnds [{} {}] [{} {}] [{} {}], ({} {} {})\n", bounds.min[0], bounds.max[0], bounds.min[1], bounds.max[1], bounds.min[2], bounds.max[2], p[0], p[1], p[2]);
            double res = 0;
            for (int i = 0; i < p.size(); ++i)
            {
                // avoids all the annoying case logic by finding
                // diff = max(bounds.min[i] - p[i], 0, p[i] - bounds.max[i])
                double diff = 0, d;

                d = (double)(bounds.min[i]) - (double)(p[i]);
                if (d > diff) diff = d;
                d = (double)(p[i]) - (double)(bounds.max[i]);
                if (d > diff) diff = d;

                // TP, 10/10/19: I reverted this back to match diy::distance()
//                 if (bounds.max[i] == p[i])
//                     res++;

                res += diff*diff;
            }
            return sqrt(res);
        }


    // Finds the neighbor(s) containing the target point
    template<class Bounds, class Point, class OutIter>
        void
        in(
            const diy::RegularLink<Bounds>& link,   // neighbors
            const Point&                    p,      // target point
            OutIter                         out,    // insert iterator for output set of neighbors
            const Bounds&                   domain, // global domain bounds
            bool                            core)   // check against core (or bounds, if false)
        {
            Bounds neigh_bounds {0}; // neighbor block bounds

            // for all neighbors of this block
            for (int n = 0; n < link.size(); n++)
            {
                if (core)
                    neigh_bounds = link.core(n);
                else
                    neigh_bounds = link.bounds(n);

                // wrap neighbor bounds, if necessary, otherwise bounds will be unchanged
                wrap_bounds(neigh_bounds, link.wrap(n), domain);
                // dprint("nbr: %d dist %f", n, utl::distance(neigh_bounds, p));
                double dist = utl::distance(neigh_bounds, p);
                // dprint("nbrid %d, %f", n, dist);
                if (dist == 0){
                    *out++ = n;
                    // dprint("infn %d %d", link.size(), n);
                }
                    
            } // for all neighbors
        }
}

#endif
