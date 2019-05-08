#ifndef UTILS_HPP
#define UTILS_HPP

#include <math.h>
#include "diy/link.hpp"
#include <cstdio>

namespace utl{
//! Find the distance between point `p` and box `bounds`.
template<class Point, class Bounds>
    double
    distance(const Bounds& bounds, const Point& p)
    {
        double res = 0;
        for (int i = 0; i < p.size(); ++i)
        {
        // avoids all the annoying case logic by finding
        // diff = max(bounds.min[i] - p[i], 0, p[i] - bounds.max[i])
            double diff = 0, d;

            d = bounds.min[i] - p[i];
            if (d > diff) diff = d;
            d = p[i] - bounds.max[i];
            if (d > diff) diff = d;

            if (bounds.max[i]==p[i]){
                res++;
            }

            res += diff*diff;
        }
        // res++;
        return sqrt(res);
    }


    //! Finds the neighbor(s) containing the target point.
    template<class Bounds, class Point, class OutIter>
    void

    in(const diy::RegularLink<Bounds>& link,  //!< neighbors
       const Point& p,                   //!< target point
       OutIter out,                      //!< insert iterator for output set of neighbors
       const Bounds& domain, 
       int flag)             //!< global domain bounds
    {
      Bounds neigh_bounds {0}; // neighbor block bounds

      // if (flag==1){
      //       fprintf(stderr, "link bounds (%d %d) (%d %d) (%d %d) ,,,", link.bounds(0).min[0], link.bounds(0).max[0], link.bounds(0).min[1], link.bounds(0).max[1], link.bounds(0).min[2], link.bounds(0).max[2]);
      //       fprintf(stderr, "link bounds (%d %d) (%d %d) (%d %d) ,,,", link.bounds(1).min[0], link.bounds(1).max[0], link.bounds(1).min[1], link.bounds(1).max[1], link.bounds(1).min[2], link.bounds(1).max[2]);
      //       fprintf(stderr, "link bounds (%d %d) (%d %d) (%d %d) ,,,", link.bounds(2).min[0], link.bounds(2).max[0], link.bounds(2).min[1], link.bounds(2).max[1], link.bounds(2).min[2], link.bounds(2).max[2]);

      // }

      // for all neighbors of this block
      for (int n = 0; n < link.size(); n++)
      {
        // wrap neighbor bounds, if necessary, otherwise bounds will be unchanged
        neigh_bounds = link.bounds(n);
        wrap_bounds(neigh_bounds, link.wrap(n), domain);

        if (utl::distance(neigh_bounds, p) == 0)
            *out++ = n;
      } // for all neighbors
  }

}

#endif 