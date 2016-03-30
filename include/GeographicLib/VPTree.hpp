/**
 * \file VPTree.hpp
 * \brief Header for GeographicLib::VPTree class
 *
 * Copyright (c) Charles Karney (2016) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_VPTREE_HPP)
#define GEOGRAPHICLIB_VPTREE_HPP 1

#if !(__cplusplus >= 201103 || (defined(_MSC_VER) && _MSC_VER >= 1700))
#error "VPTree requires C++11"
#endif

#include <algorithm>            // for nth_element, minmax_element
#include <vector>
#include <queue>                // for priority_queue
#include <utility>              // for swap + pair
#include <limits>
#include <cmath>
#include <iostream>
#include <memory>               // for unique_ptr

#if defined(_MSC_VER)
// Squelch warnings about constant conditional expressions
#  pragma warning (push)
#  pragma warning (disable: 4127)
#endif

namespace GeographicLib {

  /**
   * \brief Vantage-point tree
   *
   * <b>WARNING</b>: this requires C++11
   *
   * This class implements the vantage-point tree described by
   * - J. K. Uhlmann,
   *   <a href="doi:10.1016/0020-0190(91)90074-r">
   *   Satisfying general proximity/similarity queries with metric trees</a>,
   *   Information Processing Letters 40 175&ndash;179 (1991).
   * - P. N. Yianilos,
   *   <a href="http://pnylab.com/pny/papers/vptree/vptree/">
   *   Data structures and algorithms for nearest neighbor search in general
   *   metric spaces</a>, Proc. 4th ACM-SIAM Symposium on Discrete Algorithms,
   *   (SIAM, 1993). pp. 311&ndash;321.
   *
   * Given a set of points in some space and a distance function \e d
   * satisfying the metric conditions:
   * \f[
   * \begin{align}
   *  d(x,y) &\ge 0,\\
   *  d(x,y) &= 0, \ \text{iff $x = y$},\\
   *  d(x,y) &= d(y,x),\\
   *  d(x,z) &\ge d(x,y) + d(y,z),
   * \end{align}
   * \f]
   * the vantage-point (VP) tree provides an efficient way of determining
   * nearest neighbors.  Typically the cost of constructing a VP tree of \e N
   * points is \e N log \e N, while the cost of a query is log \e N.  Thus a VP
   * tree should be used in situations where \e N is large and at least log \e
   * N queries are to be made.
   *
   * - This implementation includes Yianilos' upper and lower bounds for the
   *   inside and outside sets.  This helps limit the number of searches
   *   (compared to just using the median distance).
   * - Rather than do a depth-first or breath-first search on the tree, the
   *   nodes to be processed are put on a priority queue with the nodes most
   *   likely to contain close points processed first.  Frequently, this allows
   *   nodes lower down on the priority queue to be skipped.
   * - This technique also allows non-exhaustive searchs to be performed (to
   *   answer questions such as "are there any points within 1km of the query
   *   point?).
   * - When building the tree, the points furthest from the parent vantage
   *   point in both the inside and outside sets are selected as the children's
   *   vantage points.  This information is already available from the
   *   computation of the upper and lower bounds of the children.  This choice
   *   seems to lead to a reasonably optimized tree.
   * - The leaf nodes contain a bucket of points (instead of just a vantage
   *   point).
   *
   * @tparam real the type used for measuring distances; it can be a real or
   *   signed integer type.
   * @tparam position the type for specifying points.
   * @tparam distance the type for a function object which takes takes two
   *   positions as arguments and returns the distance (of type \e real).
   *
   * Example of use:
   * \include example-VPTree.cpp
   **********************************************************************/
  template <typename real, typename position, class distance>
  class VPTree {
    static const unsigned bucketsize = 4;
    static const unsigned marker = unsigned(-1);
    VPTree(const VPTree&);            // copy constructor not allowed
    VPTree& operator=(const VPTree&); // copy assignment not allowed
  public:

    /**
     * Constructor for VPTree
     *
     * @param[in] pts a vector of points to include in the tree; VPTree
     *   retains a const reference to this vector.
     * @param[in] dist the distance function object; VPTree retains a
     *   const reference to this object.
     *
     * The distances computed by \e dist must satisfy the standard metric
     * conditions.  If not, the results are undefined.  Neither the data in \e
     * pts nor the query points should contain NaNs because such data violates
     * the metric conditions.
     *
     * <b>CAUTION</b>: Do not alter \e pts during the lifetime of the VPTree.
     **********************************************************************/
    VPTree(const std::vector<position>& pts, const distance& dist)
      : _pts(pts)
      , _dist(dist)
    {
      static_assert(std::numeric_limits<real>::is_signed,
                    "real must be a signed type");
      _mc = 0; _sc = 0;
      _c0 = 0; _c1 = 0; _k = 0;
      _cmin = std::numeric_limits<unsigned>::max(); _cmax = 0;
      // the pair contains distance+id
      std::vector<std::pair<real, unsigned>> ids(_pts.size());
      for (unsigned k = unsigned(ids.size()); k--;)
        ids[k] = std::make_pair(real(0), k);
      _root = init(ids, 0, unsigned(ids.size()), unsigned(ids.size() / 2));
    }

    /**
     * Search the VPTree
     *
     * @param[in] query the query point.
     * @param[out] ind a vector of indices to the closest points found.
     * @param[in] k the number of points to search for (default = 1).
     * @param[in] maxdist only return points with distances of \e maxdist or
     *   less from \e query (default is the maximum \e real).
     * @param[in] mindist only return points with distances of more than
     *   \e mindist from \e query (default = &minus;1).
     * @param[in] exhaustive whether to do an exhaustive search (default true).
     * @param[in] tol the tolerance on the results (default 0).
     * @return the distance to the closest point found (-1 if no points are
     *   found).
     *
     * The indices returned in \e ind are sorted by distance.
     *
     * With \e exhaustive = true and \e tol = 0 (their default values), this
     * finds the indices of \e k closest neighbors to \e query whose distances
     * to \e query are in (\e mindist, \e maxdist].  If these parameters have
     * their default values, then the bounds have no effect.  If \e query is
     * one of the points in the tree, then set \e mindist = 0 to prevent this
     * point from being returned.
     *
     * If \e exhaustive = false, exit as soon as \e k results satisfying the
     * distance criteria are found.  If less than \e k results are returned
     * then the search was exhaustive even if \e exhaustive = false.
     *
     * If \e tol is positive, do an approximate search; in this case the
     * results are to be interpreted as follows: if the <i>k</i>'th distance is
     * \e dk, then all results with distances less than or equal \e dk &minus;
     * \e tol are correct; all others are suspect &mdash; there may be other
     * closer results with distances greater or equal to \e dk &minus; \e tol.
     * If less than \e k results are found, then the search is exact.
     *
     * \e mindist should be used to exclude a "small" neighborhood of the query
     * point (relative to the average spacing of the data).  If \e mindist is
     * large, the efficiency of the search deteriorates.
     *
     * \e pts may contain coincident points (i.e., the distance between them
     * vanishes).  These are treated as distinct.
     **********************************************************************/
    real search(const position& query,
                std::vector<int>& ind,
                unsigned k = 1,
                real maxdist = std::numeric_limits<real>::max(),
                real mindist = -1,
                bool exhaustive = true,
                real tol = 0) const {
      std::priority_queue<std::pair<real, unsigned>> results;
      if (maxdist > mindist) {
        struct task {
          Node* n;                // the node
          real d;                 // how far query is outside boundary of node
          // -1 if on boundary or inside
          task(Node* n, real d) : n(n), d(d) {}
          bool operator<(const task& o) const {
            // sort in reverse order to process smallest d first
            return d > o.d;
          }
        };

        // distance to the kth closest point so far
        real tau = maxdist;
        std::priority_queue<task> todo;
        todo.push(task(_root.get(), -1));
        unsigned c = 0;
        while (!todo.empty()) {
          Node* current = todo.top().n;
          if (!current) continue;
          real d = todo.top().d;
          todo.pop();
          real tau1 = tau - tol;
          // compare tau and d again since tau may have become smaller.
          if (!(tau1 >= d)) continue;
          real dist = 0;   // to suppress warning about uninitialized variable
          bool exitflag = false, leaf = current->index == marker;
          for (unsigned i = 0; i < (leaf ? bucketsize : 1); ++i) {
            unsigned index = leaf ? current->leaves[i] : current->index;
            if (index == marker) break;
            dist = _dist(_pts[index], query);
            ++c;

            if (dist > mindist && dist <= tau) {
              if (results.size() == k) results.pop();
              results.push(std::make_pair(dist, index));
              if (results.size() == k) {
                if (exhaustive)
                  tau = results.top().first;
                else {
                  exitflag = true;
                  break;
                }
                if (tau <= tol) {
                  exitflag = true;
                  break;
                }
              }
            }
          }
          if (exitflag) break;

          if (current->index == marker) continue;
          tau1 = tau - tol;
          for (unsigned n = 0; n < 2; ++n) {
            if (current->child[n] && dist + current->data.upper[n] >= mindist) {
              if (dist < current->data.lower[n]) {
                d = current->data.lower[n] - dist;
                if (tau1 >= d) todo.push(task(current->child[n].get(), d));
              } else if (dist > current->data.upper[n]) {
                d = dist - current->data.upper[n];
                if (tau1 >= d) todo.push(task(current->child[n].get(), d));
              } else
                todo.push(task(current->child[n].get(), -1));
            }
          }
        }
        ++_k;
        _c1 += c;
        double omc = _mc;
        _mc += (c - omc) / _k;
        _sc += (c - omc) * (c - _mc);
        if (c > _cmax) _cmax = c;
        if (c < _cmin) _cmin = c;
      }

      real d = -1;
      ind.resize(results.size());

      for (int i = int(ind.size()); i--;) {
        ind[i] = int(results.top().second);
        if (i == 0) d = results.top().first;
        results.pop();
      }
      return d;

    }

    /**
     * @return the total number of points.
     **********************************************************************/
    int numpoints() const { return int(_pts.size()); }
    /**
     * @return a reference to the vector of points.
     **********************************************************************/
    const std::vector<position>& points() const { return _pts; }
    /**
     * @param[in] i the index of the point.
     * @return a reference to the <i>i</i>'th point.
     **********************************************************************/
    const position& point(unsigned i) const { return _pts[i]; }

    /**
     * Report accumulated statistics on the searches so far.
     *
     * @param[in,out] os the stream to write to
     * @return a reference to the stream
     **********************************************************************/
    std::ostream& report(std::ostream& os) const {
      os << "set size " << _pts.size() << "\n"
         << "setup cost " << _c0 << "\n"
         << "searches " << _k << "\n"
         << "search cost (total mean sd min max) "
         << _c1 << " "
         << unsigned(std::floor(_mc + 0.5)) << " "
         << unsigned(std::floor(std::sqrt(_sc / (_k - 1)) + 0.5)) << " "
         << _cmin << " " << _cmax << "\n";
      return os;
    }

  private:
    const std::vector<position>& _pts;
    const distance& _dist;
    mutable double _mc, _sc;
    mutable unsigned _c0, _c1, _k, _cmin, _cmax;
    struct Node {
      unsigned index;
      struct bounds {
        real lower[2], upper[2];  // bounds on inner/outer distances
      };
      union {
        bounds data;
        unsigned leaves[bucketsize];
      };
      std::unique_ptr<Node> child[2];

      Node()
        : index(0)
      {}

    };
    std::unique_ptr<Node> _root;

    std::unique_ptr<Node> init(std::vector<std::pair<real, unsigned>>& ids,
                               unsigned l, unsigned u, unsigned vp) {
      if (u == l) {
        return std::unique_ptr<Node>();
      }

      std::unique_ptr<Node> node(new Node());

      if (u - l > (bucketsize == 0 ? 1 : bucketsize)) {

        // choose a vantage point and move it to the start
        unsigned i = vp;
        std::swap(ids[l], ids[i]);

        unsigned m = (u + l + 1) / 2;

        for (unsigned k = l + 1; k < u; ++k) {
          ids[k].first = _dist(_pts[ids[l].second], _pts[ids[k].second]);
          ++_c0;
        }
        // partitian around the median distance
        std::nth_element(ids.begin() + l + 1, ids.begin() + m, ids.begin() + u);
        node->index = ids[l].second;
        if (m > l + 1) { // node->child[0] is possibly empty
          auto t = std::minmax_element(ids.begin() + l + 1, ids.begin() + m);
          node->data.lower[0] = t.first->first;
          node->data.upper[0] = t.second->first;
          // Use point with max distance as vantage point; this point act as a
          // "corner" point and leads to a good partition.
          node->child[0] = init(ids, l + 1, m,
                                unsigned(t.second - ids.begin()));
        }
        auto t = std::max_element(ids.begin() + m, ids.begin() + u);
        node->data.lower[1] = ids[m].first;
        node->data.upper[1] = t->first;
        // Use point with max distance as vantage point here too
        node->child[1] = init(ids, m, u, unsigned(t - ids.begin()));
      } else {
        if (bucketsize == 0)
          node->index = ids[l].second;
        else {
          node->index = marker;
          for (unsigned i = l; i < u; ++i)
            node->leaves[i-l] = ids[i].second;
          for (unsigned k = u - l; k < bucketsize; ++k)
            node->leaves[k] = marker;
        }
      }

      return node;
    }
  };

} // namespace GeographicLib

#if defined(_MSC_VER)
#  pragma warning (pop)
#endif

#endif  // GEOGRAPHICLIB_VPTREE_HPP
