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
   * <b>WARNING</b>: this needs C++11
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
   * the VP tree provides an efficient way of determining nearest neighbors.
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
   * - The leaf nodes can contain a bucket of points (instead of just a vantage
   *   point).
   *
   * @tparam distt the type used for measuring distances; it can be a real or
   *   signed integer type.
   * @tparam position the type for specifying points.
   * @tparam distanceop the type for a function object which takes takes two
   *   positions as arguments and returns the distance (of type distt).
   *
   * Example of use:
   * \include example-VPTree.cpp
   **********************************************************************/
  template <typename distt, typename position, class distanceop>
  class VPTree {
    static const unsigned bucketsize = 4;
    static const unsigned marker = unsigned(-1);
    VPTree(const VPTree&);            // copy constructor not allowed
    VPTree& operator=(const VPTree&); // copy assignment not allowed
  public:

    /**
     * Constructor for VPTree
     *
     * @param[in] items a vector of points to include in the tree; VPTree
     *   retains a const reference to this vector.
     * @param[in] distance the distance function object; VPTree retains a
     *   const reference to this object.
     *
     * The distances computed by \e distance must satisfy the standard metric
     * conditions.  If not, the results are undefined.
     *
     * <b>CAUTION</b>: Do not alter \e items after the VPTree has been
     * constructed.
     **********************************************************************/
    VPTree(const std::vector<position>& items,
           const distanceop& distance)
      : _items(items)
      , _distance(distance)
    {
      static_assert(std::numeric_limits<distt>::is_signed,
                    "distt must be a signed type");
      _mc = 0; _sc = 0;
      _c0 = 0; _c1 = 0; _k = 0;
      _cmin = std::numeric_limits<unsigned>::max(); _cmax = 0;
      // the pair contains distance+id
      std::vector<std::pair<distt, unsigned>> ids(_items.size());
      for (unsigned k = unsigned(ids.size()); k--;)
        ids[k] = std::make_pair(distt(0), k);
      _root = init(ids, 0, unsigned(ids.size()), unsigned(ids.size() / 2));
    }

    /**
     * Search the VPTree
     *
     * @param[in] query the query point.
     * @param[out] ind a vector of indices to the closest points found.
     * @param[in] k the number of points to search for (default = 1).
     * @param[in] maxdist only return points with distances of \e maxdist or
     *   less from \e query (default is the maximum distt).
     * @param[in] mindist only return points with distances of more than
     *   \e mindist from \e query (default = &minus;1).
     * @param[in] exhaustive if true (the default) perform an exhaustive
     *   search.  If false, exit as soon as \e k results satisfying the
     *   distance criteria are found.  If less than \e k results are returned
     *   then the search was exhaustive even if \e exhaustive = false.
     * @param[in] tol If 0 (the default), then do an exact search.  If
     *   positive, do an approximate search; in this case the results are to be
     *   interpreted as follows: if the <i>k</i>'th distance is \e dk, then all
     *   results with distances less than or equal \e dk &minus; \e tol are
     *   correct; all others are suspect &mdash; there may be other closer
     *   results with distances greater or equal to \e dk &minus; \e tol.  If
     *   less than \e k results are found, then the search is exact.
     *
     * The distance results are in (\e mindist, \e maxdist].  If these
     * parameters have their default values, then the bounds have no effect.
     **********************************************************************/
    void search(const position& query,
                std::vector<unsigned>& ind,
                unsigned k = 1,
                distt maxdist = std::numeric_limits<distt>::max(),
                distt mindist = distt(-1),
                bool exhaustive = true,
                distt tol = distt(0)) const {
      std::priority_queue<std::pair<distt, unsigned>> results;
      if (maxdist > mindist) {
        struct task {
          Node* n;                // the node
          distt d;                // how far query is outside boundary of node
          // 0 if on boundary or inside
          task(Node* n, distt d) : n(n), d(d) {}
          bool operator<(const task& o) const {
            // sort in reverse order to process smallest d first
            return d > o.d;
          }
        };

        // distance to the kth closest point so far
        distt tau = maxdist;
        std::priority_queue<task> todo;
        todo.push(task(_root.get(), 0));
        unsigned c = 0;
        while (!todo.empty()) {
          Node* current = todo.top().n;
          if (!current) continue;
          distt d = todo.top().d;
          todo.pop();
          distt tau1 = tau - tol;
          // compare tau and d again since tau may have become smaller.
          if (!(tau1 >= d)) continue;
          distt dist = 0;   // to suppress warning about uninitialized variable
          bool exitflag = false, leaf = current->index == marker;
          for (unsigned i = 0; i < (leaf ? bucketsize : 1); ++i) {
            unsigned index = leaf ? current->leaves[i] : current->index;
            if (index == marker) break;
            dist = _distance(_items[index], query);
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
          if (current->inside) {
            if (dist < current->data.inl) {
              d = current->data.inl - dist;
              if (tau1 >= d) todo.push(task(current->inside.get(), d));
            } else if (dist > current->data.inu) {
              d = dist - current->data.inu;
              if (tau1 >= d) todo.push(task(current->inside.get(), d));
            } else
              todo.push(task(current->inside.get(), distt(0)));
          }
          if (current->outside) {
            if (dist < current->data.outl) {
              d = current->data.outl - dist;
              if (tau1 >= d) todo.push(task(current->outside.get(), d));
            } else if (dist > current->data.outu) {
              d = dist - current->data.outu;
              if (tau1 >= d) todo.push(task(current->outside.get(), d));
            } else
              todo.push(task(current->outside.get(), distt(0)));
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

      ind.resize(results.size());

      for (unsigned i = unsigned(ind.size()); i--;) {
        ind[i] =results.top().second;
        results.pop();
      }

    }

    /**
     * @return the total number of points.
     **********************************************************************/
    unsigned numpoints() const { return _items.size(); }
    /**
     * @return a reference to the vector of points.
     **********************************************************************/
    const std::vector<position>& points() const { return _items; }
    /**
     * @param[in] i the index of the point.
     * @return a reference to the <i>i</i>'th point.
     **********************************************************************/
    const position& point(unsigned i) const { return _items[i]; }

    /**
     * Report acculumated statistics on the searches so far.
     *
     * @param[in,out] os the stream to write to
     * @return a reference to the stream
     **********************************************************************/
    std::ostream& report(std::ostream& os) const {
      os << "set size " << _items.size() << "\n"
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
    const std::vector<position>& _items;
    const distanceop& _distance;
    mutable double _mc, _sc;
    mutable unsigned _c0, _c1, _k, _cmin, _cmax;
    struct Node {
      unsigned index;
      struct bounds {
        distt inl, inu, outl, outu;  // bounds on inner/outer distances
      };
      union {
        bounds data;
        unsigned leaves[bucketsize];
      };
      std::unique_ptr<Node> inside;
      std::unique_ptr<Node> outside;

      Node()
        : index(0)
      {}

    };
    std::unique_ptr<Node> _root;

    std::unique_ptr<Node> init(std::vector<std::pair<distt, unsigned>>& ids,
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
          ids[k].first = _distance(_items[ids[l].second],
                                   _items[ids[k].second]);
          ++_c0;
        }
        // partitian around the median distance
        std::nth_element(ids.begin() + l + 1, ids.begin() + m, ids.begin() + u);
        node->index = ids[l].second;
        if (m > l + 1) { // node->inside is possibly empty
          auto t = std::minmax_element(ids.begin() + l + 1, ids.begin() + m);
          node->data.inl = t.first->first;
          node->data.inu = t.second->first;
          // Use point with max distance as vantage point; this point act as a
          // "corner" point and leads to a good partition.
          node->inside = init(ids, l + 1, m, unsigned(t.second - ids.begin()));
        }
        auto t = std::max_element(ids.begin() + m, ids.begin() + u);
        node->data.outl = ids[m].first;
        node->data.outu = t->first;
        // Use point with max distance as vantage point here too
        node->outside = init(ids, m, u, unsigned(t - ids.begin()));
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
