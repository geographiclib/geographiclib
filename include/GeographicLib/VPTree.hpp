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

#include <algorithm>            // for nth_element, max_element, etc.
#include <vector>
#include <queue>                // for priority_queue
#include <utility>              // for swap + pair
#include <limits>
#include <cmath>
#include <iostream>
#include <GeographicLib/Constants.hpp> // Only for GEOGRAPHICLIB_STATIC_ASSERT

// Not ready for boost serialization yet
#define HAVE_BOOST_SERIALIZATION 0
#if HAVE_BOOST_SERIALIZATION
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/array.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#endif

#if defined(_MSC_VER)
// Squelch warnings about constant conditional expressions
#  pragma warning (push)
#  pragma warning (disable: 4127)
#endif

namespace GeographicLib {

  /**
   * \brief Vantage-point tree
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
    static const int bucketsize = 10;
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
      GEOGRAPHICLIB_STATIC_ASSERT(std::numeric_limits<real>::is_signed,
                                  "real must be a signed type");
      GEOGRAPHICLIB_STATIC_ASSERT(bucketsize >= 0 && bucketsize <= 10,
                                  "bad bucketsize");
      _mc = 0; _sc = 0;
      _c0 = 0; _c1 = 0; _k = 0;
      _cmin = std::numeric_limits<int>::max(); _cmax = 0;
      // the pair contains distance+id
      std::vector<item> ids(_pts.size());
      for (int k = int(ids.size()); k--;)
        ids[k] = std::make_pair(real(0), k);
      init(ids, 0, int(ids.size()), int(ids.size() / 2));
      std::cout << sizeof(Node) << " " << _tree.size() * sizeof(Node) /
        double(_pts.size() * sizeof(position)) << "\n";
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
                int k = 1,
                real maxdist = std::numeric_limits<real>::max(),
                real mindist = -1,
                bool exhaustive = true,
                real tol = 0) const {
      std::priority_queue<item> results;
      if (k > 0 && maxdist > mindist) {
        // distance to the kth closest point so far
        real tau = maxdist;
        // first is negative of how far query is outside boundary of node
        // +1 if on boundary or inside
        // second is node index
        std::priority_queue<item> todo;
        todo.push(std::make_pair(real(1), int(_tree.size()) - 1));
        int c = 0;
        while (!todo.empty()) {
          int n = todo.top().second;
          real d = -todo.top().first;
          todo.pop();
          real tau1 = tau - tol;
          // compare tau and d again since tau may have become smaller.
          if (!(n >= 0 && tau1 >= d)) continue;
          const Node& current = _tree[n];
          real dist = 0;   // to suppress warning about uninitialized variable
          bool exitflag = false, leaf = current.index < 0;
          for (int i = 0; i < (leaf ? bucketsize : 1); ++i) {
            int index = leaf ? current.leaves[i] : current.index;
            if (index < 0) break;
            dist = _dist(_pts[index], query);
            ++c;

            if (dist > mindist && dist <= tau) {
              if (int(results.size()) == k) results.pop();
              results.push(std::make_pair(dist, index));
              if (int(results.size()) == k) {
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

          if (current.index < 0) continue;
          tau1 = tau - tol;
          for (int n = 0; n < 2; ++n) {
            if (current.data.child[n] >= 0 &&
                dist + current.data.upper[n] >= mindist) {
              if (dist < current.data.lower[n]) {
                d = current.data.lower[n] - dist;
                if (tau1 >= d)
                  todo.push(std::make_pair(-d, current.data.child[n]));
              } else if (dist > current.data.upper[n]) {
                d = dist - current.data.upper[n];
                if (tau1 >= d)
                  todo.push(std::make_pair(-d, current.data.child[n]));
              } else
                todo.push(std::make_pair(real(1), current.data.child[n]));
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
     *
     * \e i must lie in the range [0, numpoints()).  No checking is done on the
     * value of \e i.
     **********************************************************************/
    const position& point(int i) const { return _pts[i]; }

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
         << int(std::floor(_mc + 0.5)) << " "
         << int(std::floor(std::sqrt(_sc / (_k - 1)) + 0.5)) << " "
         << _cmin << " " << _cmax << "\n";
      return os;
    }

  private:
    // Package up a real and an int.  We will want to sort on the real so put
    // it first.
    typedef std::pair<real, int> item;
    const std::vector<position>& _pts;
    const distance& _dist;
    mutable double _mc, _sc;
    mutable int _c0, _c1, _k, _cmin, _cmax;

    class Node {
    public:
      struct bounds {
        real lower[2], upper[2];  // bounds on inner/outer distances
        int child[2];
      };
      union {
        bounds data;
        int leaves[bucketsize];
      };
      int index;

      Node()
        : index(-1)
      {
        for (int i = 0; i < 2; ++i) {
          data.lower[i] = data.upper[i] = 0;
          data.child[i] = -1;
        }
      }
    private:
#if HAVE_BOOST_SERIALIZATION
      friend class boost::serialization::access;
      template<class Archive> void save(Archive& ar, const unsigned int) const {
        ar & boost::serialization::make_nvp("index" , index );
        if (index < 0)
          ar & boost::serialization::make_nvp("leaves", leaves);
        else
          ar & boost::serialization::make_nvp("lower", data.lower)
            & boost::serialization::make_nvp("upper", data.upper)
            & boost::serialization::make_nvp("child", data.child);
      }
      template<class Archive> void load(Archive& ar, const unsigned int) {
        ar & boost::serialization::make_nvp("index" , index );
        if (index < 0)
          ar & boost::serialization::make_nvp("leaves", leaves);
        else
          ar & boost::serialization::make_nvp("lower", data.lower)
            & boost::serialization::make_nvp("upper", data.upper)
            & boost::serialization::make_nvp("child", data.child);
      }
      template<class Archive>
      void serialize(Archive &ar, const unsigned int file_version)
      { boost::serialization::split_member(ar, *this, file_version); }
#endif
    };

    std::vector<Node> _tree;
    int init(std::vector<item>& ids, int l, int u, int vp) {

      if (u == l)
        return -1;
      Node node;

      if (u - l > (bucketsize == 0 ? 1 : bucketsize)) {

        // choose a vantage point and move it to the start
        int i = vp;
        std::swap(ids[l], ids[i]);

        int m = (u + l + 1) / 2;

        for (int k = l + 1; k < u; ++k) {
          ids[k].first = _dist(_pts[ids[l].second], _pts[ids[k].second]);
          ++_c0;
        }
        // partitian around the median distance
        std::nth_element(ids.begin() + l + 1, ids.begin() + m, ids.begin() + u);
        node.index = ids[l].second;
        if (m > l + 1) { // node.child[0] is possibly empty
          typename std::vector<item>::iterator
            t = std::min_element(ids.begin() + l + 1, ids.begin() + m);
          node.data.lower[0] = t->first;
          t = std::max_element(ids.begin() + l + 1, ids.begin() + m);
          node.data.upper[0] = t->first;
          // Use point with max distance as vantage point; this point act as a
          // "corner" point and leads to a good partition.
          node.data.child[0] = init(ids, l + 1, m, int(t - ids.begin()));
        }
        typename std::vector<item>::iterator
          t = std::max_element(ids.begin() + m, ids.begin() + u);
        node.data.lower[1] = ids[m].first;
        node.data.upper[1] = t->first;
        // Use point with max distance as vantage point here too
        node.data.child[1] = init(ids, m, u, int(t - ids.begin()));
      } else {
        if (bucketsize == 0)
          node.index = ids[l].second;
        else {
          node.index = -1;
          // Sort the bucket entries so that the tree is independent of the
          // implementation of nth_element.
          std::sort(ids.begin() + l, ids.begin() + u);
          for (int i = l; i < u; ++i)
            node.leaves[i-l] = ids[i].second;
          for (int i = u - l; i < bucketsize; ++i)
            node.leaves[i] = -1;
        }
      }

      _tree.push_back(node);
      return int(_tree.size()) - 1;
    }
  public:
  };

} // namespace GeographicLib

#if defined(_MSC_VER)
#  pragma warning (pop)
#endif

#endif  // GEOGRAPHICLIB_VPTREE_HPP
