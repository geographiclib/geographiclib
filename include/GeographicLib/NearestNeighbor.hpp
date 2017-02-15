/**
 * \file NearestNeighbor.hpp
 * \brief Header for GeographicLib::NearestNeighbor class
 *
 * Copyright (c) Charles Karney (2016) <charles@karney.com> and licensed under
 * the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_NEARESTNEIGHBOR_HPP)
#define GEOGRAPHICLIB_NEARESTNEIGHBOR_HPP 1

#include <algorithm>            // for nth_element, max_element, etc.
#include <vector>
#include <queue>                // for priority_queue
#include <utility>              // for swap + pair
#include <cstring>
#include <limits>
#include <cmath>
#include <iostream>
#include <sstream>
// Only for GEOGRAPHICLIB_STATIC_ASSERT and GeographicLib::GeographicErr
#include <GeographicLib/Constants.hpp>

#if !defined(GEOGRAPHICLIB_HAVE_BOOST_SERIALIZATION)
#define GEOGRAPHICLIB_HAVE_BOOST_SERIALIZATION 0
#endif
#if GEOGRAPHICLIB_HAVE_BOOST_SERIALIZATION
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>
#endif

#if defined(_MSC_VER)
// Squelch warnings about constant conditional expressions
#  pragma warning (push)
#  pragma warning (disable: 4127)
#endif

namespace GeographicLib {

  /**
   * \brief Nearest-neighbor calculations
   *
   * This class implements nearest-neighbor calculations using the
   * vantage-point tree described by
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
   * N queries are to be made.  The condition, \e N is large, means that
   * \f$ N \gg 2^D \f$, where \e D is the dimensionality of the space.
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
   * This class is templated so that it can handle arbitrary metric spaces as
   * follows:
   *
   * @tparam real the type used for measuring distances; it can be a real or
   *   signed integer type.
   * @tparam position the type for specifying points.
   * @tparam distance the type for a function object which takes takes two
   *   positions as arguments and returns the distance (of type \e real).
   *
   * The \e real type must support numeric_limits queries (specifically:
   * is_signed, is_integer, max(), digits, and digits10).
   *
   * Example of use:
   * \include example-NearestNeighbor.cpp
   **********************************************************************/
  template <typename real, typename position, class distance>
  class NearestNeighbor {
    // For tracking changes to the I/O format
    static const int version = 1;
    // This is what we get "free"; but if sizeof(real) = 1 (unlikely), allow 4
    // slots (and this accommodates the default value bucket = 4).
    static const int maxbucket = (2 + ((4 * sizeof(real)) / sizeof(int) >= 2 ?
                                       (4 * sizeof(real)) / sizeof(int) : 2));
  public:

    /**
     * Default constructor for NearestNeighbor.
     *
     * This is equivalent to specifying an empty set of points.
     **********************************************************************/
    NearestNeighbor() : _numpoints(0) {}

    /**
     * Constructor for NearestNeighbor.
     *
     * @param[in] pts a vector of points to include in the set.
     * @param[in] dist the distance function object.
     * @param[in] bucket the size of the buckets at the leaf nodes; this must
     *   lie in [0, 2 + 4*sizeof(real)/sizeof(int)] (default 4).
     * @exception GeographicErr if the value of \e bucket is out of bounds or
     *   the size of \e pts is too big for an int.
     * @exception std::bad_alloc if memory for the tree can't be allocated.
     *
     * The distances computed by \e dist must satisfy the standard metric
     * conditions.  If not, the results are undefined.  Neither the data in \e
     * pts nor the query points should contain NaNs or infinities because such
     * data violates the metric conditions.
     *
     * \e pts may contain coincident points (i.e., the distance between them
     * vanishes); these are treated as distinct.
     *
     * The choice of \e bucket is a tradeoff between space and efficiency.  A
     * larger \e bucket decreases the size of the NearestNeigbor object which
     * scales as pts.size() / max(1, bucket) and reduces the number of distance
     * calculations to construct the object by log2(bucket) * pts.size().
     * However each search then requires about bucket additional distance
     * calculations.
     *
     * \warning The same arguments \e pts and \e dist must be provided
     * to the Search() function.
     **********************************************************************/
    NearestNeighbor(const std::vector<position>& pts, const distance& dist,
                    int bucket = 4) {
      Initialize(pts, dist, bucket);
    }

    /**
     * Initialize or re-initialize NearestNeighbor.
     *
     * @param[in] pts a vector of points to include in the tree.
     * @param[in] dist the distance function object.
     * @param[in] bucket the size of the buckets at the leaf nodes; this must
     *   lie in [0, 2 + 4*sizeof(real)/sizeof(int)] (default 4).
     * @exception GeographicErr if the value of \e bucket is out of bounds or
     *   the size of \e pts is too big for an int.
     * @exception std::bad_alloc if memory for the tree can't be allocated.
     *
     * See also the documentation on the constructor.
     *
     * If an exception is thrown, the state of the NearestNeighbor is
     * unchanged.
     **********************************************************************/
    void Initialize(const std::vector<position>& pts, const distance& dist,
                    int bucket = 4) {
      GEOGRAPHICLIB_STATIC_ASSERT(std::numeric_limits<real>::is_signed,
                                  "real must be a signed type");
      if (!( 0 <= bucket && bucket <= maxbucket ))
        throw GeographicLib::GeographicErr
          ("bucket must lie in [0, 2 + 4*sizeof(real)/sizeof(int)]");
      if (pts.size() > size_t(std::numeric_limits<int>::max()))
        throw GeographicLib::GeographicErr("pts array too big");
      // the pair contains distance+id
      std::vector<item> ids(pts.size());
      for (int k = int(ids.size()); k--;)
        ids[k] = std::make_pair(real(0), k);
      int cost = 0;
      std::vector<Node> tree;
      init(pts, dist, bucket, tree, ids, cost,
           0, int(ids.size()), int(ids.size()/2));
      _tree.swap(tree);
      _numpoints = int(pts.size());
      _bucket = bucket;
      _mc = _sc = 0;
      _cost = cost; _c1 = _k = _cmax = 0;
      _cmin = std::numeric_limits<int>::max();
    }

    /**
     * Search the NearestNeighbor.
     *
     * @param[in] pts the vector of points used for initialization.
     * @param[in] dist the distance function object used for initialization.
     * @param[in] query the query point.
     * @param[out] ind a vector of indices to the closest points found.
     * @param[in] k the number of points to search for (default = 1).
     * @param[in] maxdist only return points with distances of \e maxdist or
     *   less from \e query (default is the maximum \e real).
     * @param[in] mindist only return points with distances of more than
     *   \e mindist from \e query (default = &minus;1).
     * @param[in] exhaustive whether to do an exhaustive search (default true).
     * @param[in] tol the tolerance on the results (default 0).
     * @return the distance to the closest point found (&minus;1 if no points
     *   are found).
     *
     * The indices returned in \e ind are sorted by distance from \e query
     * (closest first).
     *
     * With \e exhaustive = true and \e tol = 0 (their default values), this
     * finds the indices of \e k closest neighbors to \e query whose distances
     * to \e query are in (\e mindist, \e maxdist].  If \e mindist and \e
     * maxdist have their default values, then these bounds have no effect.  If
     * \e query is one of the points in the tree, then set \e mindist = 0 to
     * prevent this point (and other coincident points) from being returned.
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
     * \warning The arguments \e pts and \e dist must be identical to
     * those used to initialize the NearestNeighbor; if not, the behavior of
     * this function is undefined (however, if the size of \e pts is wrong,
     * this function exits with no results returned).
     **********************************************************************/
    real Search(const std::vector<position>& pts, const distance& dist,
                const position& query,
                std::vector<int>& ind,
                int k = 1,
                real maxdist = std::numeric_limits<real>::max(),
                real mindist = -1,
                bool exhaustive = true,
                real tol = 0) const {
      std::priority_queue<item> results;
      if (_numpoints > 0 && _numpoints == int(pts.size()) &&
          k > 0 && maxdist > mindist) {
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
          if (!( n >= 0 && tau1 >= d )) continue;
          const Node& current = _tree[n];
          real dst = 0;     // to suppress warning about uninitialized variable
          bool exitflag = false, leaf = current.index < 0;
          for (int i = 0; i < (leaf ? _bucket : 1); ++i) {
            int index = leaf ? current.leaves[i] : current.index;
            if (index < 0) break;
            dst = dist(pts[index], query);
            ++c;

            if (dst > mindist && dst <= tau) {
              if (int(results.size()) == k) results.pop();
              results.push(std::make_pair(dst, index));
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
          for (int l = 0; l < 2; ++l) {
            if (current.data.child[l] >= 0 &&
                dst + current.data.upper[l] >= mindist) {
              if (dst < current.data.lower[l]) {
                d = current.data.lower[l] - dst;
                if (tau1 >= d)
                  todo.push(std::make_pair(-d, current.data.child[l]));
              } else if (dst > current.data.upper[l]) {
                d = dst - current.data.upper[l];
                if (tau1 >= d)
                  todo.push(std::make_pair(-d, current.data.child[l]));
              } else
                todo.push(std::make_pair(real(1), current.data.child[l]));
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
     * @return the total number of points in the set.
     **********************************************************************/
    int NumPoints() const { return _numpoints; }

    /**
     * Write the object to an I/O stream.
     *
     * @param[in,out] os the stream to write to.
     * @param[in] bin if true (the default) save in binary mode.
     * @exception std::bad_alloc if memory for the string representation of the
     *   object can't be allocated.
     *
     * The counters tracking the statistics of searches are not saved; however
     * the initializtion cost is saved.  The format of the binary saves is \e
     * not portable.
     *
     * \note <a href="http://www.boost.org/libs/serialization/doc">
     * Boost serialization</a> can also be used to save and restore a
     * NearestNeighbor object.  This requires that the
     * GEOGRAPHICLIB_HAVE_BOOST_SERIALIZATION macro be defined.
     **********************************************************************/
    void Save(std::ostream& os, bool bin = true) const {
      int realspec = std::numeric_limits<real>::digits *
        (std::numeric_limits<real>::is_integer ? -1 : 1);
      if (bin) {
        char id[] = "NearestNeighbor_";
        os.write(id, 16);
        int buf[6];
        buf[0] = version;
        buf[1] = realspec;
        buf[2] = _bucket;
        buf[3] = _numpoints;
        buf[4] = int(_tree.size());
        buf[5] = _cost;
        os.write(reinterpret_cast<const char *>(buf), 6 * sizeof(int));
        for (int i = 0; i < int(_tree.size()); ++i) {
          const Node& node = _tree[i];
          os.write(reinterpret_cast<const char *>(&node.index), sizeof(int));
          if (node.index >= 0) {
            os.write(reinterpret_cast<const char *>(node.data.lower),
                     2 * sizeof(real));
            os.write(reinterpret_cast<const char *>(node.data.upper),
                     2 * sizeof(real));
            os.write(reinterpret_cast<const char *>(node.data.child),
                     2 * sizeof(int));
          } else {
            os.write(reinterpret_cast<const char *>(node.leaves),
                     _bucket * sizeof(int));
          }
        }
      } else {
        std::stringstream ostring;
          // Ensure enough precision for type real.  If real is actually a
          // signed integer type, full precision is used anyway.  With C++11,
          // max_digits10 can be used instead.
        ostring.precision(std::numeric_limits<real>::digits10 + 2);
        ostring << version << " " << realspec << " " << _bucket << " "
                << _numpoints << " " << _tree.size() << " " << _cost;
        for (int i = 0; i < int(_tree.size()); ++i) {
          const Node& node = _tree[i];
          ostring << "\n" << node.index;
          if (node.index >= 0) {
            for (int l = 0; l < 2; ++l)
              ostring << " " << node.data.lower[l] << " " << node.data.upper[l]
                      << " " << node.data.child[l];
          } else {
            for (int l = 0; l < _bucket; ++l)
              ostring << " " << node.leaves[l];
          }
        }
        os << ostring.str();
      }
    }

    /**
     * Read the object from an I/O stream.
     *
     * @param[in,out] is the stream to read from
     * @param[in] bin if true (the default) load in binary mode.
     * @exception GeographicErr if the state read from \e is is illegal.
     *
     * The counters tracking the statistics of searches are reset by this
     * operation.  Binary data must have been saved on a machine with the same
     * architecture.  If an exception is thrown, the state of the
     * NearestNeighbor is unchanged.
     *
     * \note <a href="http://www.boost.org/libs/serialization/doc">
     * Boost serialization</a> can also be used to save and restore a
     * NearestNeighbor object.  This requires that the
     * GEOGRAPHICLIB_HAVE_BOOST_SERIALIZATION macro be defined.
     *
     * \warning The same arguments \e pts and \e dist used for
     * initialization must be provided to the Search() function.
     **********************************************************************/
    void Load(std::istream& is, bool bin = true) {
      int version1, realspec, bucket, numpoints, treesize, cost;
      if (bin) {
        char id[17];
        is.read(id, 16);
        id[16] = '\0';
        if (!(std::strcmp(id, "NearestNeighbor_") == 0))
          throw GeographicLib::GeographicErr("Bad ID");
        is.read(reinterpret_cast<char *>(&version1), sizeof(int));
        is.read(reinterpret_cast<char *>(&realspec), sizeof(int));
        is.read(reinterpret_cast<char *>(&bucket), sizeof(int));
        is.read(reinterpret_cast<char *>(&numpoints), sizeof(int));
        is.read(reinterpret_cast<char *>(&treesize), sizeof(int));
        is.read(reinterpret_cast<char *>(&cost), sizeof(int));
      } else {
        if (!( is >> version1 >> realspec >> bucket >> numpoints >> treesize
               >> cost ))
          throw GeographicLib::GeographicErr("Bad header");
      }
      if (!( version1 == version ))
        throw GeographicLib::GeographicErr("Incompatible version");
      if (!( realspec == std::numeric_limits<real>::digits *
             (std::numeric_limits<real>::is_integer ? -1 : 1) ))
        throw GeographicLib::GeographicErr("Different real types");
      if (!( 0 <= bucket && bucket <= maxbucket ))
        throw GeographicLib::GeographicErr("Bad bucket size");
      if (!( 0 <= treesize && treesize <= numpoints ))
        throw GeographicLib::GeographicErr("Bad number of points or tree size");
      std::vector<Node> tree;
      tree.reserve(treesize);
      for (int i = 0; i < treesize; ++i) {
        Node node;
        if (bin) {
          is.read(reinterpret_cast<char *>(&node.index), sizeof(int));
          if (node.index >= 0) {
            is.read(reinterpret_cast<char *>(node.data.lower),
                    2 * sizeof(real));
            is.read(reinterpret_cast<char *>(node.data.upper),
                    2 * sizeof(real));
            is.read(reinterpret_cast<char *>(node.data.child),
                    2 * sizeof(int));
          } else {
            is.read(reinterpret_cast<char *>(node.leaves),
                    bucket * sizeof(int));
            for (int l = bucket; l < maxbucket; ++l)
              node.leaves[l] = 0;
          }
        } else {
          if (!( is >> node.index ))
            throw GeographicLib::GeographicErr("Bad index");
          if (node.index >= 0) {
            for (int l = 0; l < 2; ++l) {
              if (!( is >> node.data.lower[l] >> node.data.upper[l]
                     >> node.data.child[l] ))
                throw GeographicLib::GeographicErr("Bad node data");
            }
          } else {
            // Must be at least one valid leaf followed by a sequence end
            // markers (-1).
            for (int l = 0; l < bucket; ++l) {
              if (!( is >> node.leaves[l] ))
                throw GeographicLib::GeographicErr("Bad leaf data");
            }
            for (int l = bucket; l < maxbucket; ++l)
              node.leaves[l] = 0;
          }
        }
        node.Check(numpoints, treesize, bucket);
        tree.push_back(node);
      }
      _tree.swap(tree);
      _numpoints = numpoints;
      _bucket = bucket;
      _mc = _sc = 0;
      _cost = cost; _c1 = _k = _cmax = 0;
      _cmin = std::numeric_limits<int>::max();
    }

    /**
     * The accumulated statistics on the searches so far.
     *
     * @param[out] setupcost the cost of initializing the NearestNeighbor.
     * @param[out] numsearches the number of calls to Search().
     * @param[out] searchcost the total cost of the calls to Search().
     * @param[out] mincost the minimum cost of a Search().
     * @param[out] maxcost the maximum cost of a Search().
     * @param[out] mean the mean cost of a Search().
     * @param[out] sd the standard deviation in the cost of a Search().
     *
     * Here "cost" measures the number of distance calculations needed.  Note
     * that the accumulation of statistics is \e not thread safe.
     **********************************************************************/
    void Statistics(int& setupcost, int& numsearches, int& searchcost,
                    int& mincost, int& maxcost,
                    double& mean, double& sd) const {
      setupcost = _cost; numsearches = _k; searchcost = _c1;
      mincost = _cmin; maxcost = _cmax;
      mean = _mc; sd = std::sqrt(_sc / (_k - 1));
    }

    /**
     * Reset the counters for the accumulated statistics on the searches so
     * far.
     **********************************************************************/
    void ResetStatistics() const {
      _mc = _sc = 0;
      _c1 = _k = _cmax = 0;
      _cmin = std::numeric_limits<int>::max();
    }

  private:
    // Package up a real and an int.  We will want to sort on the real so put
    // it first.
    typedef std::pair<real, int> item;
    class Node {
    public:
      struct bounds {
        real lower[2], upper[2]; // bounds on inner/outer distances
        int child[2];
      };
      union {
        bounds data;
        int leaves[maxbucket];
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

      // Sanity check on a Node
      void Check(int numpoints, int treesize, int bucket) const {
        if (!( -1 <= index && index << numpoints ))
          throw GeographicLib::GeographicErr("Bad index");
        if (index >= 0) {
          if (!( -1 <= data.child[0] && data.child[0] < treesize &&
                 -1 <= data.child[1] && data.child[1] < treesize ))
            throw GeographicLib::GeographicErr("Bad child pointers");
          if (!( 0 <= data.lower[0] && data.lower[0] <= data.upper[0] &&
                 data.upper[0] <= data.lower[1] &&
                 data.lower[1] <= data.upper[1] ))
            throw GeographicLib::GeographicErr("Bad bounds");
        } else {
          // Must be at least one valid leaf followed by a sequence end markers
          // (-1).
          bool start = true;
          for (int l = 0; l < bucket; ++l) {
            if (!( (start ?
                    ((l == 0 ? 0 : -1) <= leaves[l] && leaves[l] < numpoints) :
                    leaves[l] == -1) ))
              throw GeographicLib::GeographicErr("Bad leaf data");
            start = leaves[l] >= 0;
          }
          for (int l = bucket; l < maxbucket; ++l) {
            if (leaves[l] != 0)
              throw GeographicLib::GeographicErr("Bad leaf data");
          }
        }
      }

#if GEOGRAPHICLIB_HAVE_BOOST_SERIALIZATION
      friend class boost::serialization::access;
      template<class Archive> void save(Archive& ar, const unsigned int) const {
        ar & boost::serialization::make_nvp("index", index);
        if (index < 0)
          ar & boost::serialization::make_nvp("leaves", leaves);
        else
          ar & boost::serialization::make_nvp("lower", data.lower)
            & boost::serialization::make_nvp("upper", data.upper)
            & boost::serialization::make_nvp("child", data.child);
      }
      template<class Archive> void load(Archive& ar, const unsigned int) {
        ar & boost::serialization::make_nvp("index", index);
        if (index < 0)
          ar & boost::serialization::make_nvp("leaves", leaves);
        else
          ar & boost::serialization::make_nvp("lower", data.lower)
            & boost::serialization::make_nvp("upper", data.upper)
            & boost::serialization::make_nvp("child", data.child);
      }
      template<class Archive>
      void serialize(Archive& ar, const unsigned int file_version)
      { boost::serialization::split_member(ar, *this, file_version); }
#endif
    };
#if GEOGRAPHICLIB_HAVE_BOOST_SERIALIZATION
    friend class boost::serialization::access;
    template<class Archive> void save(Archive& ar, const unsigned) const {
      int realspec = std::numeric_limits<real>::digits *
        (std::numeric_limits<real>::is_integer ? -1 : 1);
      // Need to use version1, otherwise load error in debug mode on Linux:
      // undefined reference to GeographicLib::NearestNeighbor<...>::version.
      int version1 = version;
      ar & boost::serialization::make_nvp("version", version1)
        & boost::serialization::make_nvp("realspec", realspec)
        & boost::serialization::make_nvp("bucket", _bucket)
        & boost::serialization::make_nvp("numpoints", _numpoints)
        & boost::serialization::make_nvp("cost", _cost)
        & boost::serialization::make_nvp("tree", _tree);
    }
    template<class Archive> void load(Archive& ar, const unsigned) {
      int version1, realspec, bucket, numpoints, cost;
      ar & boost::serialization::make_nvp("version", version1);
      if (version1 != version)
        throw GeographicLib::GeographicErr("Incompatible version");
      std::vector<Node> tree;
      ar & boost::serialization::make_nvp("realspec", realspec);
      if (!( realspec == std::numeric_limits<real>::digits *
             (std::numeric_limits<real>::is_integer ? -1 : 1) ))
        throw GeographicLib::GeographicErr("Different real types");
      ar & boost::serialization::make_nvp("bucket", bucket);
      if (!( 0 <= bucket && bucket <= maxbucket ))
        throw GeographicLib::GeographicErr("Bad bucket size");
      ar & boost::serialization::make_nvp("numpoints", numpoints)
        & boost::serialization::make_nvp("cost", cost)
        & boost::serialization::make_nvp("tree", tree);
      if (!( 0 <= int(tree.size()) && int(tree.size()) <= numpoints ))
        throw GeographicLib::GeographicErr("Bad number of points or tree size");
      for (int i = 0; i < int(tree.size()); ++i)
        tree[i].Check(numpoints, int(tree.size()), bucket);
      _tree.swap(tree);
      _numpoints = numpoints;
      _bucket = bucket;
      _mc = _sc = 0;
      _cost = cost; _c1 = _k = _cmax = 0;
      _cmin = std::numeric_limits<int>::max();
    }
    template<class Archive>
    void serialize(Archive& ar, const unsigned int file_version)
    { boost::serialization::split_member(ar, *this, file_version); }
#endif

    int _numpoints, _bucket, _cost;
    std::vector<Node> _tree;
    // Counters to track stastistics on the cost of searches
    mutable double _mc, _sc;
    mutable int _c1, _k, _cmin, _cmax;

    int init(const std::vector<position>& pts, const distance& dist, int bucket,
             std::vector<Node>& tree, std::vector<item>& ids, int& cost,
             int l, int u, int vp) {

      if (u == l)
        return -1;
      Node node;

      if (u - l > (bucket == 0 ? 1 : bucket)) {

        // choose a vantage point and move it to the start
        int i = vp;
        std::swap(ids[l], ids[i]);

        int m = (u + l + 1) / 2;

        for (int k = l + 1; k < u; ++k) {
          ids[k].first = dist(pts[ids[l].second], pts[ids[k].second]);
          ++cost;
        }
        // partitian around the median distance
        std::nth_element(ids.begin() + l + 1, ids.begin() + m, ids.begin() + u);
        node.index = ids[l].second;
        if (m > l + 1) {        // node.child[0] is possibly empty
          typename std::vector<item>::iterator
            t = std::min_element(ids.begin() + l + 1, ids.begin() + m);
          node.data.lower[0] = t->first;
          t = std::max_element(ids.begin() + l + 1, ids.begin() + m);
          node.data.upper[0] = t->first;
          // Use point with max distance as vantage point; this point act as a
          // "corner" point and leads to a good partition.
          node.data.child[0] = init(pts, dist, bucket, tree, ids, cost,
                                    l + 1, m, int(t - ids.begin()));
        }
        typename std::vector<item>::iterator
          t = std::max_element(ids.begin() + m, ids.begin() + u);
        node.data.lower[1] = ids[m].first;
        node.data.upper[1] = t->first;
        // Use point with max distance as vantage point here too
        node.data.child[1] = init(pts, dist, bucket, tree, ids, cost,
                                  m, u, int(t - ids.begin()));
      } else {
        if (bucket == 0)
          node.index = ids[l].second;
        else {
          node.index = -1;
          // Sort the bucket entries so that the tree is independent of the
          // implementation of nth_element.
          std::sort(ids.begin() + l, ids.begin() + u);
          for (int i = l; i < u; ++i)
            node.leaves[i-l] = ids[i].second;
          for (int i = u - l; i < bucket; ++i)
            node.leaves[i] = -1;
          for (int i = bucket; i < maxbucket; ++i)
            node.leaves[i] = 0;
        }
      }

      tree.push_back(node);
      return int(tree.size()) - 1;
    }

  };

} // namespace GeographicLib

#if defined(_MSC_VER)
#  pragma warning (pop)
#endif

#endif  // GEOGRAPHICLIB_NEARESTNEIGHBOR_HPP
