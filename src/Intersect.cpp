#include <GeographicLib/Intersect.hpp>
#include <limits>
#include <utility>
#include <algorithm>
#include <set>
#include <cassert>
#include <iomanip>

using namespace std;

namespace GeographicLib {

  Intersect::Intersect(const Geodesic& geod)
    : _geod(geod)
    , _a(_geod.EquatorialRadius())
    , _f(_geod.Flattening())
    , _n(_f / (2 - _f))
    , _R(sqrt(_geod.EllipsoidArea() / (4 * Math::pi())))
    , _d(_R * Math::pi())       // Used to normalize intersection points
    , _tol(_d * pow(numeric_limits<real>::epsilon(), 3/real(4)))
    , _slop(_d / 1000)
    , _debug(true)
    , _comp(_slop)
    , cnt0(0)
    , cnt1(0)
    , cnt2(0)
  {
    _s1 = _s4 = _a * (1 - _f) * Math::pi();
    _s3 = 2 * distpolar(90);
    _geod.Inverse(0, 0, 90, 0, _s5); _s5 *= 2;
    if (_f > 0) {
      _s2 = distoblique();
      _s4 = _s1;
    } else {
      _s2 = _s5;
      _s4 = polarb();
      swap(_s1, _s3);
    }
    _c1 = _s3/2 + _slop;
    _c2 = 2*_s2/3 + _slop;
    _c3 = _s4 - _slop;
    if (! (2 * _c1 > _s3 && _c1 < _s4 &&
           3 * _c2 > _s2 && _c2 < _s4 && _c2 < 2 * _s1) )
      throw GeographicErr("Ellipsoid too eccentric for Closest");
  }

  Intersect::Point
  Intersect::Closest(Math::real latX, Math::real lonX, Math::real aziX,
                     Math::real latY, Math::real lonY, Math::real aziY,
                     const Intersect::Point& p0) const {
    return Closest(_geod.Line(latX, lonX, aziX, LineCaps),
                   _geod.Line(latY, lonY, aziY, LineCaps),
                   p0);
  }

  Intersect::Point
  Intersect::Closest(const GeodesicLine& lineX, const GeodesicLine& lineY,
                     const Intersect::Point& p0) const {
    int flag;
    return Solve2(lineX, lineY, XPoint(p0), flag).data();
  }

  Intersect::Point
  Intersect::Segment(Math::real latX1, Math::real lonX1,
                     Math::real latX2, Math::real lonX2,
                     Math::real latY1, Math::real lonY1,
                     Math::real latY2, Math::real lonY2,
                     int& segmode) const {
    return Segment(_geod.InverseLine(latX1, lonX1, latX2, lonX2, LineCaps),
                   _geod.InverseLine(latY1, lonY1, latY2, lonY2, LineCaps),
                   segmode);
  }

  Intersect::Point
  Intersect::Segment(const GeodesicLine& lineX,
                     const GeodesicLine& lineY, int& segmode) const {
    int flag;
    return Solve4(lineX, lineY, segmode, flag).data();
  }

  Intersect::Point
  Intersect::Next(Math::real latX, Math::real lonX,
                  Math::real aziX, Math::real aziY) const {
    return Next(_geod.Line(latX, lonX, aziX, LineCaps),
                _geod.Line(latX, lonX, aziY, LineCaps));
  }

  Intersect::Point
  Intersect::Next(const GeodesicLine& lineX, const GeodesicLine& lineY)
    const {
    int flag;
    return Solve3(lineX, lineY, flag).data();
  }

  std::vector<Intersect::Point>
  Intersect::All(Math::real latX, Math::real lonX, Math::real aziX,
                 Math::real latY, Math::real lonY, Math::real aziY,
                 Math::real maxdist, const Point& p0) const {
    return All(_geod.Line(latX, lonX, aziX, LineCaps),
               _geod.Line(latY, lonY, aziY, LineCaps),
               maxdist, p0);
  }

  std::vector<Intersect::Point>
  Intersect::All(const GeodesicLine& lineX, const GeodesicLine& lineY,
                 Math::real maxdist, const Point& p0) const {
    int flag;
    auto s = Solve5(lineX, lineY, fmax(real(0), maxdist), XPoint(p0), flag);
    std::vector<Intersect::Point> v(s.size());
    int i = 0;
    for (auto p = s.cbegin(); p != s.cend(); ++p)
      v[i++] = (*p).data();
    sort(v.begin(), v.end(), RankPoint(p0));
    return v;
  }

  Intersect::XPoint
  Intersect::Solve4(const GeodesicLine& lineX, const GeodesicLine& lineY,
                    int& segmode, int& flag) const {
    real sx = lineX.Distance()/2, sy = lineY.Distance()/2;
    XPoint p0 = XPoint(sx/2, sy/2), q = Solve2(lineX, lineY, p0, flag);
    q = fixsegment(sx, sy, q, flag);
    segmode = segmentmode(sx, sy, q);
    if (segmode != 0 && q.Dist(p0) <= (sx + sy)/2) {
      int flag1 = 0, segmode1 = 1;
      XPoint q1;
      for (int ix = 0; ix < 2 && segmode1 != 0; ++ix) {
        for (int iy = 0; iy < 2 && segmode1 != 0; ++iy) {
          XPoint t(ix * sx, iy * sy);
          if (q.Dist(t) >= 2 *_s1) {
            q1 = Solve1(lineX, lineY, t, flag1);
            q1 = fixcoincident(t, q1, flag1);
            segmode1 = segmentmode(sx, sy, q1);
          }
        }
      }
      if (segmode1 == 0) { segmode = 0; q = q1; }
    }
    return q;
  }

  Intersect::XPoint
  Intersect::Solve2(const GeodesicLine& lineX, const GeodesicLine& lineY,
                    const Intersect::XPoint& p0, int& flag) const {
    const int num = 5;
    const int ix[num] = { 0, -1,  0,  1,  0 };
    const int iy[num] = { 0,  0,  1,  0, -1 };
    bool    skip[num] = { 0,  0,  0,  0,  0 };
    XPoint q;                    // Best intersection so far
    for (int n = 0; n < num; ++n) {
      if (skip[n]) continue;
      XPoint q1 = Solve1(lineX, lineY, p0 + XPoint(ix[n] * _c1, iy[n] * _c1),
                        flag);
      q1 = fixcoincident(p0, q1, flag);
      if (_comp.eq(q, q1)) continue;
      if (q1.Dist(p0) < _s1) { q = q1; ++cnt2; break; }
      if (n == 0 || q1.Dist(p0) < q.Dist(p0)) { q = q1; ++cnt2; }
      for (int m = n + 1; m < num; ++m)
        skip[m] = skip[m] ||
          q1.Dist(p0 + XPoint(ix[m] * _c1, iy[m] * _c1))
          < 2*_s1 - _c1 - _slop;
    }
    return q;
  }

  Intersect::XPoint
  Intersect::Solve3(const GeodesicLine& lineX, const GeodesicLine& lineY,
                    int& flag)
    const {
    const int num = 8;
    const int ix[num] = { -1, -1,  1,  1, -2,  0,  2,  0 };
    const int iy[num] = { -1,  1, -1,  1,  0,  2,  0, -2 };
    bool    skip[num] = {  0,  0,  0,  0,  0,  0,  0,  0 };
    XPoint z(0,0),               // for excluding the origin
      q;                        // Best intersection so far
    for (int n = 0; n < num; ++n) {
      if (skip[n]) continue;
      XPoint q1 = Solve1(lineX, lineY, XPoint(ix[n] * _c2, iy[n] * _c2), flag);
      if (_comp.eq(z, q1)) continue;
      if (flag && n < num/2) {
        // For coincident geodesics, 2 out of first four trials will return
        // flag != 0; replace intersection with conjugate point.
        real s = ConjugateDist(lineX, q1.x, false);
        q1 = XPoint(s, flag*s);
      }
      //      if (_comp.eq(q, q1)) continue;
      if (n == 0 || q1.Dist() < q.Dist()) { q = q1; ++cnt2;}
      for (int m = n + 1; m < num; ++m)
        skip[m] = skip[m] ||
          q1.Dist(XPoint(ix[m] * _c2, iy[m] * _c2))
          < 2*_s1 - _c2 - _slop;
    }
    return q;
  }

  Intersect::XPoint
  Intersect::Solve1(const GeodesicLine& lineX, const GeodesicLine& lineY,
                    const Intersect::XPoint& p0, int& flag) const {
    ++cnt1;
    XPoint q = p0;
    for (int n = 0; n < 100; ++n) {
      ++cnt0;
      XPoint dq = Solve0(lineX, lineY, q, flag);
      q += dq;
      if (flag || !(dq.Dist() > _tol)) break; // break if nan
    }
    return q;
  }

  Intersect::XPoint
  Intersect::Solve0(const GeodesicLine& lineX, const GeodesicLine& lineY,
                    const Intersect::XPoint& p, int& flag) const {
    // threshold for coincident geodesics and intersections; this corresponds
    // to about 4.3 nm on WGS84.
    static const real eps = 3*numeric_limits<real>::epsilon();
    real latX, lonX, aziX, latY, lonY, aziY;
    lineX.Position(p.x , latX, lonX, aziX);
    lineY.Position(p.y, latY, lonY, aziY);
    real z, aziXa, aziYa;
    _geod.Inverse(latX, lonX, latY, lonY, z, aziXa, aziYa);
    real sinz = sin(z/_R), cosz = cos(z/_R);
    // X = interior angle at X, Y = exterior angle at Y
    real dX, dY, dXY,
      X = Math::AngDiff(aziX, aziXa, dX), Y = Math::AngDiff(aziY, aziYa, dY),
      XY = Math::AngDiff(X, Y, dXY);
    real s = copysign(real(1), XY + (dXY + dY - dX)); // inverted triangle
    // For z small, sinz -> z, cosz -> 1
    // ( sinY*cosX*cosz - cosY*sinX) =
    // (-sinX*cosY*cosz + cosX*sinY) -> sin(Y-X)
    // for z = pi, sinz -> 0, cosz -> -1
    // ( sinY*cosX*cosz - cosY*sinX) -> -sin(Y+X)
    // (-sinX*cosY*cosz + cosX*sinY) ->  sin(Y+X)
    real sinX, cosX; Math::sincosde(s*X, s*dX, sinX, cosX);
    real sinY, cosY; Math::sincosde(s*Y, s*dY, sinY, cosY);
    real sX, sY;
    if (z <= eps * _R) {
      sX = sY = 0;              // Already at intersection
      // Determine whether lineX and lineY are parallel or antiparallel
      if (fabs(sinX - sinY) <= eps && fabs(cosX - cosY) <= eps)
        flag = 1;
      else if (fabs(sinX + sinY) <= eps && fabs(cosX + cosY) <= eps)
        flag = -1;
      else
        flag = 0;
    } else if (fabs((sinX) <= eps && fabs(sinY) <= eps)) {
      flag = cosX * cosY > 0 ? 1 : -1;
      // Coincident geodesics, place intersection at midpoint
      sX =  cosX * z/2; sY = -cosY * z/2;
      // alt1: sX =  cosX * z; sY = 0;
      // alt2: sY = -cosY * z; sX = 0;
    } else {
      // General case.  [SKIP: Divide args by |sinz| to avoid possible
      // underflow in {sinX,sinY}*sinz; this is probably not necessary].
      // Definitely need to treat sinz < 0 (z > pi*R) correctly.  Without
      // this we have some convergence failures in Solve.
      sX = _R * atan2( sinY * sinz,  sinY*cosX*cosz - cosY*sinX);
      sY = _R * atan2( sinX * sinz, -sinX*cosY*cosz + cosX*sinY);
      flag = 0;
    }
    return XPoint(sX, sY);
  }

  set<Intersect::XPoint, Intersect::SetComp>
  Intersect::Solve5(const GeodesicLine& lineX,
                    const GeodesicLine& lineY,
                    Math::real maxdist, const XPoint& p0,
                    int& flag)  const {
    bool debug = false;
    real maxdistx = fmax(maxdist, _slop);
    const int m = int(ceil(maxdistx/_c3)), // process m x m set of tiles
      m2 = m*m + (m - 1) % 2,              // add center tile if m is even
      n = m - 1;                           // Range of i, j = [-n:2:n]
    real c3 = maxdistx/m;                  // c3 <= _c3
    if (debug) cerr << setprecision(16) << "BEGIN "
                    << _d << " " << _s4 << " " << _c3 << " " << c3 << "\n";
    vector<XPoint> start(m2);
    vector<bool> skip(m2, false);
    int k = 0, flag0 = 0;
    start[k++] = p0;
    for (int i = -n; i <= n; i += 2)
      for (int j = -n; j <= n; j += 2) {
        if (!(i == 0 && j == 0))
          start[k++] = p0 + XPoint( c3 * (i + j) / 2, c3 * (i - j) / 2);
      }
    if (debug) cerr << "0 " << m << " "  << n << " " <<  k << " " << m2 << "\n";
    assert(k == m2);
    set<XPoint, SetComp> r(_comp); // Intersections found
    set<XPoint, SetComp> c(_comp); // Closest coincident intersections
    vector<XPoint> added;
    for (int k = 0; k < m2; ++k) {
      if (skip[k]) continue;
      XPoint q = Solve1(lineX, lineY, start[k], flag);
      if (debug) cerr << "MAIN " << k << " "
                      << start[k].x << " " << start[k].y << " "
                      << q.x << " " << q.y << " " << flag << "\n";
      if (r.find(q) != r.end()  // intersection already found
          // or it's on a line of coincident intersections already processed
          || (flag0 != 0 && c.find(fixcoincident(p0, q, flag0)) != c.end()))
        continue;
      added.clear();
      if (flag != 0) {
        // This value of flag must be constitent with flag0
        assert(flag0 == 0 || flag0 == flag);
        flag0 = flag;
        if (debug) cerr << "A " << q.x << " " << q.y << "\n";
        // Process coincident intersections
        q = fixcoincident(p0, q, flag0);
        c.insert(q);
        if (debug) cerr << "B " << q.x << " " << q.y << "\n";
        // Elimate all existing intersections on this line (which
        // didn't set flag0).
        if (1) {
        for (auto qp = r.begin(); qp != r.end(); ) {
          if (_comp.eq(fixcoincident(p0, *qp, flag0), q)) {
            if (debug) cerr << "C " << qp->x << " " << qp->y << "\n";
            qp = r.erase(qp);
          }
          else
            ++qp;
        }
        }
        real s0 = q.x;
        XPoint qc;
        real t, m12, M12, M21;
        lineX.GenPosition(false, s0,
                          GeodesicLine::REDUCEDLENGTH |
                          GeodesicLine::GEODESICSCALE,
                          t, t, t, t, m12, M12, M21, t);
        // Compute line of conjugate points
        for (int sgn = -1; sgn <= 1; sgn += 2) {
          real sa = 0;
          do {
            sa = ConjugateDist(lineX, s0 + sa + sgn*_d, false, m12, M12, M21)
              - s0;
            if (debug) cerr << "X " << sgn << " " << sa-s0 << "\n";
            qc = q + XPoint(sa, flag0*sa);
            added.push_back(qc);
            r.insert(qc);
            if (debug) cerr << "D " << qc.x << " " << qc.y << "\n";
          } while (qc.Dist(p0) <= maxdistx);
        }
      }
      added.push_back(q);
      r.insert(q);
      if (1) {
        for (auto qp = added.cbegin();  qp != added.cend(); ++qp) {
          for (int l = k + 1; l < m2; ++l)
            skip[l] = skip[l] || qp->Dist(start[l]) < 2*_s1 - c3 - _slop;
        }
      }
    }
    // Trim intersections to maxdist
    for (auto qp = r.begin(); qp != r.end(); ) {
      if (!(qp->Dist(p0) <= maxdist))
        qp = r.erase(qp);
      else
        ++qp;
    }
    return r;
  }

  Math::real Intersect::distpolar(Math::real lat1, Math::real* lat2)
    const {
    GeodesicLine line = _geod.Line(lat1, 0, 0,
                                   GeodesicLine::REDUCEDLENGTH |
                                   GeodesicLine::GEODESICSCALE |
                                   GeodesicLine::DISTANCE_IN);
    real s = ConjugateDist(line, (1 + _f/2) * _a * Math::pi() / 2, true);
    if (lat2) {
      real t;
      line.GenPosition(false, s, GeodesicLine::LATITUDE,
                       *lat2, t, t, t, t, t, t, t);
    }
    return s;
  }

  Math::real Intersect::polarb(Math::real* lata, Math::real* latb) const {
    if (_f == 0) {
      if (lata) *lata = 64;
      if (latb) *latb = 90-64;
      return _d;
    }
    real
      lat0 = 63, s0 = distpolar(lat0),
      lat1 = 65, s1 = distpolar(lat1),
      lat2 = 64, s2 = distpolar(lat2),
      latx = lat2, sx = s2;
    // Solve for ds(lat)/dlat = 0 with a quadratic fit
    for (int i = 0; i < 10; ++i) {
      real den = (lat1-lat0)*s2 + (lat0-lat2)*s1 + (lat2-lat1)*s0;
      if (!(den < 0 || den > 0)) break; // Break if nan
      real latn = ((lat1-lat0)*(lat1+lat0)*s2 + (lat0-lat2)*(lat0+lat2)*s1 +
                   (lat2-lat1)*(lat2+lat1)*s0) / (2*den);
      lat0 = lat1; s0 = s1;
      lat1 = lat2; s1 = s2;
      lat2 = latn; s2 = distpolar(lat2);
      if (_f < 0 ? (s2 < sx) : (s2 > sx)) {
        sx = s2;
        latx = lat2;
      }
    }
    if (lata) *lata = latx;
    if (latb) distpolar(latx, latb);
    return 2 * sx;
  }

  // Find {semi-,}conjugate point relative to s0 which is close to s1.
  Math::real Intersect::ConjugateDist(const GeodesicLine& line, Math::real s3,
                                      bool semi, Math::real m12,
                                      Math::real M12, Math::real M21) const {
    // semi = false: solve for m23 = 0 using dm23/ds3 = M32
    // semi = true : solve for M23 = 0 using dM23/ds3 = - (1 - M23*M32)/m23
    // Here 2 is point with given m12, M12, M21 and default values s.t. point 2
    // = point 1.
    real s = s3;
    for (int i = 0; i < 100; ++i) {
      real t, m13, M13, M31;
      line.GenPosition(false, s,
                       GeodesicLine::REDUCEDLENGTH |
                       GeodesicLine::GEODESICSCALE,
                       t, t, t, t, m13, M13, M31, t);
      real
        // See "Algorithms for geodesics", eqs. 31, 32, 33.
        m23 = m13 * M12 - m12 * M13,
        // when m12 -> eps, (1 - M12 * M21) -> eps^2, I suppose.
        M23 = M13 * M21 + (m12 == 0 ? 0 : (1 - M12 * M21) * m13/m12),
        M32 = M31 * M12 + (m13 == 0 ? 0 : (1 - M13 * M31) * m12/m13);
      real ds = semi ? m23 * M23 / (1 - M23*M32) : -m23 / M32;
      s = s + ds;
      if (!(fabs(ds) > _tol)) break;
    }
    return s;
  }

  Math::real Intersect::conjdist(Math::real azi,
                                 Math::real* ds,
                                 Math::real* sp, Math::real* sm) const {
    GeodesicLine line = _geod.Line(0, 0, azi, LineCaps);
    real s = ConjugateDist(line, _d, false);
    if (ds) {
      int flag;
      XPoint p = Solve1(line, line, XPoint(s/2, -3*s/2), flag);
      if (sp) *sp = p.x;
      if (sm) *sm = p.y;
      *ds = p.Dist() - 2*s;
    }
    return s;
  }

  Math::real Intersect::distoblique(Math::real* azi,
                                    Math::real* sp,
                                    Math::real* sm) const {
    if (_f == 0) {
      if (azi) *azi = 45;
      if (sp) *sp = 0.5;
      if (sm) *sm = -1.5;
      return _d;
    }
    real sa, sb,
      azi0 = 46, ds0, s0 = conjdist(azi0, &ds0, &sa, &sb),
      azi1 = 44, ds1, s1 = conjdist(azi1, &ds1, &sa, &sb),
      azix = azi1, dsx = fabs(ds1), sx = s1, sax = sa, sbx = sb;
    // find ds(azi) = 0 by secant method
    (void) s0;
    for (int i = 0; i < 10 && ds1 != ds0; ++i) {
      real azin = (azi0*ds1-azi1*ds0)/(ds1-ds0);
      azi0 = azi1; s0 = s1; ds0 = ds1;
      azi1 = azin; s1 = conjdist(azi1, &ds1, &sa, &sb);
      if (fabs(ds1) < dsx) {
        azix = azi1, sx = s1, dsx = fabs(ds1);
        sax = sa; sbx = sb;
        if (ds1 == 0) break;
      }
    }
    if (azi) *azi = azix;
    if (sp) *sp = sax;
    if (sm) *sm = sbx;
    return sx;
  }

  Intersect::XPoint
  Intersect::fixcoincident(const Intersect::XPoint& p0,
                           const Intersect::XPoint& p, int flag) {
    if (flag == 0) return p;
    // eqs : [p0x-p1x = -f*(p0y-p1y), p1x = px+s, p1y = py+f*s]$
    // sol : solve(eqs,[s,p1x,p1y]);
    // =>
    // sol:[ s = ((p0x+f*p0y) - (px+f*py))/2,
    //       p1x = px +     ((p0x+f*p0y) - (px+f*py))/2,
    //       p1y = py + f * ((p0x+f*p0y) - (px+f*py))/2
    // ];
    real s = ((p0.x + flag * p0.y) - (p.x+ flag * p.y))/2;
    return p + XPoint(s, flag*s);
  }

  Intersect::XPoint
  Intersect::fixsegment(Math::real sx, Math::real sy,
                        const Intersect::XPoint& p, int flag) {
    if (flag == 0) return p;
    // eq0: [p1x = px+s, p1y = py+f*s]$
    // solx0:linsolve(cons(p1x=0 ,eq0),[s,p1x,p1y]);
    // solx1:linsolve(cons(p1x=sx,eq0),[s,p1x,p1y]);
    // soly0:linsolve(cons(p1y=0 ,eq0),[s,p1x,p1y]);
    // soly1:linsolve(cons(p1y=sy,eq0),[s,p1x,p1y]);
    // solx0:[s = -px      ,p1x = 0 ,p1y = py-f*px     ];
    // solx1:[s = sx-px    ,p1x = sx,p1y = py-f*(px-sx)];
    // soly0:[s = -f*py    ,p1x = px-f*py     ,p1y = 0 ];
    // soly1:[s = f*(sy-py),p1x = px-f*(py-sy),p1y = sy];
    real
      pya = p.y - flag *  p.x,     sa =             -p.x,  // pxa = 0
      pyb = p.y - flag * (p.x-sx), sb =         sx - p.x,  // pxb = sx
      pxc = p.x - flag *  p.y,     sc = flag *      -p.y,  // pyc = 0
      pxd = p.x - flag * (p.y-sy), sd = flag * (sy - p.y); // pyd = sy
    bool
      ga = 0 <= pya && pya <= sy,
      gb = 0 <= pyb && pyb <= sy,
      gc = 0 <= pxc && pxc <= sx,
      gd = 0 <= pxd && pxd <= sx;
    real s;
    // Test opposite sides of the rectangle first
    if      (ga && gb) s = (sa + sb) / 2;
    else if (gc && gd) s = (sc + sd) / 2;
    else if (ga && gc) s = (sa + sc) / 2;
    else if (ga && gd) s = (sa + sd) / 2;
    else if (gb && gc) s = (sb + sc) / 2;
    else if (gb && gd) s = (sb + sd) / 2;
    else {
      // Intersection not within segments; place intersection in smallest gap.
      if (flag > 0) {
        // distance from p to corner p0 is abs( (px - py) - (p0x - p0y) )
        // consider corners p0 = [0, sy] and p0 = [sx, 0]
        if (fabs((p.x - p.y) + sy) < fabs((p.x - p.y) - sx))
          s = (sy - (p.x + p.y))/2;
        else
          s = (sx - (p.x + p.y))/2;
      } else {
        // distance from p to corner p0 is abs( (px + p.y) - (p0x + p0y) )
        // consider corners p0 = [0, 0] and p0 = [sx, sy]
        if (fabs(p.x + p.y) < fabs((p.x + p.y) - (sx + sy)))
          s = (0 - (p.x - p.y))/2;
        else
          s = ((sx - sy) - (p.x - p.y))/2;
      }
    }
    return p + XPoint(s, flag*s);
  }

}
