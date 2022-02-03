#include <iostream>
#include <limits>
#include <string>
#include <GeographicLib/Math.hpp>
#include <GeographicLib/Utility.hpp>
#include <GeographicLib/DMS.hpp>

using namespace std;
using namespace GeographicLib;

typedef Math::real T;

bool equiv(T x, T y) {
  return (isnan(x) && isnan(y)) ||
    (x == y && signbit(x) == signbit(y));
}

// use "do { } while (false)" idiom so it can be punctuated like a statement.

#define check(expr, r) do {                        \
    T s = T(r),  t = expr;                         \
    if (!equiv(s, t)) {                            \
      cout << "Line " << __LINE__ << ": " << #expr \
           << " != " << #r << " (" << t << ")\n";   \
      ++n;                                         \
    }                                              \
  } while (false)

// use "do { } while (false)" idiom so it can be punctuated like a statement.

#define strcheck(expr, r) do {                        \
    string s = string(r),  t = expr;                         \
    if (!(s == t)) {                                           \
      cout << "Line " << __LINE__ << ": " << #expr \
           << " != " << #r << " (" << t << ")\n";   \
      ++n;                                         \
    }                                              \
  } while (false)

int main() {
  T inf = Math::infinity(),
    nan = Math::NaN(),
    eps = numeric_limits<T>::epsilon(),
    ovf = 1 / Math::sq(eps);
  int n = 0;

  check( Math::AngRound(-eps/32), -eps/32);
  check( Math::AngRound(-eps/64), -0.0   );
  check( Math::AngRound(-  T(0)), -0.0   );
  check( Math::AngRound(   T(0)), +0.0   );
  check( Math::AngRound( eps/64), +0.0   );
  check( Math::AngRound( eps/32), +eps/32);
  check( Math::AngRound((1-2*eps)/64), (1-2*eps)/64);
  check( Math::AngRound((1-eps  )/64),  T(1)    /64);
  check( Math::AngRound((1-eps/2)/64),  T(1)    /64);
  check( Math::AngRound((1-eps/4)/64),  T(1)    /64);
  check( Math::AngRound( T(1)    /64),  T(1)    /64);
  check( Math::AngRound((1+eps/2)/64),  T(1)    /64);
  check( Math::AngRound((1+eps  )/64),  T(1)    /64);
  check( Math::AngRound((1+2*eps)/64), (1+2*eps)/64);
  check( Math::AngRound((1-eps  )/32), (1-eps  )/32);
  check( Math::AngRound((1-eps/2)/32),  T(1)    /32);
  check( Math::AngRound((1-eps/4)/32),  T(1)    /32);
  check( Math::AngRound( T(1)    /32),  T(1)    /32);
  check( Math::AngRound((1+eps/2)/32),  T(1)    /32);
  check( Math::AngRound((1+eps  )/32), (1+eps  )/32);
  check( Math::AngRound((1-eps  )/16), (1-eps  )/16);
  check( Math::AngRound((1-eps/2)/16), (1-eps/2)/16);
  check( Math::AngRound((1-eps/4)/16),  T(1)    /16);
  check( Math::AngRound( T(1)    /16),  T(1)    /16);
  check( Math::AngRound((1+eps/4)/16),  T(1)    /16);
  check( Math::AngRound((1+eps/2)/16),  T(1)    /16);
  check( Math::AngRound((1+eps  )/16), (1+eps  )/16);
  check( Math::AngRound((1-eps  )/ 8), (1-eps  )/ 8);
  check( Math::AngRound((1-eps/2)/ 8), (1-eps/2)/ 8);
  check( Math::AngRound((1-eps/4)/ 8),  T(1)    / 8);
  check( Math::AngRound((1+eps/2)/ 8),  T(1)    / 8);
  check( Math::AngRound((1+eps  )/ 8), (1+eps  )/ 8);
  check( Math::AngRound( 1-eps      ),  1-eps      );
  check( Math::AngRound( 1-eps/2    ),  1-eps/2    );
  check( Math::AngRound( 1-eps/4    ),  1          );
  check( Math::AngRound( T(1)       ),  1          );
  check( Math::AngRound( 1+eps/4    ),  1          );
  check( Math::AngRound( 1+eps/2    ),  1          );
  check( Math::AngRound( 1+eps      ),  1+  eps    );
  check( Math::AngRound(T(90)-64*eps),  90-64*eps  );
  check( Math::AngRound(T(90)-32*eps),  90         );
  check( Math::AngRound(T(90)       ),  90         );

  check( Math::sind(-  inf ),  nan);
  check( Math::sind(-T(720)), -0.0);
  check( Math::sind(-T(540)), -0.0);
  check( Math::sind(-T(360)), -0.0);
  check( Math::sind(-T(180)), -0.0);
  check( Math::sind(-T(  0)), -0.0);
  check( Math::sind(+T(  0)), +0.0);
  check( Math::sind(+T(180)), +0.0);
  check( Math::sind(+T(360)), +0.0);
  check( Math::sind(+T(540)), +0.0);
  check( Math::sind(+T(720)), +0.0);
  check( Math::sind(+  inf ),  nan);

  check( Math::cosd(-  inf ),  nan);
  check( Math::cosd(-T(810)), +0.0);
  check( Math::cosd(-T(630)), +0.0);
  check( Math::cosd(-T(450)), +0.0);
  check( Math::cosd(-T(270)), +0.0);
  check( Math::cosd(-T( 90)), +0.0);
  check( Math::cosd(+T( 90)), +0.0);
  check( Math::cosd(+T(270)), +0.0);
  check( Math::cosd(+T(450)), +0.0);
  check( Math::cosd(+T(630)), +0.0);
  check( Math::cosd(+T(810)), +0.0);
  check( Math::cosd(+  inf ),  nan);

  check( Math::tand(-  inf ),  nan);
  check( Math::tand(-T(810)), -ovf);
  check( Math::tand(-T(720)), -0.0);
  check( Math::tand(-T(630)), +ovf);
  check( Math::tand(-T(540)), +0.0);
  check( Math::tand(-T(450)), -ovf);
  check( Math::tand(-T(360)), -0.0);
  check( Math::tand(-T(270)), +ovf);
  check( Math::tand(-T(180)), +0.0);
  check( Math::tand(-T( 90)), -ovf);
  check( Math::tand(-T(  0)), -0.0);
  check( Math::tand(+T(  0)), +0.0);
  check( Math::tand(+T( 90)), +ovf);
  check( Math::tand(+T(180)), -0.0);
  check( Math::tand(+T(270)), -ovf);
  check( Math::tand(+T(360)), +0.0);
  check( Math::tand(+T(450)), +ovf);
  check( Math::tand(+T(540)), -0.0);
  check( Math::tand(+T(630)), -ovf);
  check( Math::tand(+T(720)), +0.0);
  check( Math::tand(+T(810)), +ovf);
  check( Math::tand(+  inf ),  nan);

  check( Math::atan2d(+T(0), -T(0)), +180 );
  check( Math::atan2d(-T(0), -T(0)), -180 );
  check( Math::atan2d(+T(0), +T(0)), +0.0 );
  check( Math::atan2d(-T(0), +T(0)), -0.0 );
  check( Math::atan2d(+T(0), -T(1)), +180 );
  check( Math::atan2d(-T(0), -T(1)), -180 );
  check( Math::atan2d(+T(0), +T(1)), +0.0 );
  check( Math::atan2d(-T(0), +T(1)), -0.0 );
  check( Math::atan2d(-T(1), +T(0)),  -90 );
  check( Math::atan2d(-T(1), -T(0)),  -90 );
  check( Math::atan2d(+T(1), +T(0)),  +90 );
  check( Math::atan2d(+T(1), -T(0)),  +90 );
  check( Math::atan2d(+T(1),  -inf), +180 );
  check( Math::atan2d(-T(1),  -inf), -180 );
  check( Math::atan2d(+T(1),  +inf), +0.0 );
  check( Math::atan2d(-T(1),  +inf), -0.0 );
  check( Math::atan2d( +inf, +T(1)),  +90 );
  check( Math::atan2d( +inf, -T(1)),  +90 );
  check( Math::atan2d( -inf, +T(1)),  -90 );
  check( Math::atan2d( -inf, -T(1)),  -90 );
  check( Math::atan2d( +inf,  -inf), +135 );
  check( Math::atan2d( -inf,  -inf), -135 );
  check( Math::atan2d( +inf,  +inf),  +45 );
  check( Math::atan2d( -inf,  +inf),  -45 );
  check( Math::atan2d(  nan, +T(1)),  nan );
  check( Math::atan2d(+T(1),   nan),  nan );

  Math::real e;
  check( Math::sum(+T(9), -T(9), e), +0.0 );
  check( Math::sum(-T(9), +T(9), e), +0.0 );
  check( Math::sum(-T(0), +T(0), e), +0.0 );
  check( Math::sum(+T(0), -T(0), e), +0.0 );
  check( Math::sum(-T(0), -T(0), e), -0.0 );
  check( Math::sum(+T(0), +T(0), e), +0.0 );

  check( Math::AngNormalize(-T(900)), -180 );
  check( Math::AngNormalize(-T(720)), -0.0 );
  check( Math::AngNormalize(-T(540)), -180 );
  check( Math::AngNormalize(-T(360)), -0.0 );
  check( Math::AngNormalize(-T(180)), -180 );
  check( Math::AngNormalize(  -T(0)), -0.0 );
  check( Math::AngNormalize(  +T(0)), +0.0 );
  check( Math::AngNormalize( T(180)), +180 );
  check( Math::AngNormalize( T(360)), +0.0 );
  check( Math::AngNormalize( T(540)), +180 );
  check( Math::AngNormalize( T(720)), +0.0 );
  check( Math::AngNormalize( T(900)), +180 );

  check( Math::AngDiff(+T(  0), +T(  0), e), +0.0 );
  check( Math::AngDiff(+T(  0), -T(  0), e), -0.0 );
  check( Math::AngDiff(-T(  0), +T(  0), e), +0.0 );
  check( Math::AngDiff(-T(  0), -T(  0), e), +0.0 );
  check( Math::AngDiff(+T(  5), +T(365), e), +0.0 );
  check( Math::AngDiff(+T(365), +T(  5), e), -0.0 );
  check( Math::AngDiff(+T(  5), +T(185), e), +180.0 );
  check( Math::AngDiff(+T(185), +T(  5), e), -180.0 );
  check( Math::AngDiff( +eps  , +T(180), e), +180.0 );
  check( Math::AngDiff( -eps  , +T(180), e), -180.0 );
  check( Math::AngDiff( +eps  , -T(180), e), +180.0 );
  check( Math::AngDiff( -eps  , -T(180), e), -180.0 );

  check( Utility::val<T>("+0"), +0.0 );
  check( Utility::val<T>("-0"), -0.0 );
  check( Utility::val<T>("nan"), nan );
  check( Utility::val<T>("+inf"), +inf );
  check( Utility::val<T>( "inf"), +inf );
  check( Utility::val<T>("-inf"), -inf );

  strcheck( Utility::str<T>( nan, 0),  "nan" );
  strcheck( Utility::str<T>(-inf, 0), "-inf" );
  strcheck( Utility::str<T>(-T(3.5), 0),   "-4" );
  strcheck( Utility::str<T>(-T(2.5), 0),   "-2" );
  strcheck( Utility::str<T>(-T(1.5), 0),   "-2" );
  strcheck( Utility::str<T>(-T(0.5), 0),   "-0" );
  strcheck( Utility::str<T>(-T(0  ), 0),   "-0" );
  strcheck( Utility::str<T>(+T(0  ), 0),    "0" );
  strcheck( Utility::str<T>(+T(0.5), 0),    "0" );
  strcheck( Utility::str<T>(+T(1.5), 0),    "2" );
  strcheck( Utility::str<T>(+T(2.5), 0),    "2" );
  strcheck( Utility::str<T>(+T(3.5), 0),    "4" );
  strcheck( Utility::str<T>(+inf, 0),  "inf" );

  DMS::flag ind;
  check( DMS::Decode("+0", ind),  +0.0 );
  check( DMS::Decode("-0", ind),  -0.0 );
  check( DMS::Decode("nan", ind),  nan );
  check( DMS::Decode("+inf", ind),  +inf );
  check( DMS::Decode( "inf", ind),  +inf );
  check( DMS::Decode("-inf", ind),  -inf );
  check( DMS::Decode("+0N", ind),  +0.0 );
  check( DMS::Decode("-0N", ind),  -0.0 );
  check( DMS::Decode("+0S", ind),  -0.0 );
  check( DMS::Decode("-0S", ind),  +0.0 );

  if (n) {
    cout << n << " failure" << (n > 1 ? "s" : "") << "\n";
    return 1;
  }
}
