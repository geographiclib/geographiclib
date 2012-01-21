// Example of using the GeographicLib::Math class
// $Id: 04519bb67e82229a86ee23002ddc27a6b5bb2939 $

#include <iostream>
#include <exception>
#include <GeographicLib/Math.hpp>

using namespace std;
using namespace GeographicLib;

int main() {
  try {
    cout << Math::pi() << " " << Math::sq(Math::pi()) << "\n";
  }
  catch (const exception& e) {
    cerr << "Caught exception: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
