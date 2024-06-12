//! Copyright : Apache 2.0, see LICENSE 
//! 
//! Copyright : see license.txt
//!
//! \brief Basic geometric types
//
#ifndef GEOMOBJECTS_HXX_
#define GEOMOBJECTS_HXX_

#include "../StdHeaders.hxx"
#include "../Geometry/GeomConstants.hxx"

namespace sac_de_billes {

//!
using Identifier = long;

// auxiliary
using PhaseType = long;

// Geometric types
// 1D objects
template<unsigned short DIM>
using DiscPoint = std::array<long, DIM>;

template<unsigned short DIM>
using Point = std::array<double, DIM>;

// checks whether the type T can be understood as a point
template<class T, unsigned short DIM>
constexpr bool is_Point = is_base_of<Point<DIM>, T>::value or is_convertible<Point<DIM>, T>::value;
template<unsigned short DIM>
class Sphere;

} // namespace sac_de_billes

#endif /* GEOMOBJECTS_HXX_ */
