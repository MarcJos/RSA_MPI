//! Copyright : Apache 2.0, see LICENSE 
//! 
//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef SPHERE_HXX_
#define SPHERE_HXX_

#include "../StdHeaders.hxx"

#include "../Geometry/GeomTypes.hxx"


namespace sac_de_billes {
using namespace std;

template<unsigned short DIM>
class Sphere {
public:
        Sphere() :
                center{ create_array<DIM>(0.) }, radius{ 0 }, phase{ 0 } {}
        Sphere(const Point<DIM>& center_, double radius_, PhaseType phase_) :
                center(center_), radius(radius_), phase{ phase_ } {}
        Point<DIM> center;
        double radius;
        PhaseType phase;
        //! \return the volume
        double volume() const;
        //! is this point inside?
        bool isInside(const Point<DIM>& x) const;
};

using SimpleSphereFormat = array<double, 5>;

//! @return : whether two spheres are intersected
//! @tparam DIM : dimension of the ambiant space
//! @param sph1, sph2 : spheres
template<unsigned short DIM>
bool areIntersected(const Sphere<DIM>& sph1, const Sphere<DIM>& sph2);

namespace sphereTools {
//! \return the volume of a ball
template<unsigned short DIM>
double volumeSphere(const double& radius);

//! \return the intersection volume of 2 spheres
//! \param R, r : radii
//! \param d : distance between spheres
template<unsigned short DIM>
double volInter(const double& R, const double& r, const double& d);

//! from 1 format to the other
template<unsigned short DIM>
inline vector<SimpleSphereFormat> vector_fromSphere2Array(
        const vector<Sphere<DIM>>& vec_sphere);
//! from 1 format to the other
template<unsigned short DIM>
inline vector<Sphere<DIM>> vector_fromArray2Sphere(
        const vector<SimpleSphereFormat>& vec_array);
//! from 1 format to the other
template<unsigned short DIM>
SimpleSphereFormat fromSphere2Array(const Sphere<DIM>&);
//! from 1 format to the other
template<unsigned short DIM>
Sphere<DIM> fromArray2Sphere(const SimpleSphereFormat&);

//! sort the spheres by phase
template<unsigned short DIM>
void sort(vector<Sphere<DIM>>&);

//! reads from a line in a file
template<unsigned short DIM>
bool fromLine(istream& fileStream, size_t phase, Sphere<DIM>& sphere);
} // namespace sphereTools
} // namespace sac_de_billes

#include "../Geometry/Sphere.ixx"
#endif /* SPHERE_HXX_ */
