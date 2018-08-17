/****************************************************************************
          File:  pyro_types.h
       Created:  Oct 22, 2011

   Description:  Generic, empty "trait" types

*****************************************************************************/


#ifndef pyro_types_h__GUARD
#define pyro_types_h__GUARD

// Waves, distributions, filters, interpolators - just all-purpose,
// frequently-occurring identifiers that take on different meanings 
// in different contexts
struct   Constant    {}
#define  Uniform     Constant
struct   Linear      {}
struct   Triangle    {}
#define  Tri         Triangle
struct   Square      {}
struct   Sawtooth    {}
#define  Saw         Sawtooth
struct   Sine        {}
#define  Sin         Sine
struct   Cosine      {}
#define  Cos         Cosine
struct   Cubic       {}
#define  Smooth      Cubic
struct   Quadratic   {}
#define  Quad        Quadratic
struct   Quintic     {}

struct   Gaussian    {}
#define  Gauss       Gaussian
#define  Bartlett    Triangle
struct   Sinc        {}


// Geometric shapes, domains, spaces, coordinate systems
struct   Line        {}
struct   Circle      {}
struct   Disk        {}
struct   Cube        {}
#define  Box         Cube
struct   Sphere      {}
struct   Hemisphere  {}
struct   Euclidean   {}
#define  Rectilinear Euclidean
struct   Polar       {}
#define  Spherical   Sphere

// Orientation
struct XAxis {}
struct YAxis {}
struct ZAxis {}

   
#endif // End pyro_types_h

