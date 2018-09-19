/* Dan Ibanez: This code is modified from Devon Powell's
 * original r3d code, whose copyright notice is reproduced below:
 * ----
 * Routines for fast, geometrically robust clipping operations
 * and analytic volume/moment computations over polyhedra in 3D.
 *
 * Devon Powell
 * 31 August 2015
 *
 * This program was prepared by Los Alamos National Security, LLC at Los Alamos
 * National Laboratory (LANL) under contract No. DE-AC52-06NA25396 with the U.S.
 * Department of Energy (DOE). All rights in the program are reserved by the DOE
 * and Los Alamos National Security, LLC. Permission is granted to the public to
 * copy and use this software without charge, provided that this Notice and any
 * statement of authorship are reproduced on all copies. Neither the U.S.
 * Government nor LANS makes any warranty, express or implied, or assumes any
 * liability or responsibility for the use of this software.
 * ----
 * We keep our own copy due to two major changes:
 *  (1) annotation to execute on GPUs via CUDA using Kokkos.
 *      this also requires the code to be in a header file,
 *      because Kokkos is incompatible with relocatable device code
 *      in CUDA versions prior to 8.0.
 *  (2) rewritten in modern C++ with templating over
 *      dimension where easily applicable.
 */

#ifndef R3D_HPP
#define R3D_HPP

#include <cfloat>
#include <cmath>
#include <initializer_list>
#include <new>
#include <type_traits>

#ifdef R3D_USE_CUDA
#define R3D_INLINE __device__ __host__ inline
#else
#define R3D_INLINE inline
#endif

/* The MAX_VERTS constants are more than just the maximum number
   of vertices of a representable polytope: the clip() code actually
   adds new vertices prior to removing old ones, so this actually
   needs to be the maximum number of representable vertices
   plus the maximum number of representable edges intersected by a plane */

#ifndef R1D_MAX_VERTS
#define R1D_MAX_VERTS 4
#endif

#ifndef R2D_MAX_VERTS
#define R2D_MAX_VERTS 64
#endif

#ifndef R3D_MAX_VERTS
#define R3D_MAX_VERTS 128
#endif

namespace r3d {

using Real = double;
using Int = int;

constexpr Real ONE_THIRD = (1.0 / 3.0);
constexpr Real ONE_SIXTH = (1.0 / 6.0);

template <typename T>
struct ArithTraits;

template <>
struct ArithTraits<double> {
  static R3D_INLINE double max() { return DBL_MAX; }
  static R3D_INLINE double min() { return -DBL_MAX; }
};

R3D_INLINE Real cube(Real x) { return x * x * x; }

template <typename T>
constexpr R3D_INLINE T square(T x) {
  return x * x;
}

template <typename T, Int n>
class Few {
  using UninitT = typename std::aligned_storage<sizeof(T), alignof(T)>::type;
  UninitT array_[n];

 public:
  enum { size = n };
  R3D_INLINE T* data() { return reinterpret_cast<T*>(array_); }
  R3D_INLINE T const* data() const {
    return reinterpret_cast<T const*>(array_);
  }
  R3D_INLINE T volatile* data() volatile {
    return reinterpret_cast<T volatile*>(array_);
  }
  R3D_INLINE T const volatile* data() const volatile {
    return reinterpret_cast<T const volatile*>(array_);
  }
  R3D_INLINE T& operator[](Int i) { return data()[i]; }
  R3D_INLINE T const& operator[](Int i) const { return data()[i]; }
  R3D_INLINE T volatile& operator[](Int i) volatile { return data()[i]; }
  R3D_INLINE T const volatile& operator[](Int i) const volatile {
    return data()[i];
  }
  Few(std::initializer_list<T> l) {
    Int i = 0;
    for (auto it = l.begin(); it != l.end(); ++it) {
      new (data() + (i++)) T(*it);
    }
  }
  R3D_INLINE Few() {
    for (Int i = 0; i < n; ++i) new (data() + i) T();
  }
  R3D_INLINE ~Few() {
    for (Int i = 0; i < n; ++i) (data()[i]).~T();
  }
  R3D_INLINE void operator=(Few<T, n> const& rhs) volatile {
    for (Int i = 0; i < n; ++i) data()[i] = rhs[i];
  }
  R3D_INLINE void operator=(Few<T, n> const& rhs) {
    for (Int i = 0; i < n; ++i) data()[i] = rhs[i];
  }
  R3D_INLINE void operator=(Few<T, n> const volatile& rhs) {
    for (Int i = 0; i < n; ++i) data()[i] = rhs[i];
  }
  R3D_INLINE Few(Few<T, n> const& rhs) {
    for (Int i = 0; i < n; ++i) new (data() + i) T(rhs[i]);
  }
  R3D_INLINE Few(Few<T, n> const volatile& rhs) {
    for (Int i = 0; i < n; ++i) new (data() + i) T(rhs[i]);
  }
};

template <Int n>
struct Vector : public Few<Real, n> {
  R3D_INLINE Vector() {}
  inline Vector(std::initializer_list<Real> l) : Few<Real, n>(l) {}
  R3D_INLINE void operator=(Vector<n> const& rhs) volatile {
    Few<Real, n>::operator=(rhs);
  }
  R3D_INLINE Vector(Vector<n> const& rhs) : Few<Real, n>(rhs) {}
  R3D_INLINE Vector(const volatile Vector<n>& rhs) : Few<Real, n>(rhs) {}
};

R3D_INLINE Vector<2> vector_2(Real x, Real y) {
  Vector<2> v;
  v[0] = x;
  v[1] = y;
  return v;
}

R3D_INLINE Vector<3> vector_3(Real x, Real y, Real z) {
  Vector<3> v;
  v[0] = x;
  v[1] = y;
  v[2] = z;
  return v;
}

template <Int n>
R3D_INLINE Vector<n> operator+(Vector<n> a, Vector<n> b) {
  Vector<n> c;
  for (Int i = 0; i < n; ++i) c[i] = a[i] + b[i];
  return c;
}

template <Int n>
R3D_INLINE Vector<n>& operator+=(Vector<n>& a, Vector<n> b) {
  a = a + b;
  return a;
}

template <Int n>
R3D_INLINE Vector<n> operator-(Vector<n> a, Vector<n> b) {
  Vector<n> c;
  for (Int i = 0; i < n; ++i) c[i] = a[i] - b[i];
  return c;
}

template <Int n>
R3D_INLINE Vector<n>& operator-=(Vector<n>& a, Vector<n> b) {
  a = a - b;
  return a;
}

template <Int n>
R3D_INLINE Vector<n> operator-(Vector<n> a) {
  Vector<n> c;
  for (Int i = 0; i < n; ++i) c[i] = -a[i];
  return c;
}

template <Int n>
R3D_INLINE Vector<n> operator*(Vector<n> a, Real b) {
  Vector<n> c;
  for (Int i = 0; i < n; ++i) c[i] = a[i] * b;
  return c;
}

template <Int n>
R3D_INLINE Vector<n>& operator*=(Vector<n>& a, Real b) {
  a = a * b;
  return a;
}

template <Int n>
R3D_INLINE Vector<n> operator*(Real a, Vector<n> b) {
  return b * a;
}

template <Int n>
R3D_INLINE Vector<n> operator/(Vector<n> a, Real b) {
  Vector<n> c;
  for (Int i = 0; i < n; ++i) c[i] = a[i] / b;
  return c;
}

template <Int n>
R3D_INLINE Vector<n>& operator/=(Vector<n>& a, Real b) {
  a = a / b;
  return a;
}

template <Int n>
R3D_INLINE Real operator*(Vector<n> a, Vector<n> b) {
  Real c = a[0] * b[0];
  for (Int i = 1; i < n; ++i) c += a[i] * b[i];
  return c;
}

template <Int n>
R3D_INLINE Real norm_squared(Vector<n> v) {
  return v * v;
}

template <Int n>
R3D_INLINE Real norm(Vector<n> v) {
  return std::sqrt(norm_squared(v));
}

template <Int n>
R3D_INLINE Vector<n> normalize(Vector<n> v) {
  return v / norm(v);
}

R3D_INLINE Vector<3> cross(Vector<3> a, Vector<3> b) {
  return vector_3(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
      a[0] * b[1] - a[1] * b[0]);
}

template <Int dim>
struct Plane {
  Vector<dim> n; /*!< Unit-length normal vector. */
  Real d;        /*!< Signed perpendicular distance to the origin. */
};

template <Int dim>
struct Vertex {
  Int pnbrs[dim];  /*!< Neighbor indices. */
  Vector<dim> pos; /*!< Vertex position. */
};

template <Int dim>
struct MaxVerts;

template <>
struct MaxVerts<1> {
  enum { value = R1D_MAX_VERTS };
};

template <>
struct MaxVerts<2> {
  enum { value = R2D_MAX_VERTS };
};

template <>
struct MaxVerts<3> {
  enum { value = R3D_MAX_VERTS };
};

/**
 * \brief A polyhedron. Can be convex, nonconvex, even multiply-connected.
 */
template <Int dim>
struct Polytope {
  enum { max_verts = MaxVerts<dim>::value };
  Vertex<dim> verts[max_verts]; /*!< Vertex buffer. */
  Int nverts;                   /*!< Number of vertices in the buffer. */
};

template <Int dim>
R3D_INLINE Vector<dim> wav(Vector<dim> va, Real wa, Vector<dim> vb, Real wb) {
  Vector<dim> vr;
  for (Int i = 0; i < dim; ++i) vr[i] = (wa * va[i] + wb * vb[i]) / (wa + wb);
  return vr;
}

/* This class encapsulates the code that was different
 * between r3d_clip and r2d_clip
 */

template <Int dim>
struct ClipHelper;

template <>
struct ClipHelper<3> {
  R3D_INLINE static void finalize_links(Int onv, Polytope<3>& poly) {
    for (auto vstart = onv; vstart < poly.nverts; ++vstart) {
      auto vcur = vstart;
      auto vnext = poly.verts[vcur].pnbrs[0];
      do {
        Int np;
        for (np = 0; np < 3; ++np)
          if (poly.verts[vnext].pnbrs[np] == vcur) break;
        vcur = vnext;
        auto pnext = (np + 1) % 3;
        vnext = poly.verts[vcur].pnbrs[pnext];
      } while (vcur < onv);
      poly.verts[vstart].pnbrs[2] = vcur;
      poly.verts[vcur].pnbrs[1] = vstart;
    }
  }
  R3D_INLINE static void init_new_vert_link(
      Polytope<3>& poly, Int vcur, Int np) {
    (void)np;
    poly.verts[poly.nverts].pnbrs[0] = vcur;
  }
};

template <>
struct ClipHelper<2> {
  R3D_INLINE static void finalize_links(Int onv, Polytope<2>& poly) {
    for (auto vstart = onv; vstart < poly.nverts; ++vstart) {
      if (poly.verts[vstart].pnbrs[1] >= 0) continue;
      auto vcur = poly.verts[vstart].pnbrs[0];
      do {
        vcur = poly.verts[vcur].pnbrs[0];
      } while (vcur < onv);
      poly.verts[vstart].pnbrs[1] = vcur;
      poly.verts[vcur].pnbrs[0] = vstart;
    }
  }
  R3D_INLINE static void init_new_vert_link(
      Polytope<2>& poly, Int vcur, Int np) {
    poly.verts[poly.nverts].pnbrs[1 - np] = vcur;
    poly.verts[poly.nverts].pnbrs[np] = -1;
  }
};

template <>
struct ClipHelper<1> {
  R3D_INLINE static void finalize_links(Int, Polytope<1>&) {}
  R3D_INLINE static void init_new_vert_link(
      Polytope<1>& poly, Int vcur, Int np) {
    poly.verts[poly.nverts].pnbrs[np] = vcur;
  }
};

/**
 * \brief Clip a polytope against an arbitrary number of clip planes (find its
 * intersection with a set of half-spaces).
 *
 * \param [in,out] poly
 * The polygon to be clipped.
 *
 * \param [in] planes
 * An array of planes against which to clip this polygon.
 *
 */
template <Int dim>
R3D_INLINE void clip(
    Polytope<dim>& poly, Plane<dim> const* planes, Int nplanes) {
  if (poly.nverts <= 0) return;

  // variable declarations
  Int v, p, np, onv, vcur, vnext, numunclipped;

  // signed distances to the clipping plane
  Real sdists[Polytope<dim>::max_verts];
  Real smin, smax;

  // loop over each clip plane
  for (p = 0; p < nplanes; ++p) {
    // calculate signed distances to the clip plane
    onv = poly.nverts;
    smin = ArithTraits<Real>::max();
    smax = ArithTraits<Real>::min();

    // for marking clipped vertices
    Int clipped[MaxVerts<dim>::value] = {};  // all initialized to zero
    for (v = 0; v < onv; ++v) {
      sdists[v] = planes[p].d + (poly.verts[v].pos * planes[p].n);
      if (sdists[v] < smin) smin = sdists[v];
      if (sdists[v] > smax) smax = sdists[v];
      if (sdists[v] < 0.0) clipped[v] = 1;
    }

    // skip this face if the poly lies entirely on one side of it
    if (smin >= 0.0) continue;
    if (smax <= 0.0) {
      poly.nverts = 0;
      return;
    }

    // check all edges and insert new vertices on the bisected edges
    for (vcur = 0; vcur < onv; ++vcur) {
      if (clipped[vcur]) continue;
      for (np = 0; np < dim; ++np) {
        vnext = poly.verts[vcur].pnbrs[np];
        if (!clipped[vnext]) continue;
        ClipHelper<dim>::init_new_vert_link(poly, vcur, np);
        poly.verts[vcur].pnbrs[np] = poly.nverts;
        poly.verts[poly.nverts].pos = wav(poly.verts[vcur].pos, -sdists[vnext],
            poly.verts[vnext].pos, sdists[vcur]);
        ++(poly.nverts);
      }
    }

    // for each new vert, search around the poly for its new neighbors
    // and doubly-link everything
    ClipHelper<dim>::finalize_links(onv, poly);

    // go through and compress the vertex list, removing clipped verts
    // and re-indexing accordingly (reusing `clipped` to re-index everything)
    numunclipped = 0;
    for (v = 0; v < poly.nverts; ++v) {
      if (!clipped[v]) {
        poly.verts[numunclipped] = poly.verts[v];
        clipped[v] = numunclipped++;
      }
    }
    poly.nverts = numunclipped;
    for (v = 0; v < poly.nverts; ++v)
      for (np = 0; np < dim; ++np)
        poly.verts[v].pnbrs[np] = clipped[poly.verts[v].pnbrs[np]];
  }
}

template <Int dim>
R3D_INLINE void clip(Polytope<dim>& poly, Plane<dim> const& plane) {
  clip(poly, &plane, 1);
}

template <Int dim, Int nplanes>
R3D_INLINE void clip(Polytope<dim>& poly, Few<Plane<dim>, nplanes> planes) {
  clip(poly, &planes[0], nplanes);
}

/**
 * \brief Initialize a polyhedron as a tetrahedron.
 *
 * \param [in] verts
 * An array of four vectors, giving the vertices of the tetrahedron.
 *
 */
R3D_INLINE void init(Polytope<3>& poly, Few<Vector<3>, 4> verts) {
  // initialize graph connectivity
  poly.nverts = 4;
  poly.verts[0].pnbrs[0] = 1;
  poly.verts[0].pnbrs[1] = 3;
  poly.verts[0].pnbrs[2] = 2;
  poly.verts[1].pnbrs[0] = 2;
  poly.verts[1].pnbrs[1] = 3;
  poly.verts[1].pnbrs[2] = 0;
  poly.verts[2].pnbrs[0] = 0;
  poly.verts[2].pnbrs[1] = 3;
  poly.verts[2].pnbrs[2] = 1;
  poly.verts[3].pnbrs[0] = 1;
  poly.verts[3].pnbrs[1] = 2;
  poly.verts[3].pnbrs[2] = 0;

  // copy vertex coordinates
  for (Int v = 0; v < 4; ++v) poly.verts[v].pos = verts[v];
}

/**
 * \brief Initialize a triangle from a list of vertices.
 *
 * Array of length `numverts` giving the vertices of the input polygon, in
 * counterclockwise order.
 *
 */
R3D_INLINE void init(Polytope<2>& poly, Few<Vector<2>, 3> vertices) {
  constexpr Int numverts = 3;
  // init the poly
  poly.nverts = numverts;
  for (Int v = 0; v < poly.nverts; ++v) {
    poly.verts[v].pos = vertices[v];
    poly.verts[v].pnbrs[0] = (v + 1) % (numverts);
    poly.verts[v].pnbrs[1] = (numverts + v - 1) % (numverts);
  }
}

/**
 * \brief Initialize an edge from a list of vertices.
 */
R3D_INLINE void init(Polytope<1>& poly, Few<Vector<1>, 2> vertices) {
  constexpr Int numverts = 2;
  // init the poly
  poly.nverts = numverts;
  for (Int v = 0; v < poly.nverts; ++v) {
    poly.verts[v].pos = vertices[v];
    poly.verts[v].pnbrs[0] = 1 - v;
  }
}

R3D_INLINE Plane<3> tet_face_from_verts(Vector<3> a, Vector<3> b, Vector<3> c) {
  auto center = ONE_THIRD * (a + b + c);
  auto normal = normalize(cross((b - a), (c - a)));
  auto d = -(normal * center);
  return Plane<3>{normal, d};
}

/**
 * \brief Get four faces (unit normals and distances to the origin)
 * from a four-vertex description of a tetrahedron.
 *
 * \returns
 * Array of four planes defining the faces of the tetrahedron.
 *
 * \param [in] verts
 * Array of four vectors defining the vertices of the tetrahedron.
 *
 */
R3D_INLINE Few<Plane<3>, 4> faces_from_verts(Few<Vector<3>, 4> verts) {
  Few<Plane<3>, 4> faces;
  faces[0] = tet_face_from_verts(verts[3], verts[2], verts[1]);
  faces[1] = tet_face_from_verts(verts[0], verts[2], verts[3]);
  faces[2] = tet_face_from_verts(verts[0], verts[3], verts[1]);
  faces[3] = tet_face_from_verts(verts[0], verts[1], verts[2]);
  return faces;
}

/**
 * \brief Get all faces (unit normals and distances to the origin)
 * from a triangle.
 *
 * \returns
 * Array of planes of length `numverts` defining the faces of the polygon.
 *
 * \param [in] vertices
 * Array of length `numverts` giving the vertices of the input triangle, in
 * counterclockwise order.
 *
 */
R3D_INLINE Few<Plane<2>, 3> faces_from_verts(Few<Vector<2>, 3> vertices) {
  constexpr Int numverts = 3;
  Int f;
  Vector<2> p0, p1;
  Few<Plane<2>, 3> faces;
  // calculate a centroid and a unit normal for each face
  for (f = 0; f < numverts; ++f) {
    p0 = vertices[f];
    p1 = vertices[(f + 1) % numverts];
    // normal of the edge
    faces[f].n[0] = p0[1] - p1[1];
    faces[f].n[1] = p1[0] - p0[0];
    // normalize the normals and set the signed distance to origin
    faces[f].n = normalize(faces[f].n);
    faces[f].d = -(faces[f].n * p0);
  }
  return faces;
}

/**
 * \brief Get all faces (unit normals and distances to the origin)
 * from an edge
 *
 * \param [in] vertices
 * Array giving the vertices of the input edge
 *
 */
R3D_INLINE Few<Plane<1>, 2> faces_from_verts(Few<Vector<1>, 2> vertices) {
  Few<Plane<1>, 2> faces;
  faces[0].n = normalize(vertices[1] - vertices[0]);
  faces[1].n = normalize(vertices[0] - vertices[1]);
  faces[0].d = -(faces[0].n * vertices[0]);
  faces[1].d = -(faces[1].n * vertices[1]);
  return faces;
}

R3D_INLINE constexpr Int num_moments_3d(Int order) {
  return ((order + 1) * (order + 2) * (order + 3) / 6);
}

R3D_INLINE constexpr Int num_moments_2d(Int order) {
  return ((order + 1) * (order + 2) / 2);
}

template <Int dim, Int order>
struct NumMoments;

template <Int order>
struct NumMoments<3, order> {
  enum { value = num_moments_3d(order) };
};

template <Int order>
struct NumMoments<2, order> {
  enum { value = num_moments_2d(order) };
};

/**
 * \brief Integrate a polynomial density over a polyhedron using simplicial
 * decomposition.
 *
 * \details Uses the fast recursive method of Koehl (2012) to carry out the
 * integration.
 * The template parameter is the order of the polynomial density field.
 * 0 for constant (1 moment), 1 for linear (4 moments),
 * 2 for quadratic (10 moments), etc.
 *
 * \param [in] poly
 * The polyhedron over which to integrate.
 *
 * \param [in,out] moments
 * Array to be filled with the integration results, up to the specified
 * `polyorder`. Must be at least
 * `(polyorder+1)*(polyorder+2)*(polyorder+3)/6` long.
 * A conventient macro, `num_moments_3d()` is provided to compute the number
 * of moments for a given order.
 * Order of moments is row-major, i.e. `1`, `x`, `y`, `z`, `x^2`,
 * `x*y`, `x*z`, `y^2`, `y*z`, `z^2`, `x^3`, `x^2*y`...
 */

template <Int polyorder>
R3D_INLINE void reduce(Polytope<3> const& poly, Real* moments) {
  if (poly.nverts <= 0) return;

  // var declarations
  Real sixv;
  Int np, m, i, j, k, corder;
  Int vstart, pstart, vcur, vnext, pnext;
  Vector<3> v0, v1, v2;

  // zero the moments
  for (m = 0; m < num_moments_3d(polyorder); ++m) moments[m] = 0.0;

  // for keeping track of which edges have been visited
  Int emarks[Polytope<3>::max_verts][3] = {{}};  // initialized to zero

  // Storage for coefficients
  // keep two layers of the pyramid of coefficients
  // Note: Uses twice as much space as needed, but indexing is faster this way
  Int prevlayer = 0;
  Int curlayer = 1;
  Real S[polyorder + 1][polyorder + 1][2];
  Real D[polyorder + 1][polyorder + 1][2];
  Real C[polyorder + 1][polyorder + 1][2];

  // loop over all vertices to find the starting point for each face
  for (vstart = 0; vstart < poly.nverts; ++vstart)
    for (pstart = 0; pstart < 3; ++pstart) {
      // skip this face if we have marked it
      if (emarks[vstart][pstart]) continue;

      // initialize face looping
      pnext = pstart;
      vcur = vstart;
      emarks[vcur][pnext] = 1;
      vnext = poly.verts[vcur].pnbrs[pnext];
      v0 = poly.verts[vcur].pos;

      // move to the second edge
      for (np = 0; np < 3; ++np)
        if (poly.verts[vnext].pnbrs[np] == vcur) break;
      vcur = vnext;
      pnext = (np + 1) % 3;
      emarks[vcur][pnext] = 1;
      vnext = poly.verts[vcur].pnbrs[pnext];

      // make a triangle fan using edges
      // and first vertex
      while (vnext != vstart) {
        v2 = poly.verts[vcur].pos;
        v1 = poly.verts[vnext].pos;

        sixv = (-v2[0] * v1[1] * v0[2] + v1[0] * v2[1] * v0[2] +
                v2[0] * v0[1] * v1[2] - v0[0] * v2[1] * v1[2] -
                v1[0] * v0[1] * v2[2] + v0[0] * v1[1] * v2[2]);

        // calculate the moments
        // using the fast recursive method of Koehl (2012)
        // essentially building a set of trinomial pyramids, one layer at a time

        // base case
        S[0][0][prevlayer] = 1.0;
        D[0][0][prevlayer] = 1.0;
        C[0][0][prevlayer] = 1.0;
        moments[0] += ONE_SIXTH * sixv;

        // build up successive polynomial orders
        for (corder = 1, m = 1; corder <= polyorder; ++corder) {
          for (i = corder; i >= 0; --i)
            for (j = corder - i; j >= 0; --j, ++m) {
              k = corder - i - j;
              C[i][j][curlayer] = 0;
              D[i][j][curlayer] = 0;
              S[i][j][curlayer] = 0;
              if (i > 0) {
                C[i][j][curlayer] += v2[0] * C[i - 1][j][prevlayer];
                D[i][j][curlayer] += v1[0] * D[i - 1][j][prevlayer];
                S[i][j][curlayer] += v0[0] * S[i - 1][j][prevlayer];
              }
              if (j > 0) {
                C[i][j][curlayer] += v2[1] * C[i][j - 1][prevlayer];
                D[i][j][curlayer] += v1[1] * D[i][j - 1][prevlayer];
                S[i][j][curlayer] += v0[1] * S[i][j - 1][prevlayer];
              }
              if (k > 0) {
                C[i][j][curlayer] += v2[2] * C[i][j][prevlayer];
                D[i][j][curlayer] += v1[2] * D[i][j][prevlayer];
                S[i][j][curlayer] += v0[2] * S[i][j][prevlayer];
              }
              D[i][j][curlayer] += C[i][j][curlayer];
              S[i][j][curlayer] += D[i][j][curlayer];
              moments[m] += sixv * S[i][j][curlayer];
            }
          curlayer = 1 - curlayer;
          prevlayer = 1 - prevlayer;
        }

        // move to the next edge
        for (np = 0; np < 3; ++np)
          if (poly.verts[vnext].pnbrs[np] == vcur) break;
        vcur = vnext;
        pnext = (np + 1) % 3;
        emarks[vcur][pnext] = 1;
        vnext = poly.verts[vcur].pnbrs[pnext];
      }
    }

  // reuse C to recursively compute the leading multinomial coefficients
  C[0][0][prevlayer] = 1.0;
  for (corder = 1, m = 1; corder <= polyorder; ++corder) {
    for (i = corder; i >= 0; --i)
      for (j = corder - i; j >= 0; --j, ++m) {
        k = corder - i - j;
        C[i][j][curlayer] = 0.0;
        if (i > 0) C[i][j][curlayer] += C[i - 1][j][prevlayer];
        if (j > 0) C[i][j][curlayer] += C[i][j - 1][prevlayer];
        if (k > 0) C[i][j][curlayer] += C[i][j][prevlayer];
        moments[m] /=
            C[i][j][curlayer] * (corder + 1) * (corder + 2) * (corder + 3);
      }
    curlayer = 1 - curlayer;
    prevlayer = 1 - prevlayer;
  }
}

/**
 * \brief Integrate a polynomial density over a polygon using simplicial
 * decomposition.
 *
 * \details Uses the fast recursive method of Koehl (2012) to carry out the
 * integration.
 * The template parameter is the order of the polynomial density field.
 * 0 for constant (1 moment), 1 for linear (3 moments),
 * 2 for quadratic (6 moments), etc.
 *
 * \param [in] poly
 * The polygon over which to integrate.
 *
 * \param [in,out] moments
 * Array to be filled with the integration results, up to the specified
 * `polyorder`. Must be at least `(polyorder+1)*(polyorder+2)/2` long.
 * A conventient macro, `R2D_NUM_MOMENTS()` is provided to compute the
 * number of moments for a given order.
 * Order of moments is row-major, i.e. `1`, `x`, `y`, `x^2`, `x*y`, `y^2`,
 * `x^3`, `x^2*y`...
 *
 */
template <Int polyorder>
R3D_INLINE void reduce(Polytope<2> const& poly, Real* moments) {
  if (poly.nverts <= 0) return;

  // var declarations
  Int vcur, vnext, m, i, j, corder;
  Real twoa;
  Vector<2> v0, v1;

  // zero the moments
  for (m = 0; m < num_moments_2d(polyorder); ++m) moments[m] = 0.0;

  // Storage for coefficients
  // keep two layers of the triangle of coefficients
  Int prevlayer = 0;
  Int curlayer = 1;
  Real D[polyorder + 1][2];
  Real C[polyorder + 1][2];

  // iterate over edges and compute a sum over simplices
  for (vcur = 0; vcur < poly.nverts; ++vcur) {
    vnext = poly.verts[vcur].pnbrs[0];
    v0 = poly.verts[vcur].pos;
    v1 = poly.verts[vnext].pos;
    twoa = (v0[0] * v1[1] - v0[1] * v1[0]);

    // calculate the moments
    // using the fast recursive method of Koehl (2012)
    // essentially building a set of Pascal's triangles, one layer at a time

    // base case
    D[0][prevlayer] = 1.0;
    C[0][prevlayer] = 1.0;
    moments[0] += 0.5 * twoa;

    // build up successive polynomial orders
    for (corder = 1, m = 1; corder <= polyorder; ++corder) {
      for (i = corder; i >= 0; --i, ++m) {
        j = corder - i;
        C[i][curlayer] = 0;
        D[i][curlayer] = 0;
        if (i > 0) {
          C[i][curlayer] += v1[0] * C[i - 1][prevlayer];
          D[i][curlayer] += v0[0] * D[i - 1][prevlayer];
        }
        if (j > 0) {
          C[i][curlayer] += v1[1] * C[i][prevlayer];
          D[i][curlayer] += v0[1] * D[i][prevlayer];
        }
        D[i][curlayer] += C[i][curlayer];
        moments[m] += twoa * D[i][curlayer];
      }
      curlayer = 1 - curlayer;
      prevlayer = 1 - prevlayer;
    }
  }

  // reuse C to recursively compute the leading multinomial coefficients
  C[0][prevlayer] = 1.0;
  for (corder = 1, m = 1; corder <= polyorder; ++corder) {
    for (i = corder; i >= 0; --i, ++m) {
      j = corder - i;
      C[i][curlayer] = 0.0;
      if (i > 0) C[i][curlayer] += C[i - 1][prevlayer];
      if (j > 0) C[i][curlayer] += C[i][prevlayer];
      moments[m] /= C[i][curlayer] * (corder + 1) * (corder + 2);
    }
    curlayer = 1 - curlayer;
    prevlayer = 1 - prevlayer;
  }
}

template <Int dim, Int order>
struct Polynomial {
  enum { nterms = NumMoments<dim, order>::value };
  Real coeffs[nterms];
};

template <Int dim, Int order>
R3D_INLINE Real integrate(
    Polytope<dim> const& polytope, Polynomial<dim, order> polynomial) {
  Real moments[decltype(polynomial)::nterms] = {};
  reduce<order>(polytope, moments);
  Real result = 0;
  for (Int i = 0; i < decltype(polynomial)::nterms; ++i)
    result += moments[i] * polynomial.coeffs[i];
  return result;
}

/* Cop-out: instead of figuring out how to implement
   integrate() for 1D edges, just hardcode the
   length measurement */
R3D_INLINE Real measure(Polytope<1> const& polytope) {
  return std::abs(polytope.verts[1].pos[0] - polytope.verts[0].pos[0]);
}

R3D_INLINE Real measure(Polytope<2> const& polytope) {
  return integrate(polytope, Polynomial<2, 0>{{1}});
}

R3D_INLINE Real measure(Polytope<3> const& polytope) {
  return integrate(polytope, Polynomial<3, 0>{{1}});
}

template <Int dim>
R3D_INLINE void intersect_simplices(Polytope<dim>& poly,
    Few<Vector<dim>, dim + 1> verts0, Few<Vector<dim>, dim + 1> verts1) {
  init(poly, verts0);
  auto faces1 = faces_from_verts(verts1);
  clip(poly, faces1);
}

/**
 * \brief Initialize a general polyhedron from a full boundary description.
 *
 * \param [out] poly
 * The polyhedron to initialize.
 *
 * \param [in] vertices
 * Array of length `numverts` giving the vertices of the input polyhedron.
 *
 * \param [in] numverts
 * Number of vertices in the input polyhedron.
 *
 * \param [in] faceinds
 * Connectivity array, giving the indices of vertices in the
 * order they appear around each face of the input polyhedron.
 *
 * \param [in] numvertsperface
 * An array of length `numfaces` giving the number of vertices for each face
 * of the input polyhedron.
 *
 * \param [in] numfaces
 * Number of faces in the input polyhedron.
 *
 */
R3D_INLINE void init_poly(Polytope<3>& poly, Vector<3>* vertices, Int numverts,
    Int** faceinds, Int* numvertsperface, Int numfaces) {
  // dummy vars
  Int v, vprev, vcur, vnext, f, np;

  // count up the number of faces per vertex
  // and act accordingly
  Int eperv[R3D_MAX_VERTS] = {0};
  Int maxeperv = 0;
  for (f = 0; f < numfaces; ++f)
    for (v = 0; v < numvertsperface[f]; ++v) ++eperv[faceinds[f][v]];
  for (v = 0; v < numverts; ++v)
    if (eperv[v] > maxeperv) maxeperv = eperv[v];

  // clear the poly
  poly.nverts = 0;

  if (maxeperv == 3) {
    // simple case with no need for duplicate vertices

    // read in vertex locations
    poly.nverts = numverts;
    for (v = 0; v < poly.nverts; ++v) {
      poly.verts[v].pos = vertices[v];
      for (np = 0; np < 3; ++np) poly.verts[v].pnbrs[np] = R3D_MAX_VERTS;
    }

    // build graph connectivity by correctly orienting half-edges for each
    // vertex
    for (f = 0; f < numfaces; ++f) {
      for (v = 0; v < numvertsperface[f]; ++v) {
        vprev = faceinds[f][v];
        vcur = faceinds[f][(v + 1) % numvertsperface[f]];
        vnext = faceinds[f][(v + 2) % numvertsperface[f]];
        for (np = 0; np < 3; ++np) {
          if (poly.verts[vcur].pnbrs[np] == vprev) {
            poly.verts[vcur].pnbrs[(np + 2) % 3] = vnext;
            break;
          } else if (poly.verts[vcur].pnbrs[np] == vnext) {
            poly.verts[vcur].pnbrs[(np + 1) % 3] = vprev;
            break;
          }
        }
        if (np == 3) {
          poly.verts[vcur].pnbrs[1] = vprev;
          poly.verts[vcur].pnbrs[0] = vnext;
        }
      }
    }
  } else {
    // we need to create duplicate, degenerate vertices to account for more than
    // three edges per vertex. This is complicated.

    Int tface = 0;
    for (v = 0; v < numverts; ++v) tface += eperv[v];

    // need more variables
    Int v0, v1, v00, v11, numunclipped;

    // we need a few extra buffers to handle the necessary operations
    Vertex<3> vbtmp[3 * R3D_MAX_VERTS];
    Int vstart[R3D_MAX_VERTS];

    // build vertex mappings to degenerate duplicates
    // and read in vertex locations
    poly.nverts = 0;
    for (v = 0; v < numverts; ++v) {
      vstart[v] = poly.nverts;
      for (vcur = 0; vcur < eperv[v]; ++vcur) {
        vbtmp[poly.nverts].pos = vertices[v];
        for (np = 0; np < 3; ++np) vbtmp[poly.nverts].pnbrs[np] = R3D_MAX_VERTS;
        ++(poly.nverts);
      }
    }

    // fill in connectivity for all duplicates
    {
      Int util[3 * R3D_MAX_VERTS] = {0};
      for (f = 0; f < numfaces; ++f) {
        for (v = 0; v < numvertsperface[f]; ++v) {
          vprev = faceinds[f][v];
          vcur = faceinds[f][(v + 1) % numvertsperface[f]];
          vnext = faceinds[f][(v + 2) % numvertsperface[f]];
          auto vcur_old = vcur;
          vcur = vstart[vcur] + util[vcur];
          util[vcur_old]++;
          vbtmp[vcur].pnbrs[1] = vnext;
          vbtmp[vcur].pnbrs[2] = vprev;
        }
      }
    }

    // link degenerate duplicates, putting them in the correct order
    // use util to mark and avoid double-processing verts
    {
      Int util[3 * R3D_MAX_VERTS] = {0};
      for (v = 0; v < numverts; ++v) {
        for (v0 = vstart[v]; v0 < vstart[v] + eperv[v]; ++v0) {
          for (v1 = vstart[v]; v1 < vstart[v] + eperv[v]; ++v1) {
            if (vbtmp[v0].pnbrs[2] == vbtmp[v1].pnbrs[1] && !util[v0]) {
              vbtmp[v0].pnbrs[2] = v1;
              vbtmp[v1].pnbrs[0] = v0;
              util[v0] = 1;
            }
          }
        }
      }
    }

    // complete vertex pairs
    {
      Int util[3 * R3D_MAX_VERTS] = {0};
      for (v0 = 0; v0 < numverts; ++v0)
        for (v1 = v0 + 1; v1 < numverts; ++v1) {
          for (v00 = vstart[v0]; v00 < vstart[v0] + eperv[v0]; ++v00)
            for (v11 = vstart[v1]; v11 < vstart[v1] + eperv[v1]; ++v11) {
              if (vbtmp[v00].pnbrs[1] == v1 && vbtmp[v11].pnbrs[1] == v0 &&
                  !util[v00] && !util[v11]) {
                vbtmp[v00].pnbrs[1] = v11;
                vbtmp[v11].pnbrs[1] = v00;
                util[v00] = 1;
                util[v11] = 1;
              }
            }
        }
    }

    // remove unnecessary dummy vertices
    {
      Int util[3 * R3D_MAX_VERTS] = {0};
      for (v = 0; v < numverts; ++v) {
        v0 = vstart[v];
        v1 = vbtmp[v0].pnbrs[0];
        v00 = vbtmp[v0].pnbrs[2];
        v11 = vbtmp[v1].pnbrs[0];
        vbtmp[v00].pnbrs[0] = vbtmp[v0].pnbrs[1];
        vbtmp[v11].pnbrs[2] = vbtmp[v1].pnbrs[1];
        for (np = 0; np < 3; ++np)
          if (vbtmp[vbtmp[v0].pnbrs[1]].pnbrs[np] == v0) break;
        vbtmp[vbtmp[v0].pnbrs[1]].pnbrs[np] = v00;
        for (np = 0; np < 3; ++np)
          if (vbtmp[vbtmp[v1].pnbrs[1]].pnbrs[np] == v1) break;
        vbtmp[vbtmp[v1].pnbrs[1]].pnbrs[np] = v11;
        util[v0] = 1;
        util[v1] = 1;
      }

      // copy to the real vertbuffer and compress
      numunclipped = 0;
      for (v = 0; v < poly.nverts; ++v) {
        if (!util[v]) {
          poly.verts[numunclipped] = vbtmp[v];
          util[v] = numunclipped++;
        }
      }
      poly.nverts = numunclipped;
      for (v = 0; v < poly.nverts; ++v)
        for (np = 0; np < 3; ++np)
          poly.verts[v].pnbrs[np] = util[poly.verts[v].pnbrs[np]];
    }
  }
}

/**
 * \brief Initialize a (simply-connected) general polygon from a list of
 * vertices.
 *
 * \param [out] poly
 * The polygon to initialize.
 *
 * \param [in] vertices
 * Array of length `numverts` giving the vertices of the input polygon, in
 * counterclockwise order.
 *
 * \param [in] numverts
 * Number of vertices in the input polygon.
 *
 */
R3D_INLINE void init_poly(
    Polytope<2>& poly, Vector<2>* vertices, Int numverts) {
  poly.nverts = numverts;
  Int v;
  for (v = 0; v < poly.nverts; ++v) {
    poly.verts[v].pos = vertices[v];
    poly.verts[v].pnbrs[0] = (v + 1) % (poly.nverts);
    poly.verts[v].pnbrs[1] = (poly.nverts + v - 1) % (poly.nverts);
  }
}

/**
 * \brief Get all faces (unit normals and distances to the origin)
 * from a full boundary description of a polyhedron.
 *
 * \param [out] faces
 * Array of planes of length `numfaces` defining the faces of the polyhedron.
 *
 * \param [in] vertices
 * Array of length `numverts` giving the vertices of the input polyhedron.
 *
 * \param [in] faceinds
 * Connectivity array, giving the indices of vertices in the
 * order they appear around each face of the input polyhedron.
 *
 * \param [in] numvertsperface
 * An array of length `numfaces` giving the number of vertices for each face
 * of the input polyhedron.
 *
 * \param [in] numfaces
 * Number of faces in the input polyhedron.
 *
 */
R3D_INLINE void poly_faces_from_verts(Plane<3>* faces, Vector<3>* vertices,
    Int** faceinds, Int* numvertsperface, Int numfaces) {
  // dummy vars
  Int v, f;
  Vector<3> p0, p1, p2;

  // calculate a centroid and a unit normal for each face
  for (f = 0; f < numfaces; ++f) {
    auto centroid = vector_3(0, 0, 0);
    faces[f].n = vector_3(0, 0, 0);

    for (v = 0; v < numvertsperface[f]; ++v) {
      // add cross product of edges to the total normal
      p0 = vertices[faceinds[f][v]];
      p1 = vertices[faceinds[f][(v + 1) % numvertsperface[f]]];
      p2 = vertices[faceinds[f][(v + 2) % numvertsperface[f]]];
      faces[f].n = faces[f].n + cross(p1 - p0, p2 - p0);

      // add the vertex position to the centroid
      centroid = centroid + p0;
    }

    // normalize the normals and set the signed distance to origin
    centroid = centroid / Real(numvertsperface[f]);
    faces[f].n = normalize(faces[f].n);
    faces[f].d = -(faces[f].n * centroid);
  }
}

}  // end namespace r3d

#endif
