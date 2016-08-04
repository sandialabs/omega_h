#ifndef OMEGA_H_R3D_HPP
#define OMEGA_H_R3D_HPP

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
 */

#include "omega_h_math.hpp"

namespace osh {

namespace r3d {

constexpr Real ONE_THIRD = (1.0 / 3.0);
constexpr Real ONE_SIXTH = (1.0 / 6.0);

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
struct MaxVerts<2> {
  enum { value = 64 };
};

template <>
struct MaxVerts<3> {
  enum { value = 128 };
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
OSH_INLINE Vector<dim> wav(Vector<dim> va, Real wa, Vector<dim> vb, Real wb) {
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
  static void relink(Int onv, Int* nverts, Vertex<3>* vertbuffer) {
    for (auto vstart = onv; vstart < *nverts; ++vstart) {
      auto vcur = vstart;
      auto vnext = vertbuffer[vcur].pnbrs[0];
      do {
        Int np;
        for (np = 0; np < 3; ++np)
          if (vertbuffer[vnext].pnbrs[np] == vcur) break;
        vcur = vnext;
        auto pnext = (np + 1) % 3;
        vnext = vertbuffer[vcur].pnbrs[pnext];
      } while (vcur < onv);
      vertbuffer[vstart].pnbrs[2] = vcur;
      vertbuffer[vcur].pnbrs[1] = vstart;
    }
  }
  static void links_at_nverts(Int* nverts, Vertex<3>* vertbuffer, Int vcur,
                              Int np) {
    (void)np;
    vertbuffer[*nverts].pnbrs[0] = vcur;
  }
};

template <>
struct ClipHelper<2> {
  static void relink(Int onv, Int* nverts, Vertex<2>* vertbuffer) {
    for (auto vstart = onv; vstart < *nverts; ++vstart) {
      if (vertbuffer[vstart].pnbrs[1] >= 0) continue;
      auto vcur = vertbuffer[vstart].pnbrs[0];
      do {
        vcur = vertbuffer[vcur].pnbrs[0];
      } while (vcur < onv);
      vertbuffer[vstart].pnbrs[1] = vcur;
      vertbuffer[vcur].pnbrs[0] = vstart;
    }
  }
  static void links_at_nverts(Int* nverts, Vertex<3>* vertbuffer, Int vcur,
                              Int np) {
    vertbuffer[*nverts].pnbrs[1-np] = vcur;
    vertbuffer[*nverts].pnbrs[np] = -1;
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
 * \param[in] nplanes
 * The number of planes in the input array.
 *
 */
template <Int dim>
OSH_INLINE void clip(Polytope<dim>* poly, Plane<dim>* planes, Int nplanes) {
  Int* nverts = &poly->nverts;
  if (*nverts <= 0) return;

  // direct access to vertex buffer
  Vertex<dim>* vertbuffer = poly->verts;

  // variable declarations
  Int v, p, np, onv, vcur, vnext, numunclipped;

  // signed distances to the clipping plane
  Real sdists[Polytope<2>::max_verts];
  Real smin, smax;

  // loop over each clip plane
  for (p = 0; p < nplanes; ++p) {
    // calculate signed distances to the clip plane
    onv = *nverts;
    smin = ArithTraits<Real>::max();
    smax = ArithTraits<Real>::min();

    // for marking clipped vertices
    Int clipped[MaxVerts<dim>::value] = {};  // all initialized to zero
    for (v = 0; v < onv; ++v) {
      sdists[v] = planes[p].d + (vertbuffer[v].pos * planes[p].n);
      if (sdists[v] < smin) smin = sdists[v];
      if (sdists[v] > smax) smax = sdists[v];
      if (sdists[v] < 0.0) clipped[v] = 1;
    }

    // skip this face if the poly lies entirely on one side of it
    if (smin >= 0.0) continue;
    if (smax <= 0.0) {
      *nverts = 0;
      return;
    }

    // check all edges and insert new vertices on the bisected edges
    for (vcur = 0; vcur < onv; ++vcur) {
      if (clipped[vcur]) continue;
      for (np = 0; np < dim; ++np) {
        vnext = vertbuffer[vcur].pnbrs[np];
        if (!clipped[vnext]) continue;
        ClipHelper<dim>::links_at_nverts(nverts, vertbuffer, vcur, np);
        vertbuffer[vcur].pnbrs[np] = *nverts;
        vertbuffer[*nverts].pos = wav(vertbuffer[vcur].pos, -sdists[vnext],
                                      vertbuffer[vnext].pos, sdists[vcur]);
        (*nverts)++;
      }
    }

    // for each new vert, search around the poly for its new neighbors
    // and doubly-link everything
    ClipHelper<dim>::relink(onv, nverts, vertbuffer);

    // go through and compress the vertex list, removing clipped verts
    // and re-indexing accordingly (reusing `clipped` to re-index everything)
    numunclipped = 0;
    for (v = 0; v < *nverts; ++v) {
      if (!clipped[v]) {
        vertbuffer[numunclipped] = vertbuffer[v];
        clipped[v] = numunclipped++;
      }
    }
    *nverts = numunclipped;
    for (v = 0; v < *nverts; ++v)
      for (np = 0; np < dim; ++np)
        vertbuffer[v].pnbrs[np] = clipped[vertbuffer[v].pnbrs[np]];
  }
}

template <Int dim, Int nplanes>
OSH_INLINE Polytope<dim> clip(Polytope<dim> poly, Few<Plane<dim>, nplanes> planes) {
  clip(&poly, planes.data(), nplanes);
  return poly;
}

/**
 * \brief Initialize a polyhedron as a tetrahedron.
 *
 * \returns
 * The initialized polyhedron.
 *
 * \param [in] verts
 * An array of four vectors, giving the vertices of the tetrahedron.
 *
 */
OSH_INLINE Polytope<3> init_tet(Few<Vector<3>, 4> verts) {
  Polytope<3> poly;
  // direct access to vertex buffer
  Vertex<3>* vertbuffer = poly.verts;
  Int* nverts = &poly.nverts;

  // initialize graph connectivity
  *nverts = 4;
  vertbuffer[0].pnbrs[0] = 1;
  vertbuffer[0].pnbrs[1] = 3;
  vertbuffer[0].pnbrs[2] = 2;
  vertbuffer[1].pnbrs[0] = 2;
  vertbuffer[1].pnbrs[1] = 3;
  vertbuffer[1].pnbrs[2] = 0;
  vertbuffer[2].pnbrs[0] = 0;
  vertbuffer[2].pnbrs[1] = 3;
  vertbuffer[2].pnbrs[2] = 1;
  vertbuffer[3].pnbrs[0] = 1;
  vertbuffer[3].pnbrs[1] = 2;
  vertbuffer[3].pnbrs[2] = 0;

  // copy vertex coordinates
  for (Int v = 0; v < 4; ++v) vertbuffer[v].pos = verts[v];

  return poly;
}

OSH_INLINE void norm(Vector<3>& v) { v = osh::normalize(v); }

OSH_INLINE Plane<3> tet_face_from_verts(Vector<3> a, Vector<3> b, Vector<3> c) {
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
OSH_INLINE Few<Plane<3>, 4> tet_faces_from_verts(Few<Vector<3>, 4> verts) {
  Few<Plane<3>, 4> faces;
  faces[0] = tet_face_from_verts(verts[3], verts[2], verts[1]);
  faces[1] = tet_face_from_verts(verts[0], verts[2], verts[3]);
  faces[2] = tet_face_from_verts(verts[0], verts[3], verts[1]);
  faces[3] = tet_face_from_verts(verts[0], verts[1], verts[2]);
  return faces;
}

constexpr Int num_moments_3d(Int order) {
  return ((order + 1) * (order + 2) * (order + 3) / 6);
}

constexpr Int num_moments_2d(Int order) {
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
 *
 */

template <Int polyorder>
OSH_INLINE void reduce(Polytope<3>* poly, Real* moments) {
  Int* nverts = &poly->nverts;
  if (*nverts <= 0) return;

  // var declarations
  Real sixv;
  Int np, m, i, j, k, corder;
  Int vstart, pstart, vcur, vnext, pnext;
  Vector<3> v0, v1, v2;

  // direct access to vertex buffer
  Vertex<3>* vertbuffer = poly->verts;

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
  for (vstart = 0; vstart < *nverts; ++vstart)
    for (pstart = 0; pstart < 3; ++pstart) {
      // skip this face if we have marked it
      if (emarks[vstart][pstart]) continue;

      // initialize face looping
      pnext = pstart;
      vcur = vstart;
      emarks[vcur][pnext] = 1;
      vnext = vertbuffer[vcur].pnbrs[pnext];
      v0 = vertbuffer[vcur].pos;

      // move to the second edge
      for (np = 0; np < 3; ++np)
        if (vertbuffer[vnext].pnbrs[np] == vcur) break;
      vcur = vnext;
      pnext = (np + 1) % 3;
      emarks[vcur][pnext] = 1;
      vnext = vertbuffer[vcur].pnbrs[pnext];

      // make a triangle fan using edges
      // and first vertex
      while (vnext != vstart) {
        v2 = vertbuffer[vcur].pos;
        v1 = vertbuffer[vnext].pos;

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
          if (vertbuffer[vnext].pnbrs[np] == vcur) break;
        vcur = vnext;
        pnext = (np + 1) % 3;
        emarks[vcur][pnext] = 1;
        vnext = vertbuffer[vcur].pnbrs[pnext];
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
 * 0 for constant (1 moment), 1 for linear (4 moments),
 * 2 for quadratic (10 moments), etc.
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
OSH_INLINE void reduce(Polytope<2>* poly, Real* moments) {
  Int* nverts = &poly->nverts;
  if (*nverts <= 0) return;

  // var declarations
  Int vcur, vnext, m, i, j, corder;
  Real twoa;
  Vector<2> v0, v1;

  // direct access to vertex buffer
  Vertex<2>* vertbuffer = poly->verts;

  // zero the moments
  for (m = 0; m < num_moments_2d(polyorder); ++m) moments[m] = 0.0;

  // Storage for coefficients
  // keep two layers of the triangle of coefficients
  Int prevlayer = 0;
  Int curlayer = 1;
  Real D[polyorder + 1][2];
  Real C[polyorder + 1][2];

  // iterate over edges and compute a sum over simplices
  for (vcur = 0; vcur < *nverts; ++vcur) {
    vnext = vertbuffer[vcur].pnbrs[0];
    v0 = vertbuffer[vcur].pos;
    v1 = vertbuffer[vnext].pos;
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
OSH_INLINE Real integrate(Polytope<dim> polytope,
                          Polynomial<dim, order> polynomial) {
  Real moments[decltype(polynomial)::nterms] = {};
  reduce<order>(&polytope, moments);
  Real result = 0;
  for (Int i = 0; i < decltype(polynomial)::nterms; ++i)
    result += moments[i] * polynomial.coeffs[i];
  return result;
}

template <Int dim>
OSH_INLINE Real measure(Polytope<dim> polytope) {
  return integrate(polytope, Polynomial<dim, 0>{{1}});
}

}  // end namespace r3d

}  // end namespace osh

#endif
