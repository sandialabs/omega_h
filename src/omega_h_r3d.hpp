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
 * This program was prepared by Los Alamos National Security, LLC at Los Alamos National
 * Laboratory (LANL) under contract No. DE-AC52-06NA25396 with the U.S. Department of Energy (DOE).
 * All rights in the program are reserved by the DOE and Los Alamos National Security, LLC.
 * Permission is granted to the public to copy and use this software without charge, provided that
 * this Notice and any statement of authorship are reproduced on all copies.  Neither the U.S.
 * Government nor LANS makes any warranty, express or implied, or assumes any liability
 * or responsibility for the use of this software.
 * ----
 */

#include "omega_h_math.hpp"

namespace osh {

namespace r3d {

template <Int dim>
struct Plane {
  Vector<dim> n; /*!< Unit-length normal vector. */
  Real d; /*!< Signed perpendicular distance to the origin. */
};

template <Int dim>
struct Vertex {
  Int pnbrs[dim]; /*!< Neighbor indices. */
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
struct Poly {
  enum { max_verts = MaxVerts<dim>::value };
  Vertex<dim> verts[max_verts]; /*!< Vertex buffer. */
  Int nverts; /*!< Number of vertices in the buffer. */
};

OSH_INLINE Vector<3> wav(Vector<3> va, Real wa, Vector<3> vb, Real wb) {
  Vector<3> vr;
  for (Int i = 0; i < 3; ++i)
    vr[i] = (wa * va[i] + wb * vb[i]) / (wa + wb);
  return vr;
}

/**
 * \brief Clip a polyhedron against an arbitrary number of clip planes.
 *
 * \param [in, out] poly
 * The polyhedron to be clipped.
 *
 * \param [in] planes
 * An array of planes against which to clip this polyhedron.
 *
 * \param[in] nplanes
 * The number of planes in the input array.
 *
 */
void clip(Poly<3>* poly, Plane<3> planes[], Int nplanes) {

  // direct access to vertex buffer
  Vertex<3>* vertbuffer = poly->verts;
  Int* nverts = &poly->nverts;
  if(*nverts <= 0) return;

  // variable declarations
  Int v, p, np, onv, vcur, vnext, vstart,
      pnext, numunclipped;

  // signed distances to the clipping plane
  Real sdists[Poly<3>::max_verts];
  Real smin, smax;

  // loop over each clip plane
  for (p = 0; p < nplanes; ++p) {

    // calculate signed distances to the clip plane
    onv = *nverts;
    smin = ArithTraits<Real>::max();
    smax = ArithTraits<Real>::min();

    // for marking clipped vertices
    Int clipped[Poly<3>::max_verts] = {}; // all initialized to zero
    for (v = 0; v < onv; ++v) {
      sdists[v] = planes[p].d + (vertbuffer[v].pos * planes[p].n);
      if(sdists[v] < smin) smin = sdists[v];
      if(sdists[v] > smax) smax = sdists[v];
      if(sdists[v] < 0.0) clipped[v] = 1;
    }

    // skip this face if the poly lies entirely on one side of it
    if(smin >= 0.0) continue;
    if(smax <= 0.0) {
      *nverts = 0;
      return;
    }

    // check all edges and insert new vertices on the bisected edges
    for(vcur = 0; vcur < onv; ++vcur) {
      if(clipped[vcur]) continue;
      for(np = 0; np < 3; ++np) {
        vnext = vertbuffer[vcur].pnbrs[np];
        if(!clipped[vnext]) continue;
        vertbuffer[*nverts].pnbrs[0] = vcur;
        vertbuffer[vcur].pnbrs[np] = *nverts;
        vertbuffer[*nverts].pos = wav(
            vertbuffer[vcur].pos, -sdists[vnext],
            vertbuffer[vnext].pos, sdists[vcur]);
        (*nverts)++;
      }
    }

    // for each new vert, search around the faces for its new neighbors
    // and doubly-link everything
    for(vstart = onv; vstart < *nverts; ++vstart) {
      vcur = vstart;
      vnext = vertbuffer[vcur].pnbrs[0];
      do {
        for(np = 0; np < 3; ++np) if(vertbuffer[vnext].pnbrs[np] == vcur) break;
        vcur = vnext;
        pnext = (np+1)%3;
        vnext = vertbuffer[vcur].pnbrs[pnext];
      } while(vcur < onv);
      vertbuffer[vstart].pnbrs[2] = vcur;
      vertbuffer[vcur].pnbrs[1] = vstart;
    }

    // go through and compress the vertex list, removing clipped verts
    // and re-indexing accordingly (reusing `clipped` to re-index everything)
    numunclipped = 0;
    for(v = 0; v < *nverts; ++v) {
      if(!clipped[v]) {
        vertbuffer[numunclipped] = vertbuffer[v];
        clipped[v] = numunclipped++;
      }
    }
    *nverts = numunclipped;
    for(v = 0; v < *nverts; ++v)
      for(np = 0; np < 3; ++np)
        vertbuffer[v].pnbrs[np] = clipped[vertbuffer[v].pnbrs[np]];
  }
}

}  // end namespace r3d

}  // end namespace osh

#endif
