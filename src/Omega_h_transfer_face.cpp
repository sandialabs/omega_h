#include <utility>
#include <iostream>

#include "Omega_h_transfer_face.hpp"

#include "Omega_h_align.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_qr.hpp"

#include <r3d.hpp>


namespace Omega_h {

using Scalar = Real;

OMEGA_H_INLINE Scalar
                       tet4Volume(const Scalar *x, const Scalar *y, const Scalar *z) {
  const Scalar v1[] = {x[1] - x[0], y[1] - y[0], z[1] - z[0]};
  const Scalar v2[] = {x[2] - x[0], y[2] - y[0], z[2] - z[0]};
  const Scalar v3[] = {x[3] - x[0], y[3] - y[0], z[3] - z[0]};
  const Scalar x2 = v2[0];
  const Scalar x3 = v3[0];
  const Scalar y2 = v2[1];
  const Scalar y3 = v3[1];
  const Scalar z2 = v2[2];
  const Scalar z3 = v3[2];
  const Scalar volume =
      (1. / 6.) * (v1[0] * (y2 * z3 - y3 * z2) + v1[1] * (z2 * x3 - z3 * x2) +
                   v1[2] * (x2 * y3 - x3 * y2));
  return volume;
}

OMEGA_H_INLINE
void
elementPhysicalFacePolynomial( /*input*/
                              const Omega_h::Few< Omega_h::Vector<3>, 4> &nodalCoordinates,
                              const Omega_h::Few<Scalar,4> &faceFlux,
                              /*output*/
                              Omega_h::Vector<3> &constantTerm,
                              Scalar &linearFactor ) 
{
  typedef Omega_h::Vector<3> Vector;
  Scalar x[4], y[4], z[4];
  for (int node=0; node<4; ++node) {
    const Vector &coord = nodalCoordinates[node];
    x[node] = coord[0];
    y[node] = coord[1];
    z[node] = coord[2];
  }
  const Scalar tetVolume = tet4Volume(x,y,z);
  const Real Jbox = 6.0*tetVolume;

  const Vector x0 = nodalCoordinates[0];
  const Vector g1 = nodalCoordinates[1] - x0; 
  const Vector g2 = nodalCoordinates[2] - x0; 
  const Vector g3 = nodalCoordinates[3] - x0; 

  const Real a1 = -2.0*faceFlux[3]/Jbox;
  const Real a2 = -2.0*faceFlux[1]/Jbox;
  const Real a3 = -2.0*faceFlux[0]/Jbox;
  const Real a4 = +2.0*( faceFlux[0] + faceFlux[1] + faceFlux[2] + faceFlux[3] )/Jbox;

  constantTerm = a1*g1 + a2*g2 + a3*g3 - a4*x0;
  linearFactor = a4; 

  return;
}


  void transer_div_free_face_flux(      Mesh *source_mesh, 
				        Mesh *target_mesh,
				  const Int key_dim,
				  const LOs source_keys,
				  const LOs keys2prods,
				  const LOs prods2target_elements,
				  const Read<Real> sourceFluxes,
				        Write<Real> targetFluxes )
  {
    constexpr Int spaceDim  = 3;
    constexpr Int nodesPerElement = 4;
    constexpr Int nodesPerFace    = 3;
    constexpr Int facesPerElement = 4;
    constexpr Int momentOrder = 2;
    constexpr Int maxElementsPerCavity = 50;
    constexpr Int maxFacePerCavity = maxElementsPerCavity*facesPerElement;
    
    typedef r3d::Few<r3d::Vector<spaceDim>,nodesPerElement> Tet4;
    typedef r3d::Polytope<spaceDim>                         Polytope;
    typedef r3d::Polynomial<spaceDim,momentOrder>           Polynomial;

    auto const keys2source_elements = source_mesh->ask_up(key_dim,spaceDim);

    auto const targetMeshElements_to_faces = target_mesh->ask_down(spaceDim,spaceDim-1);
    auto const sourceMeshElements_to_faces = source_mesh->ask_down(spaceDim,spaceDim-1);

    auto const targetMeshFaces_to_elements = target_mesh->ask_up(spaceDim-1,spaceDim);
    auto const sourceMeshFaces_to_elements = source_mesh->ask_up(spaceDim-1,spaceDim);

    auto const targetMeshElements_to_nodes = target_mesh->ask_verts_of(spaceDim);
    auto const sourceMeshElements_to_nodes = source_mesh->ask_verts_of(spaceDim);

    auto const targetMeshFaces_to_nodes = target_mesh->ask_verts_of(spaceDim-1);
    auto const sourceMeshFaces_to_nodes = source_mesh->ask_verts_of(spaceDim-1);

    auto const targetCoordinates = target_mesh->coords();
    auto const sourceCoordinates = source_mesh->coords();

    const LOs src_keys = source_keys;
    const LOs key2prod = keys2prods;
    const LOs prd2elem = prods2target_elements;
    const Read<Real> srcflux = sourceFluxes;
    const Write<Real> trgflux = targetFluxes;

    auto functor = OMEGA_H_LAMBDA( LO const cavity ) {
      
      LO targetElements_to_MeshElements[maxElementsPerCavity];
      const LO targetBegin = key2prod[cavity];
      const LO targetEnd = key2prod[cavity+1];
      const Int numTargetElements = static_cast<int>(targetEnd - targetBegin);
      for (Int i=0; i<numTargetElements; ++i)
	targetElements_to_MeshElements[i] = prd2elem[targetBegin+i];

      LO sourceElements_to_MeshElements[maxElementsPerCavity];
      const LO key = src_keys[cavity];
      const LO sourceBegin = keys2source_elements.a2ab[key];
      const LO sourceEnd = keys2source_elements.a2ab[key+1];
      const Int numSourceElements = static_cast<int>(sourceEnd - sourceBegin);
      for (Int i=0; i<numSourceElements; ++i)
	sourceElements_to_MeshElements[i] = keys2source_elements.ab2b[sourceBegin+i];

      LO targetElementFace_to_MeshFace[maxElementsPerCavity][facesPerElement];
      for (Int e=0; e<numTargetElements; ++e) {
	LO meshElement =  targetElements_to_MeshElements[e];
	for (Int f=0; f<facesPerElement; ++f) {
	  targetElementFace_to_MeshFace[e][f] = targetMeshElements_to_faces.ab2b[meshElement*facesPerElement + f];
	}
      }

      LO sourceElementFace_to_MeshFace[maxElementsPerCavity][facesPerElement];
      for (Int e=0; e<numSourceElements; ++e) {
	LO meshElement =  sourceElements_to_MeshElements[e];
	for (Int f=0; f<facesPerElement; ++f) {
	  sourceElementFace_to_MeshFace[e][f] = sourceMeshElements_to_faces.ab2b[meshElement*facesPerElement + f];
	}
      }

      Int targetElementFace_to_targetFace[maxElementsPerCavity][facesPerElement];
      Int sourceElementFace_to_sourceFace[maxElementsPerCavity][facesPerElement];
      bool targetElementFaceIsSurface    [maxElementsPerCavity][facesPerElement];
      bool sourceElementFaceIsSurface    [maxElementsPerCavity][facesPerElement];
      bool targetElementFaceIsMeshSurface[maxElementsPerCavity][facesPerElement];
      bool sourceElementFaceIsMeshSurface[maxElementsPerCavity][facesPerElement];
      for (int elem = 0; elem<maxElementsPerCavity; ++elem) {
	for (int face = 0; face<facesPerElement; ++face) {
	  targetElementFace_to_targetFace[elem][face] = -1;
	  targetElementFaceIsSurface     [elem][face] = false;
	  targetElementFaceIsMeshSurface [elem][face] = false;
	}
      }
      for (int elem = 0; elem<maxElementsPerCavity; ++elem) {
	for (int face = 0; face<facesPerElement; ++face) {
	  sourceElementFace_to_sourceFace[elem][face] = -1;
	  sourceElementFaceIsSurface     [elem][face] = false;
	  sourceElementFaceIsMeshSurface [elem][face] = false;
	}
      }
      Int number_of_target_faces = 0;
      for (int elem = 0; elem<numTargetElements; ++elem) {
	for (int face = 0; face<facesPerElement; ++face) {
	  if ( targetElementFace_to_targetFace[elem][face] == -1) {
	    int n = 0;
	    for (int ielem = elem; ielem < numTargetElements && n < 2; ++ielem) {
	      for (int iface = 0; iface < facesPerElement && n < 2; ++iface) {
		if ( targetElementFace_to_MeshFace[elem][face] == targetElementFace_to_MeshFace[ielem][iface] ) {
		  if ( targetElementFace_to_targetFace[ielem][iface] == -1) {
		    targetElementFace_to_targetFace[ielem][iface] = number_of_target_faces; 
		    ++n;
		  }
		}
	      }
	    }
	    ++number_of_target_faces;
	    if (n==1) {
              const int f = targetElementFace_to_MeshFace[elem][face];
              LO Begin = targetMeshFaces_to_elements.a2ab[f];
              LO End   = targetMeshFaces_to_elements.a2ab[f+1];
              const Int numElements = static_cast<int>(End - Begin);
              if (1<numElements) targetElementFaceIsSurface    [elem][face] = true;
              else               targetElementFaceIsMeshSurface[elem][face] = true;
            }
	  }
	}
      }
      Int number_of_source_faces = 0;
      for (int elem = 0; elem<numSourceElements; ++elem) {
	for (int face = 0; face<facesPerElement; ++face) {
	  if ( sourceElementFace_to_sourceFace[elem][face] == -1) {
	    int n = 0;
	    for (int ielem = elem; ielem < numSourceElements && n < 2; ++ielem) {
	      for (int iface = 0; iface < facesPerElement && n < 2; ++iface) {
		if ( sourceElementFace_to_MeshFace[elem][face] == sourceElementFace_to_MeshFace[ielem][iface] ) {
		  if ( sourceElementFace_to_sourceFace[ielem][iface] == -1) {
		    sourceElementFace_to_sourceFace[ielem][iface] = number_of_source_faces; 
		    ++n;
		  }
		}
	      }
	    }
	    ++number_of_source_faces;
	    if (n==1) {
              const int f = sourceElementFace_to_MeshFace[elem][face];
              LO Begin = sourceMeshFaces_to_elements.a2ab[f];
              LO End   = sourceMeshFaces_to_elements.a2ab[f+1];
              const Int numElements = static_cast<int>(End - Begin);
              if (1<numElements) sourceElementFaceIsSurface    [elem][face] = true;
              else               sourceElementFaceIsMeshSurface[elem][face] = true;
            }
	  }
	}
      }

      int number_of_target_surface_faces = 0;
      for (int elem_trg = 0; elem_trg < numTargetElements; ++elem_trg)  
        for (int face_trg = 0; face_trg < facesPerElement; ++face_trg)  
          if (targetElementFaceIsSurface[elem_trg][face_trg]) ++number_of_target_surface_faces;

      int number_of_source_surface_faces = 0;
      for (int elem_trg = 0; elem_trg < numSourceElements; ++elem_trg)  
        for (int face_trg = 0; face_trg < facesPerElement; ++face_trg)  
          if (sourceElementFaceIsSurface[elem_trg][face_trg]) ++number_of_source_surface_faces;

      int number_of_target_mesh_surface_faces = 0;
      for (int elem_trg = 0; elem_trg < numTargetElements; ++elem_trg)  
        for (int face_trg = 0; face_trg < facesPerElement; ++face_trg)  
          if (targetElementFaceIsMeshSurface[elem_trg][face_trg]) ++number_of_target_mesh_surface_faces;

      int number_of_source_mesh_surface_faces = 0;
      for (int elem_trg = 0; elem_trg < numSourceElements; ++elem_trg)  
        for (int face_trg = 0; face_trg < facesPerElement; ++face_trg)  
          if (sourceElementFaceIsMeshSurface[elem_trg][face_trg]) ++number_of_source_mesh_surface_faces;

      if (number_of_source_surface_faces != number_of_target_surface_faces)  {
        printf(" %s %d **** ERROR ***** Internal check of surface faces failed.\n",__FILE__,__LINE__);
        for (int elem = 0; elem < numTargetElements; ++elem)  
          for (int face = 0; face < facesPerElement; ++face)  
            if (targetElementFaceIsSurface[elem][face]) {
              const int f = targetElementFace_to_MeshFace[elem][face];
              printf ("Surface Target face number:%d\n",f);
              auto const faceNodes       = Omega_h::gather_verts  <nodesPerFace>(targetMeshFaces_to_nodes,f);
              auto const faceCoordinates = Omega_h::gather_vectors<nodesPerFace,spaceDim>(targetCoordinates,faceNodes);
              printf ("Surface Target face number:%d, nodes:%d %d %d\n",f,faceNodes[0],faceNodes[1],faceNodes[2]);
              for (int n = 0; n < nodesPerFace; ++n) 
                printf ("Surface Target node number:%d, coord:%lf %lf %lf\n",faceNodes[n],faceCoordinates(0,n),faceCoordinates(1,n),faceCoordinates(2,n));

              LO Begin = targetMeshFaces_to_elements.a2ab[f];
              LO End = targetMeshFaces_to_elements.a2ab[f+1];
              const Int numElements = static_cast<int>(End - Begin);
              for (Int i=0; i<numElements; ++i)
	         printf(" Target Elements:%d\n",targetMeshFaces_to_elements.ab2b[Begin+i]);
            }

        for (int elem = 0; elem < numSourceElements; ++elem)  
          for (int face = 0; face < facesPerElement; ++face)  
            if (sourceElementFaceIsSurface[elem][face]) {
              const int f = sourceElementFace_to_MeshFace[elem][face];
              printf ("Surface Source face number:%d\n",f);
              auto const faceNodes       = Omega_h::gather_verts  <nodesPerFace>(sourceMeshFaces_to_nodes,f);
              auto const faceCoordinates = Omega_h::gather_vectors<nodesPerFace,spaceDim>(sourceCoordinates,faceNodes);
              printf ("Surface Source face number:%d, nodes:%d %d %d\n",f,faceNodes[0],faceNodes[1],faceNodes[2]);
              for (int n = 0; n < nodesPerFace; ++n) 
                printf ("Surface Source node number:%d, coord:%lf %lf %lf\n",faceNodes[n],faceCoordinates(0,n),faceCoordinates(1,n),faceCoordinates(2,n));

              LO Begin = sourceMeshFaces_to_elements.a2ab[f];
              LO End = sourceMeshFaces_to_elements.a2ab[f+1];
              const Int numElements = static_cast<int>(End - Begin);
              for (Int i=0; i<numElements; ++i)
	         printf(" Source Elements:%d\n",sourceMeshFaces_to_elements.ab2b[Begin+i]);
            }
      }

 
      Omega_h::Few<Omega_h::Few<Omega_h::Vector<spaceDim-1>,nodesPerFace>,maxFacePerCavity> targetSurfaceNodalCoordinates;
      Omega_h::Few<Omega_h::Few<Omega_h::Vector<spaceDim-1>,nodesPerFace>,maxFacePerCavity> sourceSurfaceNodalCoordinates;
      Omega_h::Few<int,maxFacePerCavity> targetSurfaceNormalDirection;
      Omega_h::Few<int,maxFacePerCavity> sourceSurfaceNormalDirection;
      Omega_h::Few<Real,maxFacePerCavity> targetSurfaceNormalPlane;
      Omega_h::Few<Real,maxFacePerCavity> sourceSurfaceNormalPlane;

      int targetElementFaceIsMeshSurfaceIndex[maxElementsPerCavity][facesPerElement];
      Omega_h::Few<int,2> sourceIndexToElementFace[maxFacePerCavity];
      for (int e = 0; e < maxElementsPerCavity; ++e)
        for (int f = 0; f < facesPerElement; ++f)
          targetElementFaceIsMeshSurfaceIndex[e][f] = -1;
      
      for (int f = 0; f < maxFacePerCavity; ++f)
        sourceIndexToElementFace[f][0] = sourceIndexToElementFace[f][1] = -1;

      int numTargetMeshSurfaceElements = 0;
      int numSourceMeshSurfaceElements = 0;
      for (int e = 0; e < numTargetElements; ++e)  {
        for (int f = 0; f < facesPerElement; ++f)  {
          if (targetElementFaceIsMeshSurface[e][f]) {
            const Real epsilon = .000000001;
            Omega_h::Few<Omega_h::Vector<spaceDim>,nodesPerFace> nodalCoordinates;
            {
              const LO meshFace = targetElementFace_to_MeshFace[e][f];
              auto const faceNodes       = Omega_h::gather_verts  <nodesPerFace>(targetMeshFaces_to_nodes,meshFace);
              auto const faceCoordinates = Omega_h::gather_vectors<nodesPerFace,spaceDim>(targetCoordinates,faceNodes);

              for (int n = 0; n < nodesPerFace; ++n) 
                for (int d = 0; d < spaceDim; ++d) 
                  nodalCoordinates[n][d] = faceCoordinates(d,n);
            }
            int normal = -1;
            for (int d = 0; d < spaceDim && normal < 0; ++d) {
              bool all_same = true;
              for (int n = 1; n < nodesPerFace && all_same; ++n) 
                all_same = all_same && std::abs(nodalCoordinates[0][d]-nodalCoordinates[n][d]) < epsilon;
              if (all_same) normal = d;
            }
            if (normal == -1) 
              printf ("%s %d ***** ERROR ***** Surface elements are not axis aligned."
                      "  Need a better algorithm. Bummer.\n",__FILE__,__LINE__);
             
            Omega_h::Few<Omega_h::Vector<spaceDim-1>,nodesPerFace> nodalCoordinates2D;
            for (int d = 0, c = 0; d < spaceDim; ++d) {
              if (d != normal) {
                for (int n = 0; n < nodesPerFace; ++n) 
                  nodalCoordinates2D[n][c] = nodalCoordinates[n][d];
                ++c;
              }
            }{
               const Omega_h::Vector<spaceDim-1> v0 = {{nodalCoordinates2D[1][0]-nodalCoordinates2D[0][0],
                                                        nodalCoordinates2D[1][1]-nodalCoordinates2D[0][1]}};
               const Omega_h::Vector<spaceDim-1> v1 = {{nodalCoordinates2D[2][0]-nodalCoordinates2D[0][0],
                                                        nodalCoordinates2D[2][1]-nodalCoordinates2D[0][1]}};
               const Real cross = v0[0]*v1[1] - v0[1]*v1[0];
               if (cross < 0) {
                  const Omega_h::Vector<spaceDim-1> t = nodalCoordinates2D[1];
                  nodalCoordinates2D[1] = nodalCoordinates2D[2];
                  nodalCoordinates2D[2] = t;
               }
            } 
            targetSurfaceNormalPlane     [numTargetMeshSurfaceElements] = nodalCoordinates[0][normal];
            targetSurfaceNormalDirection [numTargetMeshSurfaceElements] = normal;
            targetSurfaceNodalCoordinates[numTargetMeshSurfaceElements] = nodalCoordinates2D;
            targetElementFaceIsMeshSurfaceIndex[e][f] = numTargetMeshSurfaceElements;
            ++numTargetMeshSurfaceElements;
          }
        }
      }

      for (int e = 0; e < numSourceElements; ++e)  {
        for (int f = 0; f < facesPerElement; ++f)  {
          if (sourceElementFaceIsMeshSurface[e][f]) {
            const Real epsilon = .0000001;
            Omega_h::Few<Omega_h::Vector<spaceDim>,nodesPerFace> nodalCoordinates;
            {
              const LO meshFace = sourceElementFace_to_MeshFace[e][f];
              auto const faceNodes       = Omega_h::gather_verts  <nodesPerFace>(sourceMeshFaces_to_nodes,meshFace);
              auto const faceCoordinates = Omega_h::gather_vectors<nodesPerFace,spaceDim>(sourceCoordinates,faceNodes);

              for (int n = 0; n < nodesPerFace; ++n) 
                for (int d = 0; d < spaceDim; ++d) 
                  nodalCoordinates[n][d] = faceCoordinates(d,n);
            }
            int normal = -1;
            for (int d = 0; d < spaceDim && normal < 0; ++d) {
              bool all_same = true;
              for (int n = 1; n < nodesPerFace && all_same; ++n) 
                all_same = all_same && std::abs(nodalCoordinates[0][d]-nodalCoordinates[n][d]) < epsilon;
              if (all_same) normal = d;
            }
            if (normal == -1) 
              printf (" %s %d ***** ERROR ***** Surface elements are not axis aligned.  "
                      "Need a better algorithm. Bummer.\n",__FILE__,__LINE__);
             
            Omega_h::Few<Omega_h::Vector<spaceDim-1>,nodesPerFace> nodalCoordinates2D;
            for (int d = 0, c = 0; d < spaceDim; ++d) {
              if (d != normal) {
                for (int n = 0; n < nodesPerFace; ++n) 
                  nodalCoordinates2D[n][c] = nodalCoordinates[n][d];
                ++c;
              }
            } 
            {
               const Omega_h::Vector<spaceDim-1> v0 = {{nodalCoordinates2D[1][0]-nodalCoordinates2D[0][0],
                                                        nodalCoordinates2D[1][1]-nodalCoordinates2D[0][1]}};
               const Omega_h::Vector<spaceDim-1> v1 = {{nodalCoordinates2D[2][0]-nodalCoordinates2D[0][0],
                                                        nodalCoordinates2D[2][1]-nodalCoordinates2D[0][1]}};
               const Real cross = v0[0]*v1[1] - v0[1]*v1[0];
               if (cross < 0) {
                  const Omega_h::Vector<spaceDim-1> t = nodalCoordinates2D[1];
                  nodalCoordinates2D[1] = nodalCoordinates2D[2];
                  nodalCoordinates2D[2] = t;
               }
            }
            sourceSurfaceNormalPlane     [numSourceMeshSurfaceElements] = nodalCoordinates[0][normal];
            sourceSurfaceNormalDirection [numSourceMeshSurfaceElements] = normal;
            sourceSurfaceNodalCoordinates[numSourceMeshSurfaceElements] = nodalCoordinates2D;
            sourceIndexToElementFace[numSourceMeshSurfaceElements][0] = e;
            sourceIndexToElementFace[numSourceMeshSurfaceElements][1] = f;
            ++numSourceMeshSurfaceElements;
          }
        }
      }


      constexpr int maxSize = maxFacePerCavity+maxElementsPerCavity;
      const int nunknown = number_of_target_faces - number_of_target_surface_faces - number_of_target_mesh_surface_faces;
      const int nface_trg = number_of_target_faces;
      const int nelem_trg = numTargetElements;
      int nQsize = nelem_trg;

      if (4*nelem_trg != 2*nunknown + number_of_target_surface_faces + number_of_target_mesh_surface_faces) {
        printf(" %s %d **** ERROR ***** Internal check of elements/faces and unknowns failed. %d != %d + %d + %d\n",
           __FILE__,__LINE__,4*nelem_trg, 2*nunknown,number_of_target_surface_faces,number_of_target_mesh_surface_faces);
        printf("Number of target mesh surface faces:%d\n",number_of_target_mesh_surface_faces);
        printf("Number of source mesh surface faces:%d\n",number_of_source_mesh_surface_faces);
      }


      Int targetElementFaceOrientations[maxElementsPerCavity][facesPerElement] = {0};
      for (int elem = 0; elem<numTargetElements; ++elem) {
	LO meshElement =  targetElements_to_MeshElements[elem];
	for (int face = 0; face<facesPerElement; ++face) {
	  auto const code = targetMeshElements_to_faces.codes[meshElement*facesPerElement + face];
          const bool flipped = Omega_h::code_is_flipped(code);
          const int sign = flipped ? -1 : +1;
	  targetElementFaceOrientations[elem][face] = sign;
	}
      }

      Int sourceElementFaceOrientations[maxElementsPerCavity][facesPerElement] = {0};
      for (int elem = 0; elem<numSourceElements; ++elem) {
	LO meshElement =  sourceElements_to_MeshElements[elem];
	for (int face = 0; face<facesPerElement; ++face) {
	  auto const code = sourceMeshElements_to_faces.codes[meshElement*facesPerElement + face];
          const bool flipped = Omega_h::code_is_flipped(code);
          const int sign = flipped ? -1 : +1;
	  sourceElementFaceOrientations[elem][face] = sign;
	}
      }
      
      for (int elem = 0; elem<numTargetElements; ++elem) {
	LO meshElement =  targetElements_to_MeshElements[elem];
	auto const elementNodes = Omega_h::gather_verts<nodesPerElement>(targetMeshElements_to_nodes,meshElement);
	//auto const elementCoordinates = Omega_h::gather_vectors<nodesPerElement,spaceDim>(targetCoordinates,elementNodes);
	//elementCoordinates( dimension(0:2), local_node_number(0:3) );
      }
      for (int elem = 0; elem<numSourceElements; ++elem) {
	LO meshElement =  sourceElements_to_MeshElements[elem];
	auto const elementNodes = Omega_h::gather_verts<nodesPerElement>(sourceMeshElements_to_nodes,meshElement);
	//auto const elementCoordinates = Omega_h::gather_vectors<nodesPerElement,spaceDim>(sourceCoordinates,elementNodes);
	//elementCoordinates( dimension(0:2), local_node_number(0:3) );
      }

      Omega_h::Matrix<maxFacePerCavity,maxFacePerCavity>     M;
      Omega_h::Vector<maxFacePerCavity>                      f;
      Omega_h::Matrix<maxElementsPerCavity,maxFacePerCavity> Q;
      Omega_h::Vector<maxElementsPerCavity>                  q;

      for (int iface = 0; iface < maxFacePerCavity; ++iface) {
        f(iface) = 0;
        for (int jface = 0; jface < maxFacePerCavity; ++jface) {
          M(iface,jface) = 0;
        }
      }
      for (int elem = 0; elem < maxElementsPerCavity; ++elem) {
        q(elem) = 0;
        for (int face = 0; face < maxFacePerCavity; ++face) {
          Q(elem,face) = 0;
        }
      }

      for (int elem = 0; elem < nQsize; ++elem) {
        q(elem) = 0;
        for (int face = 0; face < facesPerElement; ++face) {
          const Int sign = targetElementFaceOrientations[elem][face];
          const Int iface = targetElementFace_to_targetFace[elem][face];
          Q(elem,iface) = sign;
        }
      }
      {
        Int S[maxSize][maxFacePerCavity];
        for (int i = 0; i < maxSize; ++i)
          for (int j = 0; j < maxFacePerCavity; ++j) 
             S[i][j] = 0;
        unsigned cols=0;
        for (int elem = 0; elem < nQsize; ++elem) 
          for (int face = 0; face < facesPerElement; ++face) {
            const unsigned iface = targetElementFace_to_targetFace[elem][face];
            S[elem][iface] = Q(elem,iface);
            if (cols < iface+1) cols = iface+1;
          }
        unsigned rows = nQsize;
        for (int elem_trg = 0; elem_trg < nelem_trg; ++elem_trg) 
          for (int face_trg = 0; face_trg < facesPerElement; ++face_trg) 
            if (targetElementFaceIsSurface    [elem_trg][face_trg] || 
                targetElementFaceIsMeshSurface[elem_trg][face_trg]) {
              const unsigned iface  = targetElementFace_to_targetFace[elem_trg][face_trg];
              S[rows][iface] = 1;
              ++rows;
              if (cols < iface+1) cols = iface+1;
            }
        unsigned rank=0;
        for (unsigned j=0,k=0; j<cols; ++j) {
          bool p = false;
          for (unsigned i=k; i<rows; ++i) {
            if (S[i][j]) {
              if (p) {
                const int aij = S[i][j]*S[k][j];
                for (unsigned n=j; n<cols; ++n) S[i][n] -= aij*S[k][n];
              } else {
                p = true;
                rank = k+1;
                if (i!=k) {
                  for (unsigned n=j; n<cols; ++n) std::swap(S[k][n], S[i][n]);
                }
              }
            }
          }
          if (p) ++k;
        }
        int reduce = rows-rank;
        nQsize -= reduce;
        if (1!=reduce) 
          printf ("%s:%d ************  ERROR Rank:%d != Rows:%d  ************\n",__FILE__,__LINE__,rank,rows);


// Check
        for (int i = 0; i < maxSize; ++i)
          for (int j = 0; j < maxFacePerCavity; ++j) 
             S[i][j] = 0;
        cols=0;
        for (int elem = 0; elem < nQsize; ++elem) 
          for (int face = 0; face < facesPerElement; ++face) {
            const unsigned iface = targetElementFace_to_targetFace[elem][face];
            S[elem][iface] = Q(elem,iface);
            if (cols < iface+1) cols = iface+1;
          }
        rows = nQsize;
        for (int elem_trg = 0; elem_trg < nelem_trg; ++elem_trg) 
          for (int face_trg = 0; face_trg < facesPerElement; ++face_trg) 
            if (targetElementFaceIsSurface    [elem_trg][face_trg] || 
                targetElementFaceIsMeshSurface[elem_trg][face_trg]) {
              const unsigned iface  = targetElementFace_to_targetFace[elem_trg][face_trg];
              S[rows][iface] = 1;
              ++rows;
              if (cols < iface+1) cols = iface+1;
            }
        rank=0;
        for (unsigned j=0,k=0; j<cols; ++j) {
          bool p = false;
          for (unsigned i=k; i<rows; ++i) {
            if (S[i][j]) {
              if (p) {
                const int aij = S[i][j]*S[k][j];
                for (unsigned n=j; n<cols; ++n) S[i][n] -= aij*S[k][n];
              } else {
                p = true;
                rank = k+1;
                if (i!=k) {
                  for (unsigned n=j; n<cols; ++n) std::swap(S[k][n], S[i][n]);
                }
              }
            }
          }
          if (p) ++k;
        }
        reduce = rows-rank;
        if (reduce) 
          printf ("%s:%d ************  ERROR Rank:%d != Rows:%d  ************\n",__FILE__,__LINE__,rank,rows);
      }

      int nsize = nface_trg+nQsize;

      for (int elem = 0; elem < nelem_trg; ++elem) {
        Scalar moments[Polynomial::nterms];
        for (int i=0; i<Polynomial::nterms; ++i) moments[i]=0;
        Omega_h::Few<Omega_h::Vector<spaceDim>,nodesPerElement> nodalCoordinates;

        {
          LO meshElement =  targetElements_to_MeshElements[elem];
          auto const elementNodes = Omega_h::gather_verts<nodesPerElement>(targetMeshElements_to_nodes,meshElement);
          auto const elementCoordinates = Omega_h::gather_vectors<nodesPerElement,spaceDim>(targetCoordinates,elementNodes);

          Tet4 nodalCoordinates_trg;
          for (int n = 0; n < nodesPerElement; ++n) 
             for (int d = 0; d < spaceDim; ++d) 
               nodalCoordinates_trg[n][d] = elementCoordinates(d,n);

          Polytope poly_trg;
          r3d::init(poly_trg, nodalCoordinates_trg);
          r3d::reduce<momentOrder>(poly_trg, moments);

          for (int n=0; n<nodesPerElement; ++n)
             for (int d=0; d<spaceDim; ++d)
                nodalCoordinates[n][d] = nodalCoordinates_trg[n][d];
        }
        for (int iface = 0; iface < facesPerElement; ++iface) {
  
          Omega_h::Few<Scalar,facesPerElement> ifacesomething;
          for (int i=0; i<facesPerElement; ++i) ifacesomething[i]=0;
          ifacesomething[iface]= targetElementFaceOrientations[elem][iface];
          Omega_h::Vector<3> iconstantTerm;
          Scalar ilinearFactor;
          elementPhysicalFacePolynomial( /*input*/
                                         nodalCoordinates,
                                         ifacesomething,
                                         /*output*/
                                         iconstantTerm,
                                         ilinearFactor );

          for (int jface = 0; jface < facesPerElement; ++jface) {
            Omega_h::Few<Scalar,facesPerElement> jfacesomething;
            for (int i=0; i<facesPerElement; ++i) jfacesomething[i]=0;
            jfacesomething[jface]= targetElementFaceOrientations[elem][jface];
            Omega_h::Vector<3> jconstantTerm;
            Scalar jlinearFactor;
            elementPhysicalFacePolynomial( /*input*/
                                         nodalCoordinates,
                                         jfacesomething,
                                         /*output*/
                                         jconstantTerm,
                                         jlinearFactor );

            Scalar src_poly[Polynomial::nterms];
            src_poly[0] = iconstantTerm*jconstantTerm; //inner_product(iconstantTerm,jconstantTerm);
            src_poly[1] = iconstantTerm[0]*jlinearFactor + jconstantTerm[0]*ilinearFactor;
            src_poly[2] = iconstantTerm[1]*jlinearFactor + jconstantTerm[1]*ilinearFactor;
            src_poly[3] = iconstantTerm[2]*jlinearFactor + jconstantTerm[2]*ilinearFactor;
            src_poly[4] = ilinearFactor*jlinearFactor;
            src_poly[5] = 0;
            src_poly[6] = 0;
            src_poly[7] = ilinearFactor*jlinearFactor;
            src_poly[8] = 0;
            src_poly[9] = ilinearFactor*jlinearFactor;
  
            Scalar integral = 0;
            for (int i=0; i<Polynomial::nterms; ++i) {
              integral += moments[i]*src_poly[i];
            }
            const Int i = targetElementFace_to_targetFace[elem][iface];
            const Int j = targetElementFace_to_targetFace[elem][jface];
            M(i,j) += integral;
          }
        }
      }


      for (int elem_trg = 0; elem_trg < nelem_trg; ++elem_trg) {
        Tet4 nodalCoordinates_trg;
        {
          LO meshElement =  targetElements_to_MeshElements[elem_trg];
          auto const elementNodes = Omega_h::gather_verts<nodesPerElement>(targetMeshElements_to_nodes,meshElement);
          auto const elementCoordinates = Omega_h::gather_vectors<nodesPerElement,spaceDim>(targetCoordinates,elementNodes);
          for (int n = 0; n < nodesPerElement; ++n) 
            for (int d = 0; d < spaceDim; ++d)
              nodalCoordinates_trg[n][d] = elementCoordinates(d,n);
        }

        for (int elem_src = 0; elem_src < numSourceElements; ++elem_src) {
          Tet4 nodalCoordinates_src;
          {
            LO meshElement =  sourceElements_to_MeshElements[elem_src];
            auto const elementNodes = Omega_h::gather_verts<nodesPerElement>(sourceMeshElements_to_nodes,meshElement);
            auto const elementCoordinates = Omega_h::gather_vectors<nodesPerElement,spaceDim>(sourceCoordinates,elementNodes);
            for (int n = 0; n < nodesPerElement; ++n) 
              for (int d = 0; d < spaceDim; ++d)
                nodalCoordinates_src[n][d] = elementCoordinates(d,n);
          }
          Scalar moments[Polynomial::nterms] = {0.}; 
          for (int i=0; i<Polynomial::nterms; ++i) moments[i] = 0;  
          Polytope intersection;
          r3d::intersect_simplices(intersection, nodalCoordinates_src, nodalCoordinates_trg);
          r3d::reduce<momentOrder>(intersection, moments);

          Omega_h::Vector<3> jconstantTerm;
          Scalar jlinearFactor;
          {
            Omega_h::Few<Omega_h::Vector<spaceDim>,nodesPerElement> nodalCoordinates;
            for (int n=0; n<nodesPerElement; ++n)
               for (int d=0; d<spaceDim; ++d)
                  nodalCoordinates[n][d] = nodalCoordinates_src[n][d];

            Omega_h::Vector<facesPerElement> faceFluxSource;
            for (int face = 0; face < facesPerElement; ++face) {
              const int sign = sourceElementFaceOrientations[elem_src][face];
	      LO meshFace = sourceElementFace_to_MeshFace[elem_src][face];
	      const Real faceFlux = srcflux[meshFace];
              faceFluxSource[face] = sign * faceFlux;
            }    
            elementPhysicalFacePolynomial( /*input*/
                                           nodalCoordinates,
                                           faceFluxSource,
                                           /*output*/
                                           jconstantTerm,
                                           jlinearFactor );
          }
          for (int iface = 0; iface < facesPerElement; ++iface) {

            Omega_h::Vector<3> iconstantTerm;
            Scalar ilinearFactor;
            {
              Omega_h::Few<Omega_h::Vector<spaceDim>,nodesPerElement> nodalCoordinates;
              for (int n=0; n<nodesPerElement; ++n)
                 for (int d=0; d<spaceDim; ++d)
                    nodalCoordinates[n][d] = nodalCoordinates_trg[n][d];
              Omega_h::Few<Scalar,facesPerElement> ifacesomething;
              for (int i=0; i<facesPerElement; ++i) ifacesomething[i]=0;
              ifacesomething[iface] = targetElementFaceOrientations[elem_trg][iface];
              elementPhysicalFacePolynomial( /*input*/
                                             nodalCoordinates,
                                             ifacesomething,
                                             /*output*/
                                             iconstantTerm,
                                             ilinearFactor );
            }

            Scalar src_poly[Polynomial::nterms];
            for (int i=0; i<Polynomial::nterms; ++i) src_poly[i]=0;
            src_poly[0] = iconstantTerm*jconstantTerm; //inner_product(iconstantTerm,jconstantTerm);
            src_poly[1] = iconstantTerm[0]*jlinearFactor + jconstantTerm[0]*ilinearFactor;
            src_poly[2] = iconstantTerm[1]*jlinearFactor + jconstantTerm[1]*ilinearFactor;
            src_poly[3] = iconstantTerm[2]*jlinearFactor + jconstantTerm[2]*ilinearFactor;
            src_poly[4] = ilinearFactor*jlinearFactor;
            src_poly[5] = 0;
            src_poly[6] = 0;
            src_poly[7] = ilinearFactor*jlinearFactor;
            src_poly[8] = 0;
            src_poly[9] = ilinearFactor*jlinearFactor;

            Scalar integral = 0;
            for (int i=0; i<Polynomial::nterms; ++i) {
              integral += moments[i]*src_poly[i];
            }

            const Int i = targetElementFace_to_targetFace[elem_trg][iface];
            f(i) += integral;
          }
        }
      }


      Omega_h::Matrix<maxSize,maxSize> A;
      Omega_h::Vector<maxSize>         b;
      for (int i=0; i<maxSize; ++i) {
        b(i) = 0;
        for (int j=0; j<maxSize; ++j) {
          A(i,j) = 0;
        }
      }
      for (int i=0; i<nface_trg; ++i) {
        for (int j=0; j<nface_trg; ++j) {
          A(i,j) = M(i,j);
        }
      }
      for (int i=0; i<nQsize; ++i) {
        for (int j=0; j<nface_trg; ++j) {
          A(nface_trg+i,j) = Q(i,j);
          A(j,nface_trg+i) = Q(i,j);
        }
      }

      for (int i=0; i<nface_trg; ++i) b(i) = f(i);
      for (int i=0; i<nQsize; ++i) b(i+nface_trg) = q(i);

      for (int elem_trg = 0; elem_trg < nelem_trg; ++elem_trg) {
        for (int face_trg = 0; face_trg < facesPerElement; ++face_trg) {
          if (targetElementFaceIsSurface[elem_trg][face_trg]) {
            const Int i  = targetElementFace_to_targetFace[elem_trg][face_trg];
            A(nsize,i) = 1;
            A(i,nsize) = 1;
            b(nsize)   = trgflux[targetElementFace_to_MeshFace[elem_trg][face_trg]];
            ++nsize;
          } else if (targetElementFaceIsMeshSurface[elem_trg][face_trg]) {
            const Real epsilon = .0000001;
            const int target_index = targetElementFaceIsMeshSurfaceIndex[elem_trg][face_trg]; 

            const Omega_h::Few<Omega_h::Vector<spaceDim-1>,nodesPerFace> targetNodalCoordinates2D =
              targetSurfaceNodalCoordinates[target_index];
            const Real targetNormalPlane = targetSurfaceNormalPlane     [target_index];
            const int    targetNormalDir   = targetSurfaceNormalDirection [target_index];

            Scalar integral = 0;
            for (int source_index=0; source_index < numSourceMeshSurfaceElements; ++source_index) {
              const Real sourceNormalPlane = sourceSurfaceNormalPlane     [source_index];
              const int  sourceNormalDir   = sourceSurfaceNormalDirection [source_index];
              const int  elem_src = sourceIndexToElementFace[source_index][0];
              const int  face_src = sourceIndexToElementFace[source_index][1];
              if (targetNormalDir == sourceNormalDir && 
                  std::abs(targetNormalPlane-sourceNormalPlane) < epsilon) {
                 
                constexpr Int spaceDim2D = spaceDim - 1;
                typedef r3d::Few<r3d::Vector<spaceDim2D>,nodesPerFace> Tri3;
                typedef r3d::Polytope<spaceDim2D>                 Polytope2D;
                Real measure;
                {
                  Polytope2D intersection;

                  const Omega_h::Few<Omega_h::Vector<spaceDim2D>,nodesPerFace> sourceNodalCoordinates2D =
                    sourceSurfaceNodalCoordinates[source_index];
                  Tri3 nodalCoordinates_src;
                  for (int n = 0; n < nodesPerFace; ++n) 
                    for (int d = 0; d < spaceDim2D; ++d)
                      nodalCoordinates_src[n][d] = sourceNodalCoordinates2D[n][d];
  
                  Tri3 nodalCoordinates_trg;
                  for (int n = 0; n < nodesPerFace; ++n) 
                    for (int d = 0; d < spaceDim2D; ++d)
                      nodalCoordinates_trg[n][d] = targetNodalCoordinates2D[n][d];
  
                  r3d::init(intersection, nodalCoordinates_src);
                  const Real measure_src = r3d::measure(intersection);

                  r3d::intersect_simplices(intersection, nodalCoordinates_src, nodalCoordinates_trg);
                  measure = r3d::measure(intersection);

                  measure /= measure_src;

                  r3d::init(intersection, nodalCoordinates_trg);
                }

                const int sourceMeshFace = sourceElementFace_to_MeshFace[elem_src][face_src];
                const Real sourceFlux = srcflux[sourceMeshFace];

                integral += measure*sourceFlux; 

                const int sourceSign = sourceElementFaceOrientations[elem_src][face_src];
                const int targetSign = targetElementFaceOrientations[elem_trg][face_trg];
                if (sourceSign != targetSign) 
                  printf (" %s %d ***** ERROR ***** Unexpected source and target mesh flux orientations sign flip.\n",
                          __FILE__,__LINE__);
              }
            }
            const Int i  = targetElementFace_to_targetFace[elem_trg][face_trg];
            A(nsize,i) = 1;
            A(i,nsize) = 1;
            b(nsize)   = integral;
            ++nsize;
          }
        }
      }
      if(maxSize<nsize)printf("%s %d Maximum cavity size exceeded. Increase max number of elements.",__FILE__,__LINE__);

      Omega_h::Vector<maxSize> X = Omega_h::solve_using_qr(nsize, nsize, A, b);

      for (int elem_trg = 0; elem_trg < nelem_trg; ++elem_trg) {
        for (int face_trg = 0; face_trg < facesPerElement; ++face_trg) {
          const auto faceID_trg = targetElementFace_to_targetFace[elem_trg][face_trg];
	  LO meshFace = targetElementFace_to_MeshFace[elem_trg][face_trg];
          trgflux[meshFace] = X(faceID_trg);
        }
      }
    };//end functor

    const LO nCavities = source_keys.size(); 
    parallel_for( nCavities, std::move(functor) );
    return;
  }

}//end namespace Omega_h

