#include "Omega_h_transfer_face.hpp"

#include "Omega_h_align.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_qr.hpp"

#include <r3d.hpp>


namespace Omega_h {

using Scalar = double;

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
  const double Jbox = 6.0*tetVolume;

  const Vector x0 = nodalCoordinates[0];
  const Vector g1 = nodalCoordinates[1] - x0; 
  const Vector g2 = nodalCoordinates[2] - x0; 
  const Vector g3 = nodalCoordinates[3] - x0; 

  const double a1 = -2.0*faceFlux[3]/Jbox;
  const double a2 = -2.0*faceFlux[1]/Jbox;
  const double a3 = -2.0*faceFlux[0]/Jbox;
  const double a4 = +2.0*( faceFlux[0] + faceFlux[1] + faceFlux[2] + faceFlux[3] )/Jbox;

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

    auto const targetMeshElements_to_nodes = target_mesh->ask_verts_of(spaceDim);
    auto const sourceMeshElements_to_nodes = source_mesh->ask_verts_of(spaceDim);

    auto const targetCoordinates = target_mesh->coords();
    auto const sourceCoordinates = source_mesh->coords();

    const LOs src_keys = source_keys;
    const LOs key2prod = keys2prods;
    const LOs prd2elem = prods2target_elements;
    const Read<Real> srcflux = sourceFluxes;

    auto functor = OMEGA_H_LAMBDA( LO const cavity ) {
      
      LO targetElements_to_MeshElements[maxElementsPerCavity];
      LO targetBegin = key2prod[cavity];
      LO targetEnd = key2prod[cavity+1];
      const Int numTargetElements = static_cast<int>(targetEnd - targetBegin);
      for (Int i=0; i<numTargetElements; ++i)
	targetElements_to_MeshElements[i] = prd2elem[targetBegin+i];

      LO sourceElements_to_MeshElements[maxElementsPerCavity];
      LO key = src_keys[cavity];
      LO sourceBegin = keys2source_elements.a2ab[key];
      LO sourceEnd = keys2source_elements.a2ab[key+1];
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
      bool targetElementFaceIsSurface    [maxElementsPerCavity][facesPerElement];
      for (int elem = 0; elem<numTargetElements; ++elem) {
	for (int face = 0; face<facesPerElement; ++face) {
	  targetElementFace_to_targetFace[elem][face] = -1;
	  targetElementFaceIsSurface[elem][face] = false;
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
	    if (n==1) targetElementFaceIsSurface[elem][face] = true;
	  }
	}
      }

      int number_of_surface_faces = 0;
      for (int elem_trg = 0; elem_trg < numTargetElements; ++elem_trg)  
        for (int face_trg = 0; face_trg < facesPerElement; ++face_trg)  
          if (targetElementFaceIsSurface[elem_trg][face_trg]) ++number_of_surface_faces;
      const int nunknown = number_of_target_faces - number_of_surface_faces;

      const int nface_trg = number_of_target_faces;
      const int nelem_trg = numTargetElements;
      const int nQsize = nunknown < nelem_trg ? nunknown : nelem_trg;
      int nsize = nface_trg+nQsize;

      if (4*nelem_trg != 2*nunknown + number_of_surface_faces) 
        printf("Internal check of elements/faces and unknowns failed.");


      Int targetElementFaceOrientations[maxElementsPerCavity][facesPerElement];
      for (int elem = 0; elem<numTargetElements; ++elem) {
	LO meshElement =  targetElements_to_MeshElements[elem];
	for (int face = 0; face<facesPerElement; ++face) {
	  auto const code = targetMeshElements_to_faces.codes[meshElement*facesPerElement + face];
          const bool flipped = Omega_h::code_is_flipped(code);
          const int sign = flipped ? -1 : +1;
	  targetElementFaceOrientations[elem][face] = sign;
	}
      }

      Int sourceElementFaceOrientations[maxElementsPerCavity][facesPerElement];
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
          const size_t iface = targetElementFace_to_targetFace[elem][face];
          Q(elem,iface) = sign;
        }
      }

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
            const size_t i = targetElementFace_to_targetFace[elem][iface];
            const size_t j = targetElementFace_to_targetFace[elem][jface];
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
	      const double faceFlux = srcflux[meshFace];
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

            const size_t i = targetElementFace_to_targetFace[elem_trg][iface];
            f(i) += integral;
          }
        }
      }


      constexpr int maxSize = maxFacePerCavity+maxElementsPerCavity;
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
            const int sign = targetElementFaceOrientations[elem_trg][face_trg];
            const int i  = targetElementFace_to_targetFace[elem_trg][face_trg];
            A(nsize,i) = sign;
            A(i,nsize) = sign;
            b(nsize)   = targetFluxes[targetElementFace_to_MeshFace[elem_trg][face_trg]];
            ++nsize;
          }
        }
      }
      if(maxSize<nsize)printf("Maximum cavity size exceeded. Increase max number of elements.");

      Omega_h::Vector<maxSize> X = Omega_h::solve_using_qr(nsize, nsize, A, b);

      for (int elem_trg = 0; elem_trg < nelem_trg; ++elem_trg) {
        for (int face_trg = 0; face_trg < facesPerElement; ++face_trg) {
          const size_t faceID_trg = targetElementFace_to_targetFace[elem_trg][face_trg];
	  LO meshFace = targetElementFace_to_MeshFace[elem_trg][face_trg];
          targetFluxes[meshFace] = X(faceID_trg);
        }
      }
    };//end functor

    const LO nCavities = source_keys.size(); 
    parallel_for( nCavities, std::move(functor) );
    return;
  }

}//end namespace Omega_h

