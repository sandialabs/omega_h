#include "Omega_h_transfer_face.hpp"

#include "Omega_h_affine.hpp"
#include "Omega_h_conserve.hpp"
#include "Omega_h_fit.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_metric.hpp"
#include "Omega_h_quality.hpp"
#include "Omega_h_shape.hpp"

namespace Omega_h {
  void transer_div_free_face_flux(Mesh *source_mesh, 
				  Mesh *target_mesh,
				  Int key_dim,
				  LOs source_keys,
				  LOs keys2prods,
				  LOs prods2target_elements,
				  Read<Real> sourceFluxes,
				  Write<Real> targetFluxes )
  {
    constexpr Int spaceDim  = 3;
    constexpr Int nodesPerElement = 4;
    constexpr Int facesPerElement = 4;
    constexpr Int momentOrder = 2;
    constexpr Int maxElementsPerCavity = 50;
    constexpr Int maxFacePerCavity = maxElementsPerCavity*facesPerElement;
    
    auto const keys2source_elements = source_mesh->ask_up(key_dim,spaceDim);

    auto const targetMeshElements_to_faces = target_mesh->ask_down(spaceDim,spaceDim-1);
    auto const sourceMeshElements_to_faces = source_mesh->ask_down(spaceDim,spaceDim-1);

    auto const targetMeshElements_to_nodes = target_mesh->ask_verts_of(spaceDim);
    auto const sourceMeshElements_to_nodes = source_mesh->ask_verts_of(spaceDim);

    auto const targetCoordinates = target_mesh->coords();
    auto const sourceCoordinates = source_mesh->coords();

    auto functor = OMEGA_H_LAMBDA( LO const cavity ) {
      
      LO targetElements_to_MeshElements[maxElementsPerCavity];
      LO targetBegin = keys2prods[cavity];
      LO targetEnd = keys2prods[cavity+1];
      const Int numTargetElements = static_cast<int>(targetEnd - targetBegin);
      for (Int i=0; i<numTargetElements; ++i)
	targetElements_to_MeshElements[i] = prods2target_elements[targetBegin+i];

      LO sourceElements_to_MeshElements[maxElementsPerCavity];
      LO key = source_keys[cavity];
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
	  sourceElementFace_to_MeshFace[e][f] = soureMeshElements_to_faces.ab2b[meshElement*facesPerElement + f];
	}
      }

      Int targetElementFace_to_targetFace[maxElementsPerCavity][facesPerElement];
      bool targetElementFaceIsSurface[maxElementsPerCavity][facesPerElement];
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

      Int targetElementFaceOrientations[maxElementsPerCavity][facesPerElement];
      for (int elem = 0; elem<numTargetElements; ++elem) {
	LO meshElement =  targetElements_to_MeshElements[elem];
	for (int face = 0; face<facesPerElement; ++face) {
	  auto const code = targetMeshElements_to_faces.codes[meshElement*facesPerElement + face];
	  targetElementFaceOrientations[elem][face] = Omega_h::code_is_flipped(code);
	}
      }

      Int sourceElementFaceOrientations[maxElementsPerCavity][facesPerElement];
      for (int elem = 0; elem<numSourceElements; ++elem) {
	LO meshElement =  sourceElements_to_MeshElements[elem];
	for (int face = 0; face<facesPerElement; ++face) {
	  auto const code = sourceMeshElements_to_faces.codes[meshElement*facesPerElement + face];
	  sourceElementFaceOrientations[elem][face] = Omega_h::code_is_flipped(code);
	}
      }
      
      for (int elem = 0; elem<numTargetElements; ++elem) {
	LO meshElement =  targetElements_to_MeshElements[elem];
	auto const elementNodes = Omega_h::gather_verts<nodesPerElement>(targetMeshElements_to_nodes,meshElement);
	auto const elementCoordinates = Omega_h::gather_vectors<nodesPerElement,spaceDim>(targetCoordinates,elementNodes);
	//elementCoordinates( dimension(0:2), local_node_number(0:3) );
      }
      for (int elem = 0; elem<numSourceElements; ++elem) {
	LO meshElement =  sourceElements_to_MeshElements[elem];
	auto const elementNodes = Omega_h::gather_verts<nodesPerElement>(sourceMeshElements_to_nodes,meshElement);
	auto const elementCoordinates = Omega_h::gather_vectors<nodesPerElement,spaceDim>(sourceCoordinates,elementNodes);
	//elementCoordinates( dimension(0:2), local_node_number(0:3) );
      }

      for (int elem = 0; elem<numSourceElements; ++elem) {
	for (int face = 0; face<facesPerElement; ++face) {
	  LO meshFace = sourceElementFace_to_MeshFace[elem][face];
	  const double faceFlux = sourceFluxes[meshFace];
	}
      }
      for (int elem = 0; elem<numTargetElements; ++elem) {
	for (int face = 0; face<facesPerElement; ++face) {
	  LO meshFace = targetElementFace_to_MeshFace[elem][face];
	  double &faceFlux = targetFluxes[meshFace];
	}

    };//end functor

    const LO nCavities = source_keys.size(); 
    parallel_for( nCavities, std::move(functor) );
    return;
  }

}//end namespace Omega_h

