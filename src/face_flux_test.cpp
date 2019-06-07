#include <Omega_h_build.hpp>
#include <Omega_h_coarsen.hpp>
#include <Omega_h_refine.hpp>
#include <Omega_h_compare.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_reduce.hpp>
#include <Omega_h_assoc.hpp>
#include <Omega_h_expr.hpp>

#include <cmath>
#include <vector>
#include <string>
#include <stdexcept>

using namespace Omega_h;




void FunctionInitialCondition(
    const MeshDimSets &facesets,
    const MeshDimSets &elemsets,
    const Reals       &coords,
    const LOs         &face_node_id,
    const LOs         &elem_face_id,
    Write<Real>        field) {

  constexpr int N = 3;
  constexpr int F = 4;
  constexpr int D = 3;
  const std::vector<std::string> 
    face_sets={"x-","y-","z-","x+","y+","z+"};
  const std::vector<std::string>
    element_blocks={"body"};
  const std::string expr_string="vector(2*z, 5*x, y)";
  const std::string name="Face flux test";

  for (const std::string &faceset : face_sets) {
    auto fsIter = facesets.find(faceset);
    if(fsIter == facesets.end())
       fail("Faceset block doesn't exist!");
    const auto     faceLids = (fsIter->second);
    const auto   nset_faces = faceLids.size();
    // Edge Centers
    Write<double> Cx(N*nset_faces);
    Write<double> Cy(N*nset_faces);
    Write<double> Cz(N*nset_faces);
    // Edge Vectors
    Write<double> Vx(N*nset_faces);
    Write<double> Vy(N*nset_faces);
    Write<double> Vz(N*nset_faces);

    const Reals cord = coords;
    const LOs    ids = face_node_id;
    auto prepare = OMEGA_H_LAMBDA(int set_face) {
      const auto face = faceLids[set_face];
      int nodes[N];
      for (int i = 0; i < N; ++i)
        nodes[i] = ids[N*face+i];
      double X[N][3];
      for (int i = 0; i < N; ++i)
        for (int j = 0; j < D; ++j)
           X[i][j] = cord[D*nodes[i]+j];

      const auto offset = set_face*N;
      for (int i = 0; i < N; ++i) {
        const int j = (i+1)%N;
        Vx[offset+i] =  X[j][0] - X[i][0];
        Vy[offset+i] =  X[j][1] - X[i][1];
        Vz[offset+i] =  X[j][2] - X[i][2];
        Cx[offset+i] = (X[j][0] + X[i][0])/2;
        Cy[offset+i] = (X[j][1] + X[i][1])/2;
        Cz[offset+i] = (X[j][2] + X[i][2])/2;
      }
    };
    parallel_for(nset_faces, prepare);
    const std::vector<std::string> S={"x","y","z"};
    ExprReader reader(N*nset_faces, D);
    reader.register_variable(S[0], any(Reals(Cx)));
    reader.register_variable(S[1], any(Reals(Cy)));
    reader.register_variable(S[2], any(Reals(Cz)));
    auto result = reader.read_string(expr_string, name);
    reader.repeat(result);
    Reals field_osh = any_cast<Reals>(result);

    Write<Real> fld=field;
    auto save = OMEGA_H_LAMBDA(int set_face) {
      double Fl[N][3];
      double v [N][3];
      const auto face = faceLids[set_face];
      const auto offset = set_face*N;
      for (int i = 0; i < N; ++i)
        for (int s = 0; s < D; ++s)
          Fl[i][s] = field_osh[(offset+i) * D + s];
      for (int i = 0; i < N; ++i) {
        v[i][0] = Vx[offset+i];
        v[i][1] = Vy[offset+i];
        v[i][2] = Vz[offset+i];
      }
      double flux=0;
      for (int i = 0; i < N; ++i)
        for (int j = 0; j < D; ++j)
          flux += Fl[i][j]*v[i][j];

      fld[face] = flux;
    };
    parallel_for(nset_faces, save);
  }

  for (const std::string &blockname : element_blocks) {
    auto esIter = elemsets.find(blockname);
    if(esIter == elemsets.end())
       fail("Element block doesn't exist!");
    const auto elementLids = (esIter->second);
    const int   nset_elems = elementLids.size();
    const int  nset_faces = F*nset_elems;
    // Edge Centers
    Write<double> Cx(N*nset_faces);
    Write<double> Cy(N*nset_faces);
    Write<double> Cz(N*nset_faces);
    // Edge Vectors
    Write<double> Vx(N*nset_faces);
    Write<double> Vy(N*nset_faces);
    Write<double> Vz(N*nset_faces);
    const Reals cord = coords;
    const LOs   ids  = face_node_id;
    const LOs   fids = elem_face_id;
    auto prepare = OMEGA_H_LAMBDA(int set_face) {
      const int  el = elementLids[set_face/F];
      const int   f = set_face%F;
      const auto face = fids[F*el+f];
      int nodes[N];
      for (int i = 0; i < N; ++i)
        nodes[i] = ids[N*face+i];
      double X[N][3];
      for (int i = 0; i < N; ++i)
        for (int j = 0; j < D; ++j)
           X[i][j] = cord[D*nodes[i]+j];

      const auto offset = set_face*N;
      for (int i = 0; i < N; ++i) {
        const int j = (i+1)%N;
        Vx[offset+i] =  X[j][0] - X[i][0];
        Vy[offset+i] =  X[j][1] - X[i][1];
        Vz[offset+i] =  X[j][2] - X[i][2];
        Cx[offset+i] = (X[j][0] + X[i][0])/2;
        Cy[offset+i] = (X[j][1] + X[i][1])/2;
        Cz[offset+i] = (X[j][2] + X[i][2])/2;
      }
    };
    parallel_for(nset_faces, prepare);
    const std::vector<std::string> S={"x","y","z"};
    ExprReader reader(N*nset_faces, D);
    reader.register_variable(S[0], any(Reals(Cx)));
    reader.register_variable(S[1], any(Reals(Cy)));
    reader.register_variable(S[2], any(Reals(Cz)));
    auto result = reader.read_string(expr_string, name);
    reader.repeat(result);
    Reals field_osh = any_cast<Reals>(result);

    Write<Real> fld=field;
    auto save = OMEGA_H_LAMBDA(int set_face) {
      double Fl[N][3];
      double v [N][3];
      const int  el = elementLids[set_face/F];
      const int   f = set_face%F;
      const auto face = fids[f*el+f];
      const auto offset = set_face*N;
      for (int i = 0; i < N; ++i)
        for (int s = 0; s < D; ++s)
          Fl[i][s] = field_osh[(offset+i) * D + s];
      for (int i = 0; i < N; ++i) {
        v[i][0] = Vx[offset+i];
        v[i][1] = Vy[offset+i];
        v[i][2] = Vz[offset+i];
      }
      double flux=0;
      for (int i = 0; i < N; ++i)
        for (int j = 0; j < D; ++j)
          flux += Fl[i][j]*v[i][j];

      fld[face] = flux;
    };
    parallel_for(nset_faces, save);
  }
}

int main(int argc, char** argv) {
  Library lib(&argc, &argv);
  CommPtr world = lib.world();
  Mesh mesh = build_box(world, OMEGA_H_SIMPLEX, 1., 2., 3., 2, 2, 2);
  AdaptOpts opts(&mesh);
  mesh.add_tag<Real>(VERT, "metric", 1);


  Reals coords  = mesh.coords();
  Write<Real> flux_w(mesh.nfaces());

  Assoc assoc = get_box_assoc(3);
  const MeshSets mesh_sets = invert(&mesh, assoc);
  const MeshDimSets facesets = mesh_sets[SIDE_SET];
  const MeshDimSets elemsets = mesh_sets[ELEM_SET];
  LOs face_node_id = mesh.get_adj(FACE, VERT).ab2b;
  LOs elem_face_id = mesh.get_adj(REGION, FACE).ab2b;

  FunctionInitialCondition(
    facesets,
    elemsets,
    coords,
    face_node_id,
    elem_face_id,
    flux_w);

  {
    auto f = OMEGA_H_LAMBDA(LO v) {
      Vector<3> x = get_vector<3>(coords, v);
      Vector<3> y = {{2*x[0], 5*x[1], x[2]}};
      for (int i=0; i<3; ++i) flux_w[v] = y[0];
    };
    parallel_for(mesh.nfaces(), f);
  }

  const std::string name = "magnetic face flux";
  Reals flux(flux_w);
  mesh.add_tag<Real>(FACE, name, 1, flux);

  mesh.set_tag(
      VERT, "metric", Reals(mesh.nverts(), metric_eigenvalue_from_length(1.3)));
  while (coarsen_by_size(&mesh, opts))
    ;

  mesh.set_tag(
      VERT, "metric", Reals(mesh.nverts(), metric_eigenvalue_from_length(1.3)));
  while (refine_by_size(&mesh, opts))
    ;

  Read<Real> flux_r=mesh.get_array<Real>(FACE, name); 

  Write<Real> OK(1.,0);
  {
    auto f = OMEGA_H_LAMBDA(LO v) {
      Real x = flux_r[v];
      Vector<3> y = {{2*x, 5*x, x}};
      if (x != y[0]) OK[0] += 1;
      if (x != y[1]) OK[0] += 1;
      if (x != y[2]) OK[0] += 1;
    };
    parallel_for(mesh.nfaces(), f);
  }
  const bool ok = 0.==Reals(OK)[0];

  mesh.ask_qualities();
  if (ok) return 2;
  return 0;
}
