#include <Omega_h_file.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_class.hpp>
#include <Omega_h_cmdline.hpp>
#include <Omega_h_for.hpp>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  Omega_h::CmdLine cmdline;
  cmdline.add_arg<std::string>("input.osh");
  cmdline.add_arg<std::string>("output.osh");
  if (!cmdline.parse_final(lib.world(), &argc, argv)) return -1;
  auto path_in = cmdline.get<std::string>("input.osh");
  auto path_out = cmdline.get<std::string>("output.osh");
  auto mesh = Omega_h::binary::read(path_in, &lib, /*strict=*/true);
  mesh.set_parting(OMEGA_H_GHOSTED);
  auto faces_to_elements = mesh.ask_up(2, 3);
  auto faces_to_face_elements = faces_to_elements.a2ab;
  auto face_elements_to_elements = faces_to_elements.ab2b;
  auto element_class_ids = mesh.get_array<Omega_h::ClassId>(3, "class_id");
  auto face_class_dims = Omega_h::deep_copy(mesh.get_array<Omega_h::Byte>(2, "class_dim"));
  auto face_class_ids = Omega_h::deep_copy(mesh.get_array<Omega_h::ClassId>(2, "class_id"));
  auto f = OMEGA_H_LAMBDA(Omega_h::LO face) {
    auto first = faces_to_face_elements[face];
    auto last = faces_to_face_elements[face + 1];
    int touches_air = 0;
    int touches_object = 0;
    for (auto face_element = first; face_element < last;
        ++face_element) {
      auto element = face_elements_to_elements[face_element];
      auto class_id = element_class_ids[element];
      if (class_id == 0) ++touches_air;
      if (class_id == 1) ++touches_object;
    }
    int dim = 3;
    int id = -1;
    if (touches_air == 1 && touches_object == 1) {
      dim = 2;
      id = 3;
    }
    if (touches_air == 1 && touches_object == 0) {
      dim = 2;
      id = 1;
    }
    if (touches_air == 0 && touches_object == 1) {
      dim = 2;
      id = 2;
    }
    face_class_dims[face] = dim;
    face_class_ids[face] = id;
  };
  Omega_h::parallel_for(mesh.nfaces(), f, "classify_Morph_faces");
  mesh.add_tag<Omega_h::Byte>(2, "class_dim", 1, face_class_dims);
  mesh.add_tag<Omega_h::ClassId>(2, "class_id", 1, face_class_ids);
  Omega_h::finalize_classification(&mesh);
  Omega_h::binary::write(path_out, &mesh);
}
