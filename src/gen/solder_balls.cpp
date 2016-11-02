#include <gmodel.hpp>
#include <cstdlib>
#include <cstring>
#include <sstream>

static gmod::ObjPtr new_solder_ball(gmod::Vector center,
    double radius, double height,
    double mid_res, double edge_res) {
  double half_height = height / 2.0;
  double disk_radius = sqrt(radius * radius - half_height * half_height);
  gmod::default_size = edge_res;
  auto bot_disk = gmod::new_disk(center - gmod::Vector{0,0,half_height},
      gmod::Vector{0,0,1}, gmod::Vector{disk_radius,0,0});
  auto top_disk = gmod::new_disk(center + gmod::Vector{0,0,half_height},
      gmod::Vector{0,0,1}, gmod::Vector{disk_radius,0,0});
  gmod::default_size = mid_res;
  auto mid_circ = gmod::new_circle(center, gmod::Vector{0,0,1},
      gmod::Vector{radius,0,0});
  auto bot_loop = gmod::face_loop(bot_disk);
  auto top_loop = gmod::face_loop(top_disk);
  auto bot_pts = gmod::loop_points(bot_loop);
  auto top_pts = gmod::loop_points(top_loop);
  auto mid_pts = gmod::loop_points(mid_circ);
  auto center_pt = gmod::new_point2(center);
  gmod::ObjPtr verticals[2][4];
  for (size_t i = 0; i < 4; ++i) {
    verticals[0][i] = gmod::new_arc2(bot_pts[i], center_pt, mid_pts[i]);
    verticals[1][i] = gmod::new_arc2(mid_pts[i], center_pt, top_pts[i]);
  }
  auto shell = gmod::new_shell();
  gmod::add_use(shell, gmod::REVERSE, bot_disk);
  gmod::add_use(shell, gmod::FORWARD, top_disk);
  for (size_t i = 0; i < 4; ++i) {
    auto loop = gmod::new_loop();
    gmod::add_use(loop, gmod::FORWARD, bot_loop->used[i].obj);
    gmod::add_use(loop, gmod::FORWARD, verticals[0][(i + 1) % 4]);
    gmod::add_use(loop, gmod::REVERSE, verticals[0][i]);
    gmod::add_use(loop, gmod::REVERSE, mid_circ->used[i].obj);
    auto side = gmod::new_ruled2(loop);
    gmod::add_use(shell, gmod::FORWARD, side);
  }
  for (size_t i = 0; i < 4; ++i) {
    auto loop = gmod::new_loop();
    gmod::add_use(loop, gmod::FORWARD, mid_circ->used[i].obj);
    gmod::add_use(loop, gmod::FORWARD, verticals[1][(i + 1) % 4]);
    gmod::add_use(loop, gmod::REVERSE, verticals[1][i]);
    gmod::add_use(loop, gmod::REVERSE, top_loop->used[i].obj);
    auto side = gmod::new_ruled2(loop);
    gmod::add_use(shell, gmod::FORWARD, side);
  }
  return gmod::new_volume2(shell);
}

static gmod::ObjPtr solder_ball_bot(gmod::ObjPtr solder_ball) {
  return solder_ball->used[0].obj->used[0].obj;
}

static gmod::ObjPtr solder_ball_top(gmod::ObjPtr solder_ball) {
  return solder_ball->used[0].obj->used[1].obj;
}

int main(int argc, char** argv) {
  int balls_per_side = 1;
  double length_unit = 1.0;
  double block_height = 1.2;
  double ball_height = 0.6;
  double ball_diameter = 0.8;
  double block_length = length_unit * balls_per_side;
  double outer_res = length_unit;
  double inner_res = length_unit / 6.0;
  double edge_res = 0.06;
  double mid_res = 0.12;
  const char* out_path = nullptr;
  for (int i = 1; i < argc; ++i) {
    if (!strcmp("-n", argv[i])) {
      balls_per_side = atoi(argv[++i]);
    } else
    if (!strcmp("-u", argv[i])) {
      inner_res = outer_res = edge_res = mid_res = atof(argv[++i]);
      block_height = ball_height;
    } else
    if (!strcmp("-o", argv[i])) {
      out_path = argv[++i];
    }
  }
  gmod::default_size = outer_res;
  auto bot_block = gmod::new_cube(
      gmod::Vector{0, 0, 0},
      gmod::Vector{block_length, 0, 0},
      gmod::Vector{0, block_length, 0},
      gmod::Vector{0, 0, block_height});
  auto bot_face = get_cube_face(bot_block, gmod::TOP);
  auto top_block = gmod::new_cube(
      gmod::Vector{0, 0, block_height + ball_height},
      gmod::Vector{block_length, 0, 0},
      gmod::Vector{0, block_length, 0},
      gmod::Vector{0, 0, block_height});
  auto top_face = get_cube_face(top_block, gmod::BOTTOM);
  auto bot_pts = loop_points(face_loop(bot_face));
  auto top_pts = loop_points(face_loop(top_face));
  for (auto pt : bot_pts) pt->size = inner_res;
  for (auto pt : top_pts) pt->size = inner_res;
  auto model = gmod::new_group();
  gmod::add_to_group(model, top_block);
  gmod::add_to_group(model, bot_block);
  for (int i = 0; i < balls_per_side; ++i)
  for (int j = 0; j < balls_per_side; ++j) {
    auto sb = new_solder_ball(
        gmod::Vector{length_unit / 2.0 + i * length_unit,
                     length_unit / 2.0 + j * length_unit,
                     block_height + ball_height / 2.0},
        ball_diameter / 2.0,
        ball_height,
        mid_res, edge_res);
    gmod::weld_volume_face_into(bot_block, sb,
        bot_face, solder_ball_bot(sb));
    gmod::weld_volume_face_into(top_block, sb,
        top_face, solder_ball_top(sb));
    gmod::add_to_group(model, sb);
  }
  std::stringstream stream;
  if (out_path) {
    stream << out_path;
  } else {
    stream << "solder_balls_" << balls_per_side
      << "x" << balls_per_side << ".geo";
  }
  auto s = stream.str();
  gmod::write_closure_to_geo(model, s.c_str());
}

