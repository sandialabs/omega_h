Int const simplex_degrees[DIMS][DIMS] = {
  {1,-1,-1,-1},
  {2, 1,-1,-1},
  {3, 3,-1,-1},
  {4, 6, 4, 1}
};

char const* const singular_names[DIMS] = {
  "vertex",
  "edge",
  "triangle",
  "tet"
};

char const* const plural_names[DIMS] = {
  "vertices",
  "edges",
  "triangles",
  "tets"
};
