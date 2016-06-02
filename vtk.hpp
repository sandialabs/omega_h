namespace vtk {
void write_vtu(std::ostream& stream, Mesh& mesh, Int cell_dim);
void write_vtu(std::string const& filename, Mesh& mesh, Int cell_dim);
void write_pvtu(std::ostream& stream, Mesh& mesh, Int cell_dim,
    std::string const& piecepath);
void write_pvtu(std::string const& filename, Mesh& mesh, Int cell_dim,
    std::string const& piecepath);
void write_parallel(std::string const& path, Mesh& mesh, Int cell_dim);
}
