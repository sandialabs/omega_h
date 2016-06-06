namespace vtk {
void write_pvtu(std::ostream& stream, Mesh& mesh, Int cell_dim,
    std::string const& piecepath);
void write_pvtu(std::string const& filename, Mesh& mesh, Int cell_dim,
    std::string const& piecepath);
void write_pvd(std::string const& root_path, std::vector<Real> const& times);
}
