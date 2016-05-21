namespace vtk {

namespace {

/* start of C++ ritual dance to print a string based on
   type properties */

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-member-function"
#endif

template <bool is_signed, std::size_t size>
struct IntTraits;

template <>
struct IntTraits<true,1> {
  inline static char const* name() { return "Int8"; }
};

template <>
struct IntTraits<true,4> {
  inline static char const* name() { return "Int32"; }
};

template <>
struct IntTraits<false,4> {
  inline static char const* name() { return "UInt32"; }
};

template <>
struct IntTraits<true,8> {
  inline static char const* name() { return "Int64"; }
};

template <>
struct IntTraits<false,8> {
  inline static char const* name() { return "UInt64"; }
};

template <std::size_t size>
struct FloatTraits;

template <>
struct FloatTraits<4> {
  inline static char const* name() { return "Float32"; }
};

template <>
struct FloatTraits<8> {
  inline static char const* name() { return "Float64"; }
};

template <typename T, typename Enable = void>
struct Traits;

template <typename T>
struct Traits<T,
  typename std::enable_if<std::is_integral<T>::value>::type> {
  inline static char const* name() {
    return IntTraits<std::is_signed<T>::value,sizeof(T)>::name();
  }
};

template <typename T>
struct Traits<T,
  typename std::enable_if<std::is_floating_point<T>::value>::type> {
  inline static char const* name() {
    return FloatTraits<sizeof(T)>::name();
  }
};

#ifdef __clang__
#pragma clang diagnostic pop
#endif

/* end of C++ ritual dance to get a string based on type properties */

template <typename T>
void describe_array(std::ostream& stream, std::string const& name,
    Int ncomps)
{
  stream << "type=\"" << Traits<T>::name() << "\"";
  stream << " Name=\"" << name << "\"";
  stream << " NumberOfComponents=\"" << ncomps << "\"";
  stream << " format=\"binary\"";
}

template <typename T>
void write_array(std::ostream& stream, std::string const& name,
    Int ncomps, Read<T> array)
{
  stream << "<DataArray ";
  describe_array<T>(stream, name, ncomps);
  stream << ">\n";
  HostRead<T> uncompressed(array);
  std::size_t uncompressed_bytes = sizeof(T) *
    static_cast<std::size_t>(array.size());
#ifdef USE_ZLIB
  uLong source_bytes = uncompressed_bytes;
  uLong dest_bytes = ::compressBound(source_bytes);
  Bytef* compressed = new Bytef[dest_bytes];
  int ret = ::compress2(
      compressed,
      &dest_bytes,
      reinterpret_cast<const Bytef*>(&uncompressed[0]),
      source_bytes,
      9); //compress to the max
  CHECK(ret == Z_OK);
  std::string encoded = base64::encode(compressed, dest_bytes);
  delete [] compressed;
  std::size_t header[4] = {
    1,
    uncompressed_bytes,
    uncompressed_bytes,
    dest_bytes
  };
  std::string enc_header = base64::encode(header, sizeof(header));
#else
  std::string enc_header = base64::encode(&uncompressed_bytes,
      sizeof(std::size_t));
  std::string encoded = base64::encode(&uncompressed[0],
      uncompressed_bytes);
#endif
  stream << enc_header << encoded << '\n';
  stream << "</DataArray>\n";
}

void write_tag(std::ostream& stream, TagBase const* tag, Int space_dim)
{
  if (is<I8>(tag))
    write_array(stream, tag->name(), tag->ncomps(), to<I8>(tag)->array());
  else if (is<I32>(tag))
    write_array(stream, tag->name(), tag->ncomps(), to<I32>(tag)->array());
  else if (is<I64>(tag))
    write_array(stream, tag->name(), tag->ncomps(), to<I64>(tag)->array());
  else if (is<Real>(tag)) {
    Reals array = to<Real>(tag)->array();
    // VTK / Paraview expect vector fields to have 3 components
    // regardless of whether this is a 2D mesh or not.
    // this ad-hoc filter adds a 3rd zero component to any
    // fields with 2 components for 2D meshes
    if (space_dim == 2 && tag->ncomps() == space_dim) {
      Int np = array.size() / 2;
      Write<Real> tmp(np * 3);
      auto f = LAMBDA(Int i) {
        tmp[i * 3 + 0] = array[i * 2 + 0];
        tmp[i * 3 + 1] = array[i * 2 + 1];
        tmp[i * 3 + 2] = 0.0;
      };
      parallel_for(np, f);
      write_array(stream, tag->name(), 3, Reals(tmp));
    } else {
      write_array(stream, tag->name(), tag->ncomps(), array);
    }
  } else
    fail("unknown tag type in write_tag");
}

enum Type {
  VTK_VERTEX         = 1,
  VTK_POLY_VERTEX    = 2,
  VTK_LINE           = 3,
  VTK_POLY_LINE      = 4,
  VTK_TRIANGLE       = 5,
  VTK_TRIANGLE_STRIP = 6,
  VTK_POLYGON        = 7,
  VTK_PIXEL          = 8,
  VTK_QUAD           = 9,
  VTK_TETRA          =10,
  VTK_VOXEL          =11,
  VTK_HEXAHEDRON     =12,
  VTK_WEDGE          =13,
  VTK_PYRAMID        =14
};

static I8 const vtk_types[DIMS] = {
  VTK_VERTEX,
  VTK_LINE,
  VTK_TRIANGLE,
  VTK_TETRA
};

static void write_vtkfile_start_tag(std::ostream& stream)
{
  stream << "<VTKFile type=\"UnstructuredGrid\" byte_order=\"";
  if (is_little_endian_cpu())
    stream << "LittleEndian";
  else
    stream << "BigEndian";
  stream << "\" header_type=\"";
  stream << Traits<std::size_t>::name();
  stream << "\"";
#ifdef USE_ZLIB
  stream << " compressor=\"vtkZLibDataCompressor\"";
#endif
  stream << ">\n";
}

void write_piece_start_tag(std::ostream& stream, Mesh const& mesh, Int cell_dim)
{
  stream << "<Piece NumberOfPoints=\"" << mesh.nverts() << "\"";
  stream << " NumberOfCells=\"" << mesh.nents(cell_dim) << "\">\n";
}

void write_connectivity(std::ostream& stream, Mesh& mesh, Int cell_dim)
{
  LOs ev2v = mesh.ask_adj(cell_dim, VERT).ab2b;
  LOs ends(mesh.nents(cell_dim), simplex_degrees[cell_dim][VERT],
                                 simplex_degrees[cell_dim][VERT]);
  write_array(stream, "connectivity", 1, ev2v);
  write_array(stream, "offsets", 1, ends);
  Read<I8> types(mesh.nents(cell_dim), vtk_types[cell_dim]);
  write_array(stream, "types", 1, types);
}

}//end anonymous namespace

void write_vtu(std::ostream& stream, Mesh& mesh, Int cell_dim) {
  write_vtkfile_start_tag(stream);
  stream << "<UnstructuredGrid>\n";
  write_piece_start_tag(stream, mesh, cell_dim);
  stream << "<Points>\n";
  write_tag(stream, &(mesh.get_tag<Real>(VERT,"coordinates")), mesh.dim());
  stream << "</Points>\n";
  stream << "<Cells>\n";
  write_connectivity(stream, mesh, cell_dim);
  stream << "</Cells>\n";
  stream << "</Piece>\n";
  stream << "</UnstructuredGrid>\n";
  stream << "</VTKFile>\n";
}

void write_vtu(std::ostream& stream, Mesh& mesh) {
  write_vtu(stream, mesh, mesh.dim());
}

}//end namespace vtk
