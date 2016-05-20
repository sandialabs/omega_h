/* start of C++ ritual dance to print a string based on
   type properties */
template <bool is_signed, std::size_t size>
struct VtkIntTraits;

template <>
struct VtkIntTraits<true,1> {
  inline static char const* name() { return "Int8"; }
};

template <>
struct VtkIntTraits<true,4> {
  inline static char const* name() { return "Int32"; }
};

template <>
struct VtkIntTraits<false,4> {
  inline static char const* name() { return "UInt32"; }
};

template <>
struct VtkIntTraits<true,8> {
  inline static char const* name() { return "Int64"; }
};

template <>
struct VtkIntTraits<false,8> {
  inline static char const* name() { return "UInt64"; }
};

template <std::size_t size>
struct VtkFloatTraits;

template <>
struct VtkFloatTraits<4> {
  inline static char const* name() { return "Float32"; }
};

template <>
struct VtkFloatTraits<8> {
  inline static char const* name() { return "Float64"; }
};

template <typename T, typename Enable = void>
struct VtkTraits;

template <typename T>
struct VtkTraits<T,
  typename std::enable_if<std::is_integral<T>::value>::type> {
  inline static char const* name() {
    return VtkIntTraits<std::is_signed<T>::value,sizeof(T)>::name();
  }
};

template <typename T>
struct VtkTraits<T,
  typename std::enable_if<std::is_floating_point<T>::value>::type> {
  inline static char const* name() {
    return VtkFloatTraits<sizeof(T)>::name();
  }
};

/* end of C++ ritual dance to get a string based on type properties */

static void write_vtkfile_start_tag(std::ostream& stream)
{
  stream << "<VTKFile type=\"UnstructuredGrid\" byte_order=\"";
  if (is_little_endian_cpu())
    stream << "LittleEndian";
  else
    stream << "BigEndian";
  stream << "\" header_type=\"";
  stream << VtkTraits<std::size_t>::name();
  stream << "\"";
#ifdef USE_ZLIB
  stream << " compressor=\"vtkZLibDataCompressor\"";
#endif
  stream << ">\n";
}

void write_vtu(std::ostream& stream, Mesh const& mesh) {
  write_vtkfile_start_tag(stream);
  (void) mesh;
}
