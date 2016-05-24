bool is_little_endian_cpu();

namespace binary {

template <typename T>
void write_value(std::ostream& stream, T val);
template <typename T>
void read_value(std::istream& stream, T& val);
template <typename T>
void write_array(std::ostream& stream, Read<T> array);
template <typename T>
void read_array(std::istream& stream, Read<T>& array,
    bool is_compressed);
void write(std::ostream& stream, std::string const& val);
void read(std::istream& stream, std::string& val);

void write(std::ostream& stream, Mesh& mesh);
void read(std::istream& stream, Mesh& mesh);

}
