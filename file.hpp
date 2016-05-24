bool is_little_endian_cpu();

namespace file {

template <typename T>
void write_binary(std::ostream& stream, T val);
template <typename T>
void read_binary(std::istream& stream, T& val);
template <typename T>
void write_binary(std::ostream& stream, Read<T> array);
template <typename T>
void read_binary(std::istream& stream, Read<T>& array,
    bool is_compressed);

}
