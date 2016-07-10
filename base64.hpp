#ifndef BASE64_HPP
#define BASE64_HPP

#include <istream>
#include <string>

namespace osh {

namespace base64 {

std::size_t encoded_size(std::size_t size);
std::string encode(void const* data, std::size_t size);
void decode(std::string const& text, void* data, std::size_t size);
std::string read_encoded(std::istream& f);
}

}  // end namespace osh

#endif
