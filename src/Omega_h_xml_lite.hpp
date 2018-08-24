#ifndef XML_LITE_HPP
#define XML_LITE_HPP

#include <iosfwd>
#include <map>
#include <string>

namespace Omega_h {

namespace xml_lite {

struct Tag {
  std::string elem_name;
  enum { START, END, SELF_CLOSING } type;
  std::map<std::string, std::string> attribs;
};

bool parse_tag(std::string const& line, xml_lite::Tag* tag_out);

xml_lite::Tag read_tag(std::istream& stream);
}  // namespace xml_lite

}  // end namespace Omega_h

#endif
