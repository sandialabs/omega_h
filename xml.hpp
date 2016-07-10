#ifndef XML_HPP
#define XML_HPP

namespace osh {

namespace xml {

struct Tag {
  std::string elem_name;
  enum { START, END, SELF_CLOSING } type;
  std::map<std::string, std::string> attribs;
};

bool parse_tag(std::string const& line, xml::Tag* tag_out);

xml::Tag read_tag(std::istream& stream);
}

} //end namespace osh

#endif
