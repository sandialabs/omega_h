#include "Omega_h_xml_lite.hpp"

#include "Omega_h_defines.hpp"
#include "Omega_h_fail.hpp"

#include <istream>

namespace Omega_h {

namespace xml_lite {

bool parse_tag(std::string const& line, xml_lite::Tag* tag_out) {
  /* state machine transitions:
    </ Foo   bar="quux" />
  --||-|--|--|--||----|-||-
  0 12 3  4  5  67    3 89   */
  xml_lite::Tag tag;
  Int state = 0;
  std::string att_nm;
  std::string att_val;
  char quot = '\0';
  tag.type = xml_lite::Tag::START;
  for (auto c : line) {
    switch (state) {
      case 0:
        if (c == '<')
          state = 1;
        else if (c != ' ')
          return false;
        break;
      case 1:
        if (c == '/') {
          tag.type = xml_lite::Tag::END;
          state = 2;
        } else if (c != ' ') {
          tag.elem_name.push_back(c);
          state = 3;
        }
        break;
      case 2:
        if (c != ' ') {
          state = 3;
          tag.elem_name.push_back(c);
        }
        break;
      case 3:
        if (c == ' ')
          state = 4;
        else if (c == '/')
          state = 8;
        else if (c == '>')
          state = 9;
        else
          tag.elem_name.push_back(c);
        break;
      case 4:
        if (c == '/') {
          tag.type = xml_lite::Tag::SELF_CLOSING;
          state = 8;
        } else if (c == '>') {
          state = 9;
        } else if (c != ' ') {
          state = 5;
          att_nm.push_back(c);
        }
        break;
      case 5:
        if (c == '=')
          state = 6;
        else
          att_nm.push_back(c);
        break;
      case 6:
        if (c == '\"' || c == '\'') {
          quot = c;
          state = 7;
        }
        break;
      case 7:
        if (c == quot) {
          tag.attribs[att_nm] = att_val;
          att_nm = att_val = std::string();
          state = 4;
        } else {
          att_val.push_back(c);
        }
        break;
      case 8:
        if (c == '>') {
          state = 9;
        }
        break;
    }
  }
  if (state == 9) {
    *tag_out = tag;
    return true;
  }
  return false;
}

xml_lite::Tag read_tag(std::istream& stream) {
  std::string line;
  std::getline(stream, line);
  xml_lite::Tag st;
  OMEGA_H_CHECK(parse_tag(line, &st));
  return st;
}
}  // namespace xml_lite

}  // end namespace Omega_h
