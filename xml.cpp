namespace xml {

bool parse_start_tag(std::string const& line,
    StartTag* tag_out) {
/* state machine transitions:
  < Foo   bar="quux" >
--|-|--|--|--||----|-|--
0 1 2  3  4  56    3 7         */
  StartTag tag;
  Int state = 0;
  std::string att_nm;
  std::string att_val;
  char quot = '\0';
  for (auto c : line) {
    switch (state) {
    case 0:
      if (c == '<') state = 1;
      else if (c != ' ') return false;
      break;
    case 1:
      if (c != ' ') {
        state = 2;
        tag.elem_name.push_back(c);
      }
      break;
    case 2:
      if (c == ' ') state = 3;
      else tag.elem_name.push_back(c);
      break;
    case 3:
      if (c == '/' || c == '>') {
        state = 7;
      } else if (c != ' ') {
        state = 4;
        att_nm.push_back(c);
      }
      break;
    case 4:
      if (c == '=') state = 5;
      else att_nm.push_back(c);
      break;
    case 5:
      if (c == '\"' || c == '\'') {
        quot = c;
        state = 6;
      }
      break;
    case 6:
      if (c == quot) {
        tag.attribs[att_nm] = att_val;
        att_nm = att_val = std::string();
        state = 3;
      } else {
        att_val.push_back(c);
      }
      break;
    }
  }
  if (state == 7) {
    *tag_out = tag;
    return true;
  }
  return false;
}

}
