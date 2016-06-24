namespace xml {

struct StartTag {
  std::string elem_name;
  std::map<std::string, std::string> attribs;
};

bool parse_start_tag(std::string const& line,
    StartTag* tag_out);

}
