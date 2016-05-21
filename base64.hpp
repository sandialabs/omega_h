namespace base64 {

std::string encode(void const* data, std::size_t size);
void decode(std::string const& text, void* data, std::size_t size);
std::string read_encoded(std::istream& f);

}

