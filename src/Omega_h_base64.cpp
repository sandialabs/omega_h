#include "Omega_h_base64.hpp"

#include "Omega_h_fail.hpp"

namespace Omega_h {

namespace base64 {

namespace {

char const value_to_char[] = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',
    'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y',
    'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n',
    'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '0', '1', '2',
    '3', '4', '5', '6', '7', '8', '9', '+', '/'};

unsigned char const char_to_value[256] = {255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 62, 255, 255, 255, 63, 52, 53, 54, 55, 56, 57,
    58, 59, 60, 61, 255, 255, 255, 0, 255, 255,

    255, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
    20, 21, 22, 23, 24, 25, 255, 255, 255, 255, 255, 255, 26, 27, 28, 29, 30,
    31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
    50, 51, 255, 255, 255, 255, 255,

    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255,

    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255};

void encode_3(unsigned char const* in, char* out) {
  out[0] = value_to_char[in[0] >> 2];
  out[1] = value_to_char[((in[0] << 4) & 0x30) | (in[1] >> 4)];
  out[2] = value_to_char[((in[1] << 2) & 0x3C) | (in[2] >> 6)];
  out[3] = value_to_char[in[2] & 0x3F];
}

void encode_2(unsigned char const* in, char* out) {
  out[0] = value_to_char[in[0] >> 2];
  out[1] = value_to_char[((in[0] << 4) & 0x30) | (in[1] >> 4)];
  out[2] = value_to_char[(in[1] << 2) & 0x3C];
  out[3] = '=';
}

void encode_1(unsigned char const* in, char* out) {
  out[0] = value_to_char[in[0] >> 2];
  out[1] = value_to_char[(in[0] << 4) & 0x30];
  out[2] = '=';
  out[3] = '=';
}

template <class T>
unsigned U(T val) {
  return static_cast<unsigned>(val);
}

template <class T>
unsigned char UC(T val) {
  return static_cast<unsigned char>(val);
}

void decode_4(char const* in, unsigned char* out, std::size_t nout = 3) {
  unsigned char val[4];
  for (unsigned i = 0; i < 4; ++i) {
    val[i] = char_to_value[U(in[i])];
    OMEGA_H_CHECK(val[i] < 64);
  }
  /* cast it all !
   * (technically this should all be run as
   *  unsigned char, but apparently it gets pushed
   *  up to int if we try that, so run it all as unsigned
   *  and cast back which we know is safe) */
  out[0] =
      UC(((U(val[0]) << U(2)) & U(0xFC)) | ((U(val[1]) >> U(4)) & U(0x03)));
  if (nout < 2) return;
  out[1] =
      UC(((U(val[1]) << U(4)) & U(0xF0)) | ((U(val[2]) >> U(2)) & U(0x0F)));
  if (nout < 3) return;
  out[2] =
      UC(((U(val[2]) << U(6)) & U(0xC0)) | ((U(val[3]) >> U(0)) & U(0x3F)));
}

}  // end anonymous namespace

std::size_t encoded_size(std::size_t size) {
  auto quot = size / 3;
  auto rem = size % 3;
  auto nunits = rem ? (quot + 1) : quot;
  auto nchars = nunits * 4;
  return nchars;
}

std::string encode(void const* data, std::size_t size) {
  auto quot = size / 3;
  auto rem = size % 3;
  auto nunits = rem ? (quot + 1) : quot;
  auto nchars = nunits * 4;
  std::string out(nchars, '\0');
  unsigned char const* in = static_cast<unsigned char const*>(data);
  for (std::size_t i = 0; i < quot; ++i) encode_3(&in[i * 3], &out[i * 4]);
  switch (rem) {
    case 0:
      break;
    case 1:
      encode_1(&in[quot * 3], &out[quot * 4]);
      break;
    case 2:
      encode_2(&in[quot * 3], &out[quot * 4]);
      break;
  }
  return out;
}

void decode(std::string const& text, void* data, std::size_t size) {
  std::size_t quot = size / 3;
  std::size_t rem = size % 3;
  unsigned char* out = static_cast<unsigned char*>(data);
  for (std::size_t i = 0; i < quot; ++i) decode_4(&text[i * 4], &out[i * 3]);
  if (rem) decode_4(&text[quot * 4], &out[quot * 3], rem);
}

std::string read_encoded(std::istream& f) {
  std::string out;
  while (true) {
    int c = f.get();
    if (c < 0 || c > 127) break;
    unsigned char val = char_to_value[c];
    if (val > 63) break;
    out.push_back(static_cast<char>(c));
  }
  return out;
}

}  // end namespace base64

}  // end namespace Omega_h
