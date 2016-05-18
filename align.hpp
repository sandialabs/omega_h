INLINE I8 make_code(bool is_flipped, I8 rotation, I8 which_down) {
  return static_cast<I8>((which_down << 3) | (rotation << 1) | is_flipped);
}

INLINE bool code_is_flipped(I8 code) {
  return code & 1;
}

INLINE I8 code_rotation(I8 code) {
  return (code >> 1) & 3;
}

INLINE I8 code_which_down(I8 code) {
  return (code >> 3);
}

