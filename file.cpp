bool is_little_endian_cpu()
{
  static std::uint16_t const endian_canary = 0x1;
  std::uint8_t const* p = reinterpret_cast<std::uint8_t const*>(&endian_canary);
  return *p == 0x1;
}
