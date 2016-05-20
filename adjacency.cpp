Adj::Adj() {
}

Adj::Adj(LOs ab2b_):
  ab2b(ab2b_) {
}

Adj::Adj(LOs ab2b_, Read<I8> codes_):
  ab2b(ab2b_),codes(codes_) {
}

Adj::Adj(LOs a2ab_, LOs ab2b_, Read<I8> codes_):
  a2ab(a2ab_),ab2b(ab2b_),codes(codes_) {
}

Adj::Adj(LOs a2ab_, LOs ab2b_):
  a2ab(a2ab_),ab2b(ab2b_) {
}
