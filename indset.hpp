namespace indset {

enum {
  UNKNOWN,
  IN,
  NOT_IN
};

Read<I8> iteration(
    LOs xadj, LOs adj,
    Reals quality,
    Read<GO> global,
    Read<I8> old_state);

}
