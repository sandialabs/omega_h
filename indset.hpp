namespace indset {

enum {
  NOT_IN,
  IN,
  UNKNOWN
};

Read<I8> local_iteration(
    LOs xadj, LOs adj,
    Reals quality,
    Read<GO> global,
    Read<I8> old_state);

Read<I8> iteration(
    Dist owners2copies,
    LOs xadj, LOs adj,
    Reals quality,
    Read<GO> global,
    Read<I8> old_state);

Read<I8> find(
    Dist owners2copies,
    LOs xadj, LOs adj,
    Reals quality,
    Read<GO> global,
    Read<I8> candidates);

}

Read<I8> find_indset(
    Mesh& mesh,
    Int ent_dim,
    Reals quality,
    Read<I8> candidates);
