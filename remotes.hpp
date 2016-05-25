struct Remotes {
  Remotes() {}
  Remotes(Read<I32> ranks_, LOs idxs_):
    ranks(ranks_),idxs(idxs_) {
  }
  Read<I32> ranks;
  LOs idxs;
};
