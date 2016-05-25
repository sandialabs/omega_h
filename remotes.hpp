struct Remotes {
  Remotes() {}
  Remotes(Read<I32> ranks_, Read<I32> idxs_):
    ranks(ranks_),idxs(idxs_) {
  }
  Read<I32> ranks;
  Read<I32> idxs;
};
