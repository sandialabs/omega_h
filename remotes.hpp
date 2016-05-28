struct Remotes {
  Remotes() {}
  Remotes(Read<I32> ranks_, LOs idxs_):
    ranks(ranks_),idxs(idxs_) {
  }
  Read<I32> ranks;
  LOs idxs;
};

Remotes expand(Remotes a2c, LOs a2b);
Remotes unmap(LOs a2b, Remotes b2c);
