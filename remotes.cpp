Remotes expand(Remotes a_copies2a_owners, LOs a2b) {
  return Remotes(expand(a_copies2a_owners.ranks, a2b, 1),
                 expand(a_copies2a_owners.idxs, a2b, 1));
}
