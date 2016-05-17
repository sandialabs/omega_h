namespace indset {

Read<I8> iteration(
    LOs xadj, LOs adj,
    Reals quality,
    Read<GO> global,
    Read<I8> old_state) {
  LO n = global.size();
  Write<I8> new_state(old_state.size());
  auto f = LAMBDA(LO v) {
    if (old_state[v] != UNKNOWN)
      return;
    LO begin = xadj[v];
    LO end = xadj[v + 1];
    // nodes adjacent to chosen ones are rejected
    for (LO j = begin; j < end; ++j) {
      LO u = adj[j];
      if (old_state[u] == IN) {
        new_state[v] = NOT_IN;
        return;
      }
    }
    // check if node is a local maximum
    Real v_qual = quality[v];
    for (LO j = begin; j < end; ++j) {
      LO u = adj[j];
      // neighbor was rejected, ignore its presence
      if (old_state[u] == NOT_IN) continue;
      Real u_qual = quality[u];
      // neighbor has higher quality
      if (u_qual > v_qual) return;
      // neighbor has equal quality, tiebreaker by global ID
      if (u_qual == v_qual && global[u] > global[v]) return;
    }
    // only local maxima reach this line
    new_state[v] = IN;
  };
  parallel_for(n, f);
  return new_state;
}

}
