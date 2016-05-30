Dist bi_partition(CommPtr comm, Read<I8> marks) {
  CHECK(comm->size() % 2 == 0);
  Write<I32> dest_ranks(marks.size());
  Write<LO> dest_idxs(marks.size());
  LO linsize = -1;
  auto halfsize = comm->size() / 2;
  I32 rank_start = 0;
  for (Int half = 0; half < 2; ++half) {
    marks = invert_marks(marks);
    auto marked = collect_marked(marks);
    auto total = comm->allreduce(GO(marked.size()), SUM);
    auto start = comm->exscan(GO(marked.size()), SUM);
    Read<GO> globals(marked.size(), start, 1);
    auto owners = globals_to_linear_owners(globals, total, halfsize);
    auto write = LAMBDA(LO i) {
      dest_ranks[marked[i]] = owners.ranks[i] + rank_start;
      dest_idxs[marked[i]] = owners.idxs[i];
    };
    parallel_for(marked.size(), write);
    if (rank_start <= comm->rank() && comm->rank() < (rank_start + halfsize)) {
      linsize = linear_partition_size(total, halfsize, comm->rank() - rank_start);
    }
    rank_start += halfsize;
  }
  auto dests = Remotes(Read<I32>(dest_ranks), Read<LO>(dest_idxs));
  return Dist(comm, dests, linsize);
}
