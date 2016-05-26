class Dist {
  CommPtr parent_comm_;
  LOs roots2items_[2];
  LOs items2content_[2];
  LOs msgs2content_[2];
  Read<I32> msgs2ranks_[2];
  CommPtr comm_[2];
public:
  Dist();
  Dist(CommPtr comm, Remotes fitems2rroots, LO nrroots);
  void set_parent_comm(CommPtr parent_comm);
  void set_dest_ranks(Read<I32> items2ranks);
  void set_dest_idxs(LOs fitems2rroots, LO nrroots);
  void set_roots2items(LOs froots2fitems);
  Dist invert() const;
  template <typename T>
  Read<T> exch(Read<T> data, Int width) const;
  CommPtr parent_comm() const;
  CommPtr comm() const;
  LOs content2msgs() const;
  LOs items2msgs() const;
  LOs roots2items() const;
  Read<I32> msgs2ranks() const;
  Read<I32> items2ranks() const;
  LO nitems() const;
  LO nroots() const;
private:
  enum { F, R };
};
