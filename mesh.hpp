class Mesh {
  public:
    Mesh();
    void set_dim(I8 dim);
    void set_verts(LO nverts);
    I8 dim() const;
    LO nents(I8 dim) const;
    template <typename T>
    void add_tag(I8 dim, std::string const& name, I8 ncomps, Read<T> data);
    void remove_tag(I8 dim, std::string const& name);
    bool has_tag(I8 dim, std::string const& name) const;
    template <typename T>
    Tag<T> const& get_tag(I8 dim, std::string const& name) const;
    I8 count_tags(I8 dim) const;
    TagBase const* get_tag(I8 dim, I8 i) const;
    bool has_ents(I8 dim) const;
    void set_ents(I8 dim, Adj down);
    bool has_adj(I8 from, I8 to) const;
    Adj get_adj(I8 from, I8 to) const;
    Adj ask_adj(I8 from, I8 to);
    Read<GO> ask_globals(I8 dim);
  private:
    typedef std::shared_ptr<TagBase> TagPtr;
    typedef std::shared_ptr<Adj> AdjPtr;
    typedef std::vector<TagPtr> TagVector;
    typedef TagVector::iterator TagIter;
    typedef TagVector::const_iterator TagCIter;
    template <typename T>
    void add_tag_priv(I8 dim, std::string const& name, I8 ncomps, Read<T> data);
    void remove_tag_priv(I8 dim, std::string const& name);
    TagIter tag_iter(I8 dim, std::string const& name);
    TagCIter tag_iter(I8 dim, std::string const& name) const;
    void check_dim(I8 dim) const;
    void check_dim2(I8 dim) const;
    void add_adj(I8 from, I8 to, Adj adj);
    Adj derive_adj(I8 from, I8 to);
    I8 dim_;
    LO nents_[DIMS];
    TagVector tags_[DIMS];
    AdjPtr adjs_[DIMS][DIMS];
};
