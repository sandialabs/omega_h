class Mesh {
  public:
    Mesh();
    void set_dim(Int dim);
    void set_verts(LO nverts);
    void set_ents(Int dim, Adj down);
    Int dim() const;
    LO nents(Int dim) const;
    LO nelems() const;
    LO nverts() const;
    template <typename T>
    void add_tag(Int dim, std::string const& name, Int ncomps);
    template <typename T>
    void set_tag(Int dim, std::string const& name, Read<T> data);
    template <typename T>
    Tag<T> const& get_tag(Int dim, std::string const& name) const;
    void remove_tag(Int dim, std::string const& name);
    bool has_tag(Int dim, std::string const& name) const;
    Int count_tags(Int dim) const;
    TagBase const* get_tag(Int dim, Int i) const;
    bool has_ents(Int dim) const;
    bool has_adj(Int from, Int to) const;
    Adj get_adj(Int from, Int to) const;
    Adj ask_adj(Int from, Int to);
    Read<GO> ask_globals(Int dim);
  private:
    typedef std::shared_ptr<TagBase> TagPtr;
    typedef std::shared_ptr<Adj> AdjPtr;
    typedef std::vector<TagPtr> TagVector;
    typedef TagVector::iterator TagIter;
    typedef TagVector::const_iterator TagCIter;
    TagIter tag_iter(Int dim, std::string const& name);
    TagCIter tag_iter(Int dim, std::string const& name) const;
    void check_dim(Int dim) const;
    void check_dim2(Int dim) const;
    void add_adj(Int from, Int to, Adj adj);
    Adj derive_adj(Int from, Int to);
    Int dim_;
    LO nents_[DIMS];
    TagVector tags_[DIMS];
    AdjPtr adjs_[DIMS][DIMS];
};
