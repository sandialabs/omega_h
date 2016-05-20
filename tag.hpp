class TagBase {
  public:
    TagBase(std::string const& name, Int ncomps);
    virtual ~TagBase();
    std::string const& name() const;
    Int ncomps() const;
  private:
    std::string name_;
    Int ncomps_;
};

template <typename T>
class Tag : public TagBase {
  public:
    Tag(std::string const& name, Int ncomps);
    Read<T> data() const;
    void set_data(Read<T> data);
  private:
    Read<T> data_;
};

template <typename T>
bool is(TagBase const* t);

template <typename T>
Tag<T> const* to(TagBase const* t);
template <typename T>
Tag<T>* to(TagBase* t);
