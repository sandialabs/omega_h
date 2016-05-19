class TagBase {
  public:
    TagBase(std::string const& name, I8 ncomps);
    virtual ~TagBase();
    std::string const& name() const;
    I8 ncomps() const;
  private:
    std::string name_;
    I8 ncomps_;
};

template <typename T>
class Tag : public TagBase {
  public:
    Tag(std::string const& name, I8 ncomps, Read<T> data);
    Read<T> data() const;
  private:
    Read<T> data_;
};

template <typename T>
bool is(TagBase const* t);

template <typename T>
Tag<T> const* to(TagBase const* t);
