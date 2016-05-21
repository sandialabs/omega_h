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
    Read<T> array() const;
    void set_array(Read<T> array);
  private:
    Read<T> array_;
};

template <typename T>
bool is(TagBase const* t);

template <typename T>
Tag<T> const* to(TagBase const* t);
template <typename T>
Tag<T>* to(TagBase* t);
