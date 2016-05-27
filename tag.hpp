enum TagType {
  OSH_I8  = 0,
  OSH_I32 = 2,
  OSH_I64 = 3,
  OSH_F64 = 5,
};

class TagBase {
  public:
    TagBase(std::string const& name, Int ncomps);
    virtual ~TagBase();
    std::string const& name() const;
    Int ncomps() const;
    virtual TagType type() const = 0;
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
    virtual TagType type() const override;
  private:
    Read<T> array_;
};

template <typename T>
bool is(TagBase const* t);

template <typename T>
Tag<T> const* to(TagBase const* t);
template <typename T>
Tag<T>* to(TagBase* t);
