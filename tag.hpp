enum TagType {
  OSH_I8  = 0,
  OSH_I32 = 2,
  OSH_I64 = 3,
  OSH_F64 = 5,
};

enum Xfer {
  OSH_DONT_TRANSFER,
  OSH_INHERIT,
  OSH_LINEAR_INTERP,
  OSH_POINTWISE,
  OSH_CONSERVE,
  OSH_GLOBAL,
  OSH_LENGTH,
  OSH_QUALITY,
  OSH_METRIC
};

class TagBase {
  public:
    TagBase(std::string const& name, Int ncomps, Xfer xfer);
    virtual ~TagBase();
    std::string const& name() const;
    Int ncomps() const;
    Xfer xfer() const;
    virtual TagType type() const = 0;
  private:
    std::string name_;
    Int ncomps_;
    Xfer xfer_;
};

template <typename T>
class Tag : public TagBase {
  public:
    Tag(std::string const& name, Int ncomps, Xfer xfer);
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
