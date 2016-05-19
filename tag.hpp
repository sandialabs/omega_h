class TagBase {
  public:
    TagBase(std::string const& name, I8 ncomps);
    virtual ~TagBase();
    std::string const& name() const;
    I8 ncomps() const;
    virtual TagBase* clone() const = 0;
  private:
    std::string name_;
    I8 ncomps_;
};

template <typename T>
bool is(TagBase const* t);

template <typename T>
bool to(TagBase const* t);

template <typename T>
class Tag : public TagBase {
  public:
    Tag(std::string const& name, I8 ncomps, Read<T> data);
    Read<T> data() const;
    virtual TagBase* clone() const override;
  private:
    Read<T> data_;
};
