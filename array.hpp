template <class T>
bool operator==(Read<T> a, Read<T> b);

template <class T>
Write<T> deep_copy(Read<T> a);

template <typename T>
typename StandinTraits<T>::type sum(Read<T> a);
template <typename T>
T min(Read<T> a);
template <typename T>
T max(Read<T> a);

template <typename T>
typename StandinTraits<T>::type sum(CommPtr comm, Read<T> a);
template <typename T>
T min(CommPtr comm, Read<T> a);
template <typename T>
T max(CommPtr comm, Read<T> a);

bool are_close(Reals a, Reals b, Real tol = EPSILON, Real floor = EPSILON);

template <typename T>
std::ostream& operator<<(std::ostream& o, Read<T> a);

template <typename T>
Read<T> multiply_each_by(T factor, Read<T> a);
template <typename T>
Read<T> multiply_each(Read<T> a, Read<T> b);
template <typename T>
Read<T> divide_each(Read<T> a, Read<T> b);
template <typename T>
Read<T> add_each(Read<T> a, Read<T> b);
template <typename T>
Read<T> subtract_each(Read<T> a, Read<T> b);
template <typename T>
Read<T> add_to_each(Read<T> a, T b);
template <typename T>
Read<I8> each_geq_to(Read<T> a, T b);
template <typename T>
Read<I8> each_gt(Read<T> a, T b);
template <typename T>
Read<I8> each_lt(Read<T> a, T b);
template <typename T>
Read<I8> gt_each(Read<T> a, Read<T> b);
template <typename T>
Read<I8> each_neq_to(Read<T> a, T b);
template <typename T>
Read<I8> each_eq_to(Read<T> a, T b);
Read<I8> land_each(Read<I8> a, Read<I8> b);
Read<I8> lor_each(Read<I8> a, Read<I8> b);
