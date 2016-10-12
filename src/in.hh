#ifndef snip_in_hh
#define snip_in_hh

template <typename T1, typename T2>
inline bool in(T1 x, T2 a) noexcept(noexcept(x==a)) { return x==a; }

template <typename T1, typename T2, typename... TT>
inline bool in(T1 x, T2 a, TT... bb) { return in(x,a) || in(x,bb...); }

#endif
