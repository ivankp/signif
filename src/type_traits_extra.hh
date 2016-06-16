#ifndef snip_type_traits_extra_hh
#define snip_type_traits_extra_hh

#include <type_traits>
#include <tuple>

template <bool Cond, class Type = void>
using enable_t = typename std::enable_if<Cond,Type>::type;

// ========================================================

template <size_t I, typename... TT>
using type_at = typename std::tuple_element<I, std::tuple<TT...>>::type;

// operator support =======================================
// http://stackoverflow.com/questions/36296301/c-sfinae-detect-increment-operator

template <typename T>
struct has_op_pre_increment {
  template <typename U, typename = decltype(++std::declval<U&>())>
  static long A(U&&);
  static char A(...);

  static constexpr bool value = sizeof(A(std::declval<T>())) == sizeof(long);
};

template <typename T>
struct has_op_post_increment {
  template <typename U, typename = decltype(std::declval<U&>()++)>
  static long A(U&&);
  static char A(...);

  static constexpr bool value = sizeof(A(std::declval<T>())) == sizeof(long);
};

template <typename T1, typename T2>
struct has_op_plus_eq {
  template <typename U1, typename U2,
            typename = decltype(std::declval<U1&>()+=std::declval<U2>())>
  static long A(U1&&, U2&&);
  static char A(...);

  static constexpr bool value =
    sizeof(A(std::declval<T1>(),std::declval<T2>())) == sizeof(long);
};

template <typename T1, typename... TT>
struct has_op_parenth {
  template <typename U1, typename... UU,
            typename = decltype(std::declval<U1&>()(std::declval<UU>()...))>
  static long A(U1&&, UU&&...);
  static char A(...);

  static constexpr bool value =
    sizeof(A(std::declval<T1>(),std::declval<TT>()...)) == sizeof(long);
};

#endif
