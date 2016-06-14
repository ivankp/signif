#ifndef binner_h
#define binner_h

#include <limits>
#include <type_traits>

/*
 * Same convention as in ROOT TH1:
 * bin = 0;       underflow bin
 * bin = 1;       first bin with low-edge xlow INCLUDED
 * bin = nbins;   last bin with upper-edge xup EXCLUDED
 * bin = nbins+1; overflow bin
 */

template <typename Bin, typename Edge = double>
class binner {
public:
  typedef Bin   bin_t;
  typedef Edge edge_t;
  typedef typename std::vector<edge_t>::iterator  edge_iter;
  typedef typename std::vector< bin_t>::iterator   bin_iter;
  typedef typename std::vector<edge_t>::size_type size_type;

protected:
  std::vector<edge_t> _edges;
  std::vector< bin_t> _bins;

  template <bool Cond, class Type = void>
  using enable_if_t = typename std::enable_if<Cond,Type>::type;

  template <typename T1, typename T2>
  inline enable_if_t<std::is_arithmetic<T1>::value> _fill(T1& x, T2 y) noexcept
  { x += y; }

  template <typename T1, typename... TT>
  inline enable_if_t< !std::is_arithmetic<T1>::value || (sizeof...(TT) > 1) >
  _fill(T1& x, TT... args) noexcept(noexcept(x(args...)))
  { x(args...); }

public:
  binner(size_type nbins, edge_t xlow, edge_t xup)
  : _edges(nbins+1), _bins(nbins+2)
  {
    const edge_t step = (xup-xlow)/nbins;
    for (size_type i=0; i<=nbins; ++i)
      _edges[i] = xlow + i*step;
  }

  binner(std::vector<edge_t>&& edges)
  : _edges(std::move(edges)), _bins(_edges.size()+1)
  { }

  template <typename InputIterator>
  binner(InputIterator first, InputIterator last)
  : _edges(first,last), _bins(_edges.size()+1)
  { }

  size_type fill(edge_t e) {
    size_type i = _edges.size()-1;
    for (;;--i) {
      if (e >= _edges[i]) {
        ++i;
        break;
      }
      if (i==0) break;
    }
    _bins[i]++;
    return i;
  }

  template <typename... TT>
  size_type fill(edge_t e, const TT&... xx) {
    size_type i = _edges.size()-1;
    for (;;--i) {
      if (e >= _edges[i]) {
        ++i;
        break;
      }
      if (i==0) break;
    }
    _fill(_bins[i], xx...);
    return i;
  }

  template <typename... TT>
  void fill_bin(size_type i, const TT&... xx) {
    _fill(_bins.at(i), xx...);
  }

  inline bin_t& operator[](size_type i) noexcept {
    return _bins[i];
  }
  inline const bin_t& operator[](size_type i) const noexcept {
    return _bins[i];
  }

  inline bin_t& bin(size_type i) {
    return _bins.at(i);
  }
  inline const bin_t& bin(size_type i) const {
    return _bins.at(i);
  }

  edge_t ledge(size_type i) const {
    return i ? _edges.at(i-1) : -std::numeric_limits<edge_t>::infinity();
  }
  edge_t redge(size_type i) const {
    return i<_edges.size() ? _edges.at(i)
                           : std::numeric_limits<edge_t>::infinity();
  }

  size_type nbins() const { return _edges.size()-1; }

  const std::vector<edge_t>& edges() const { return _edges; }
  const std::vector< bin_t>&  bins() const { return _bins;  }
};

#endif
