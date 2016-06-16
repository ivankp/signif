#ifndef binner_hh
#define binner_hh

#include <limits>
#include <utility>
#include <algorithm>
#include <type_traits>

/*
 * Same convention as in ROOT TH1:
 * bin = 0;       underflow bin
 * bin = 1;       first bin with low-edge xlow INCLUDED
 * bin = nbins;   last bin with upper-edge xup EXCLUDED
 * bin = nbins+1; overflow bin
 */

template <bool Cond, class Type = void>
using enable_t = typename std::enable_if<Cond,Type>::type;

template <typename Bin, typename Enable = void>
struct binner_default_filler {
  template <typename... TT>
  // inline enable_t<(sizeof...(TT) > 1)>
  inline void operator()(Bin& bin, TT&&... args)
  noexcept(noexcept(bin(std::forward<TT>(args)...)))
  { bin(std::forward<TT>(args)...); }
};

template <typename Bin>
struct binner_default_filler<Bin, enable_t<std::is_arithmetic<Bin>::value>> {
  template <typename T>
  inline void operator()(Bin& bin, T x) noexcept { bin += x; }
  inline void operator()(Bin& bin) noexcept { ++bin; }
};

// TODO: overload for += operator with 1 argument
// http://stackoverflow.com/questions/36296301/c-sfinae-detect-increment-operator

//===================================================================

template <typename Bin, typename Edge = double,
          typename Filler = binner_default_filler<Bin>>
class binner {
public:
  typedef Bin    bin_t;
  typedef Edge   edge_t;
  typedef Filler filler_t;
  typedef typename std::vector<edge_t>::iterator  edge_iter;
  typedef typename std::vector< bin_t>::iterator   bin_iter;
  // typedef typename std::vector<edge_t>::size_type size_type;
  typedef unsigned size_type;

protected:
  std::vector<edge_t> _edges;
  std::vector< bin_t> _bins;

public:
  binner() { }

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

  binner(const std::vector<edge_t>& edges)
  : _edges(edges), _bins(_edges.size()+1)
  { }

  template <typename InputIterator>
  binner(InputIterator first, InputIterator last)
  : _edges(first,last), _bins(_edges.size()+1)
  { }

  template <typename InputIterator>
  void init(InputIterator first, InputIterator last)
  {
    _edges = {first,last};
    _bins = decltype(_bins)(_edges.size()+1);
  }

  //---------------------------------------------

  size_type find_bin(edge_t e) noexcept {
    size_type i = _edges.size()-1;
    for (;;--i) {
      if (e >= _edges[i]) {
        ++i;
        break;
      }
      if (i==0) break;
    }
    return i;
  }

  template <typename... TT>
  size_type fill(edge_t e, TT&&... args) {
    size_type i = find_bin(e);
    filler_t()(_bins[i], std::forward<TT>(args)...);
    return i;
  }

  //---------------------------------------------

  template <typename... TT>
  void fill_bin(size_type i) { ++_bins.at(i); }

  template <typename... TT>
  void fill_bin(size_type i, TT&&... args) {
    filler_t()(_bins.at(i), std::forward<TT>(args)...);
  }

  //---------------------------------------------

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

  //---------------------------------------------

  edge_t ledge(size_type i) const {
    return i ? _edges.at(i-1) : -std::numeric_limits<edge_t>::infinity();
  }
  edge_t redge(size_type i) const {
    return i<_edges.size() ? _edges.at(i)
                           : std::numeric_limits<edge_t>::infinity();
  }

  //---------------------------------------------

  size_type nbins() const { return _edges.size()-1; }

  const std::vector<edge_t>& edges() const { return _edges; }
  const std::vector< bin_t>&  bins() const { return _bins;  }
  std::vector<edge_t>& edges() { return _edges; }
  std::vector< bin_t>&  bins() { return _bins;  }
};

#endif
