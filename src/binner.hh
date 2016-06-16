#ifndef snip_binner_hh
#define snip_binner_hh

#include <limits>
#include <utility>

#include "type_traits_extra.hh"

/*
 * Same convention as in ROOT TH1:
 * bin = 0;       underflow bin
 * bin = 1;       first bin with low-edge xlow INCLUDED
 * bin = nbins;   last bin with upper-edge xup EXCLUDED
 * bin = nbins+1; overflow bin
 */

template <typename Bin>
struct binner_filler_plus_eq {
  template <typename T>
  inline void operator()(Bin& bin, T&& x)
  noexcept(noexcept(bin += std::forward<T>(x)))
  { bin += std::forward<T>(x); }
};

template <typename Bin, typename Enable = void>
struct binner_filler_default: public binner_filler_plus_eq<Bin> {
  template <typename... TT>
  inline enable_t< sizeof...(TT)!=1 ||
                  !has_op_plus_eq<Bin,type_at<0,TT...>>::value >
  operator()(Bin& bin, TT&&... args)
  noexcept(noexcept(bin(std::forward<TT>(args)...)))
  { bin(std::forward<TT>(args)...); }

  using binner_filler_plus_eq<Bin>::operator();
};

template <typename Bin>
struct binner_filler_default<Bin, enable_t<
  has_op_pre_increment<Bin>::value
>>: public binner_filler_plus_eq<Bin> {
  inline void operator()(Bin& bin) noexcept(noexcept(++bin)) { ++bin; }

  using binner_filler_plus_eq<Bin>::operator();
};

template <typename Bin>
struct binner_filler_default<Bin, enable_t<
  !has_op_pre_increment<Bin>::value && has_op_post_increment<Bin>::value
>>: public binner_filler_plus_eq<Bin> {
  inline void operator()(Bin& bin) noexcept(noexcept(bin++)) { bin++; }

  using binner_filler_plus_eq<Bin>::operator();
};

//===================================================================

template <typename Bin, typename Edge = double,
          typename Filler = binner_filler_default<Bin>>
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
  binner(const binner& o): _edges(o._edges), _bins(o._bins) { }
  binner(binner&& o) noexcept
  : _edges(std::move(o._edges)), _bins(std::move(o._bins)) { }

  binner& operator=(const binner& o) {
    _edges = o._edges;
    _bins = o._bins;
    return *this;
  }
  binner& operator=(binner&& o) noexcept {
    _edges = std::move(o._edges);
    _bins = std::move(o._bins);
    return *this;
  }

  binner(size_type nbins, edge_t xlow, edge_t xup)
  : _edges(nbins+1), _bins(nbins+2)
  {
    const edge_t step = (xup-xlow)/nbins;
    for (size_type i=0; i<=nbins; ++i)
      _edges[i] = xlow + i*step;
  }

  binner(const std::vector<edge_t>& edges)
  : _edges(edges), _bins(_edges.size()+1)
  { }
  binner& operator=(const std::vector<edge_t>& edges) {
    _edges = edges;
    _bins = std::vector<bin_t>(_edges.size()+1);
    return *this;
  }

  binner(std::vector<edge_t>&& edges) noexcept
  : _edges(std::move(edges)), _bins(_edges.size()+1)
  { }
  binner& operator=(std::vector<edge_t>&& edges) noexcept {
    _edges = std::move(edges);
    _bins = std::vector<bin_t>(_edges.size()+1);
    return *this;
  }

  binner(std::initializer_list<edge_t> il)
  : _edges(il), _bins(_edges.size()+1)
  { }
  binner& operator=(std::initializer_list<edge_t> il) {
    _edges = il;
    _bins = std::vector<bin_t>(_edges.size()+1);
    return *this;
  }

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
  size_type fill(edge_t e, TT&&... args)
  noexcept(noexcept( filler_t()(std::declval<bin_t&>(),
                                std::forward<TT>(args)...) ))
  {
    size_type i = find_bin(e);
    filler_t()(_bins[i], std::forward<TT>(args)...);
    return i;
  }

  template <typename... TT>
  inline size_type operator()(edge_t e, TT&&... args)
  noexcept(noexcept(
    std::declval<binner&>().fill(e, std::forward<TT>(args)...) ))
  {
    return fill(e, std::forward<TT>(args)...);
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

  inline size_type nbins() const { return _bins.size()-2; }

  inline const std::vector<edge_t>& edges() const noexcept { return _edges; }
  inline const std::vector< bin_t>&  bins() const noexcept { return _bins;  }
  inline std::vector<edge_t>& edges() noexcept { return _edges; }
  inline std::vector< bin_t>&  bins() noexcept { return _bins;  }
};

#endif
