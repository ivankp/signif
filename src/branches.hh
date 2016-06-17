// Developed by Ivan Pogrebnyak

#ifndef ttree_branches_hh
#define ttree_branches_hh

#include <TTree.h>

template<typename T>
inline void branches_impl(TTree* tree, const char* name, T* ptr) {
  tree->SetBranchAddress(name, ptr);
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus(name,1);
}

template<typename T, typename... TT>
inline void branches_impl(
  TTree* tree, const char* name, T* ptr, TT... bb
) {
  tree->SetBranchAddress(name, ptr);
  branches_impl(tree,bb...);
  tree->SetBranchStatus(name,1);
}

template<typename... TT>
inline void branches(TTree* tree, TT... bb) {
  branches_impl(tree,bb...);
}

template<typename... TT>
inline void branches_reset(TTree* tree) {
  tree->SetBranchStatus("*",1);
}

template<typename... TT>
inline void branches_reset(TTree* tree, bool on) {
  tree->SetBranchStatus("*",on);
}

template<typename... TT>
inline void branches_reset(TTree* tree, TT... bb) {
  branches_reset(tree);
  branches_impl(tree,bb...);
}

//===============================================

template<typename T>
inline void branches_set(TTree* tree, const char* name, T* ptr) {
  tree->SetBranchAddress(name, ptr);
}

template<typename T, typename... TT>
inline void branches_set(TTree* tree, const char* name, T* ptr, TT... bb) {
  branches_set(tree, name, ptr);
  branches_set(tree,bb...);
}

//===============================================

template<typename T>
inline void branches_set_on(TTree* tree, const char* name, T* ptr) {
  tree->SetBranchStatus(name,1);
  tree->SetBranchAddress(name, ptr);
}

template<typename T, typename... TT>
inline void branches_set_on(
  TTree* tree, const char* name, T* ptr, TT... bb
) {
  branches_set_on(tree, name, ptr);
  branches_set_on(tree,bb...);
}


#endif
