#ifndef OTCETERA_NEWICK_H
#define OTCETERA_NEWICK_H
#include <iostream>
#include <fstream>
#include "otc/otcetera.h"
#include "otc/tree.h"

namespace otc {

/* TREE IO */
#if 0
// non wide string versions
template<typename T>
std::unique_ptr<RootedTree<T> > readNextNewick(std::istream &inp);
template<typename T>
unsigned int readNewickStream(std::istream &inp, bool (*callback)(std::unique_ptr<RootedTree<T> >));

template<typename T>
inline std::unique_ptr<RootedTree<T> > readNextNewick(std::istream &inp) {
	return std::unique_ptr<RootedTree<T> > (new RootedTree<T>());
}

template<typename T>
std::unique_ptr<RootedTree<T> > readNextWNewick(std::wistream &inp);
template<typename T>
inline unsigned int readNewickStream(std::istream &inp, bool (*callback)(std::unique_ptr<RootedTree<T> >)) {
	auto c = 0U;
	for (;;) {
		std::unique_ptr<RootedTree<T> > nt = readNextNewick<T>(inp);
		if (nt == nullptr) {
			return c;
		}
		++c;
		auto cbr = callback(std::move(nt));
		if (!cbr) {
			return c;
		}
	}
}
#endif

template<typename T>
inline std::unique_ptr<RootedTree<T> > readNextWNewick(std::wistream &inp) {
	return std::unique_ptr<RootedTree<T> > (new RootedTree<T>());
}



}// namespace otc
#endif
