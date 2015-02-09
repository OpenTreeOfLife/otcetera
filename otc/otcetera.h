#if !defined OTCETERA_H
#define OTCETERA_H
#include <memory>
#include <iostream>

#include "otc/tree.h"
namespace otc {

template<typename T>
std::unique_ptr<RootedTree<T> > readNextNewick(std::iostream &inp);

template<typename T>
inline unsigned int readNewickStream(std::iostream &inp, bool (*callback)(std::unique_ptr<RootedTree<T> >)) {
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





}
#endif
