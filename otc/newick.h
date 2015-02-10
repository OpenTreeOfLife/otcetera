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

class NewickTokenizer {
	public:
		NewickTokenizer(std::wistream &inp, const std::string & filepath)
			:inputStream(inp),
			inpFilepath(filepath) {
		}
		class iterator;
		class Token {
			public:
				const std::wstring & content() const {
					return this->tokenContent;
				}
			private:
				Token(const std::wstring &content, const std::string & fn, std::size_t p)
					:tokenContent(content), 
					inpFilepath(fn),
					pos(p) {
				}
				const std::wstring tokenContent;
				const std::string inpFilepath;
				const std::size_t pos;
				friend class NewickTokenizer::iterator;
		};

		class iterator : std::forward_iterator_tag {
			public:
				bool operator==(const iterator & other) const {
					if (other.atEnd) {
						return this->atEnd;
					}
					if (this->atEnd) {
						return false;
					}
					return (&(this->inputStream) == &(other.inputStream)) && (this->pos == other.pos);
				}
				bool operator!=(const iterator & other) const {
						return !(*this == other);
					}
				Token operator*() const {
					return Token(currWord, inpFilepath, pos);
				}
				iterator & operator++() {
					if (this->atEnd) {
						throw std::out_of_range("Incremented a dead NewickTokenizer::iterator");
					}
					wchar_t c;
					this->inputStream.get(c);
					if (this->inputStream.eof()) {
						this->atEnd = true;
					} else {
						currWord.assign(1, c);
					}
				}
			private:
				iterator(std::wistream &inp, const std::string & filepath)
					:inputStream(inp),
					inpFilepath(filepath),
					atEnd(inp.good()),
					pos(0) {
				}
				iterator(std::wistream &inp) // USE in end() ONLY!
					:inputStream(inp),
					atEnd(true),
					pos(0) {
				}
				std::wistream & inputStream;
				const std::string inpFilepath;
				bool atEnd;
				size_t pos;
				std::wstring currWord;
				friend class NewickTokenizer;
		};
		iterator begin() {
			return iterator(this->inputStream, this->inpFilepath);
		}
		iterator end() {
			return iterator(this->inputStream, this->inpFilepath);
		}

	private:
		std::wistream & inputStream;
		std::string inpFilepath;
};

//Takes wide istream and (optional) filepath (just used for error reporting if not empty)
template<typename T>
std::unique_ptr<RootedTree<T> > readNextWNewick(std::wistream &inp, const std::string & filepath);

template<typename T>
inline std::unique_ptr<RootedTree<T> > readNextWNewick(std::wistream &inp, const std::string & filepath) {
	NewickTokenizer tokenizer(inp, filepath);
	for (auto token : tokenizer) {
		std::cout << "token = \"" << token->content() << "\"\n"; 
	}
	return std::unique_ptr<RootedTree<T> > (new RootedTree<T>());
}



}// namespace otc
#endif
