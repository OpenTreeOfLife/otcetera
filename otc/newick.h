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
		NewickTokenizer(std::istream &inp, const std::string & filepath)
			:inputStream(inp),
			inpFilepath(filepath) {
		}
		class iterator;
		class Token {
			public:
				const std::string & content() const {
					return this->tokenContent;
				}
			private:
				Token(const std::string &content, const std::string & fn, std::size_t p)
					:tokenContent(content), 
					inpFilepath(fn),
					pos(p) {
				}
				const std::string tokenContent;
				const std::string inpFilepath;
				const std::size_t pos;
				friend class NewickTokenizer::iterator;
		};

		class iterator : std::forward_iterator_tag {
			public:
				bool operator==(const iterator & other) const {
					LOG(TRACE) << "Equality test atEnd = " << this->atEnd << " other.atEnd = " << other.atEnd << '\n';
					if (other.atEnd) {
						return this->atEnd;
					}
					if (this->atEnd) {
						return false;
					}
					return (&(this->inputStream) == &(other.inputStream)) && (this->pos == other.pos);
				}
				bool operator!=(const iterator & other) const {
						LOG(TRACE) << "Inequality test atEnd = " << this->atEnd << " other.atEnd = " << other.atEnd << '\n';
						return !(*this == other);
					}
				Token operator*() const {
					LOG(TRACE) << "* operator\n";
					return Token(currWord, inpFilepath, pos);
				}
				iterator & operator++() {
					LOG(TRACE) << "increment\n";
					if (this->atEnd) {
						throw std::out_of_range("Incremented a dead NewickTokenizer::iterator");
					}
					char c;
					this->inputStream.get(c);
					if (this->inputStream.eof()) {
						this->atEnd = true;
					} else {
						currWord.assign(1, c);
					}
					return *this;
				}
			private:
				iterator(std::istream &inp, const std::string & filepath)
					:inputStream(inp),
					inpFilepath(filepath),
					atEnd(!inp.good()),
					pos(0) {
					LOG(TRACE) << "create live\n";
					++(*this);
				}
				iterator(std::istream &inp) // USE in end() ONLY!
					:inputStream(inp),
					atEnd(true),
					pos(0) {
					LOG(TRACE) << "create dead\n";
					
				}
				std::istream & inputStream;
				const std::string inpFilepath;
				bool atEnd;
				size_t pos;
				std::string currWord;
				friend class NewickTokenizer;
		};
		iterator begin() {
			iterator b(this->inputStream, this->inpFilepath);
			return b;
		}
		iterator end() {
			return iterator(this->inputStream);
		}

	private:
		std::istream & inputStream;
		std::string inpFilepath;
};

//Takes wide istream and (optional) filepath (just used for error reporting if not empty)
template<typename T>
std::unique_ptr<RootedTree<T> > readNextNewick(std::istream &inp, const std::string & filepath);

template<typename T>
inline std::unique_ptr<RootedTree<T> > readNextNewick(std::istream &inp, const std::string & filepath) {
	assert(inp.good());
	NewickTokenizer tokenizer(inp, filepath);
	auto i = 0U;
	for (auto token : tokenizer) {
		std::cout << "token = \"" << token.content() << "\"\n"; 
		if (i > 1000U) {
			break;
		}
	}
	return std::unique_ptr<RootedTree<T> > (new RootedTree<T>());
}



}// namespace otc
#endif
