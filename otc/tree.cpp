#include "otc/tree.h"
using namespace otc;


/*
template<typename T> class TreeNode;
template<typename T>
class TreeEdge {
	public:
		typedef class TreeNode<T> Node_t;
		const Node_t * GetParent() const {
			return parent;
		}
		const Node_t * GetChild() const {
			return child;
		}
		T & GetData() const {
			return child->GetData();
		}
	private:
		void SetParent(Node_t * p) {
			this->parent = p;
		}
		void WriteAsNewick(std::ostream &out, bool nhx) const;
		TreeEdge<T>(Node_t * par, Node_t  * des)
			:parent(par),
			child(des) {
		}
		Node_t * parent;
		Node_t * child;
		friend class TreeNode<T>;
};
*/
