
#ifndef GUARD_edge
#define GUARD_edge

#include <vector>

class Edge;
enum VertexType {vrow, vcol};

struct Vertex 
{
    int index;
    VertexType vtype;
    bool visited;
    std::vector<Edge*> edges;
    Vertex(int idx, VertexType a) : 
        index(idx), vtype(a), visited(false) {}
};

class Edge 
{
public:
    Edge(Vertex* const tail, Vertex* const head, std::vector<Edge*>* edge_list, double const weight, int const fixed, double const p, double const lambda);
    int fixed() const { return fixed_; }
    double weight() const {return weight_; }
    double p() const {return p_; }
    double lambda() const {return lambda_; }
    void weight(double w);
    void setPos(int const pos, VertexType vtype) { (vtype == vrow) ? tail_pos_ = pos : head_pos_ = pos; };
    void setEdgeListPos(int const edge_list_pos) { edge_list_pos_ = edge_list_pos; }
    Vertex* ends(VertexType vtype) const { return (vtype == vrow) ? tail_ : head_; };
private:
    void add();
    void remove();
    Vertex * const tail_, * const head_;
    std::vector<Edge*>* const edge_list_;
    double weight_, p_, lambda_;
    int head_pos_, tail_pos_, edge_list_pos_;
    int const fixed_;
};

#endif