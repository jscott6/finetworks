
#ifndef GUARD_edge
#define GUARD_edge

#include <vector>

class Edge;

struct Vertex 
{
    int index;
    bool visited;
    std::vector<Edge*> edges;
    Vertex() : visited(false) {}
};

class Edge 
{
public:
    Edge(Vertex* const head, Vertex* const tail, std::vector<Edge*>* edge_list, double const weight, int const fixed, double const p, double const lambda);
    int fixed() const { return fixed_; }
    double weight() const {return weight_; }
    double p() const {return p_; }
    double lambda() const {return lambda_; }
    void weight(double w);
    void setHeadPos(int const head_pos) { head_pos_ = head_pos; }
    void setTailPos(int const tail_pos) { tail_pos_ = tail_pos; }
    void setEdgeListPos(int const edge_list_pos) { edge_list_pos_ = edge_list_pos; }
    Vertex* head() const { return head_; }
    Vertex* tail() const { return tail_; }
private:
    void add();
    void remove();
    Vertex* const head_, * const tail_;
    std::vector<Edge*>* const edge_list_;
    double weight_, p_, lambda_;
    int head_pos_, tail_pos_, edge_list_pos_;
    int const fixed_;
};

#endif