
#ifndef GUARD_edge
#define GUARD_edge

#include <vector>

class Edge;

struct Vertex 
{
    int index;
    int pos;
    std::vector<Edge*> edges;
    Vertex() : pos(-1) {}
};

class Edge 
{
public:
    Edge(Vertex* const head, Vertex* const tail, double const weight, int const fixed, double const p, double const lambda);
    int fixed() const { return fixed_; }
    double weight() const {return weight_; }
    double p() const {return p_; }
    double lambda() const {return lambda; }
    void weight(double w);
    void setHeadPos(int const head_pos) { head_pos_ = head_pos; }
    void setTailPos(int const tail_pos) { tail_pos_ = tail_pos; }
    Vertex* head() const { return head_; }
    Vertex* tail() const { return tail_; }
    void moveBack();
private:
    void add();
    void remove();
    Vertex* const head_, * const tail_;
    double weight_, p_, lambda_;
    int head_pos_, tail_pos_;
    int const fixed_;
};

#endif