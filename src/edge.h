
#ifndef GUARD_edge
#define GUARD_edge

#include <vector>

class Edge;

struct Vertex 
{
    int index;
    int pos;
    std::vector<Edge*> edges;
};

class Edge 
{
public:
    Edge(Vertex* const head, Vertex* const tail, double const weight, int const fixed);
    int fixed() const { return fixed_; }
    double weight() const {return weight_; }
    void setHeadPos(int const head_pos) { head_pos_ = head_pos; }
    void setTailPos(int const tail_pos) { tail_pos_ = tail_pos; }
    int head() const { return head_->index; }
    int tail() const { return tail_->index; }
    void moveBack();
private:
    void add();
    void remove();
    Vertex* const head_, * const tail_;
    double weight_;
    int head_pos_, tail_pos_;
    int const fixed_;
};

#endif