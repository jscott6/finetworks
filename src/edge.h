
#ifndef GUARD_edge
#define GUARD_edge

#include <vector>

struct Vertex 
{
    int const index;
    int pos;
    std::vector<Edge*> edges;
};

class Edge 
{
public:
    Edge(Vertex* const head, Vertex* const tail, int const weight, int const fixed);
    int fixed() const { return fixed_; }
    float weight() const {return weight_; }
    void setHeadPos(int const head_pos) { head_pos_ = head_pos; }
    void setTailPos(int const tail_pos) { tail_pos_ = tail_pos; }
    void moveBack();
private:
    void add();
    void remove();
    Vertex* const head_, * const tail_;
    float weight_;
    int head_pos_, tail_pos_;
    int const fixed_;
};

#endif