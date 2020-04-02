#ifndef GUARD_vertex
#define GUARD_vertex

#include <vector>

enum VertexType {vrow, vcol};
class Edge;

struct Vertex 
{
    int index;
    VertexType vtype;
    bool visited;
    std::vector<Edge*> edges;
    Vertex(int idx, VertexType a) : 
        index(idx), vtype(a), visited(false) {}
};

#endif