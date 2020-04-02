#ifndef GUARD_edge
#define GUARD_edge

#include <vector>
#include <RcppCommon.h>
#include "vertex.h"

class Edge 
{
public:
    Edge(Vertex* const row, Vertex* const col, std::vector<Edge*>* edge_list, double const weight, int const fixed, double const p, double const lambda, double const tol = 1e-12);
    int fixed() const { return fixed_; }
    double weight() const {return weight_; }
    double p() const {return p_; }
    double lambda() const {return lambda_; }
    void weight(double w);
    void setPos(int const pos, VertexType vtype) { (vtype == vrow) ? row_pos_ = pos : col_pos_ = pos; };
    void setEdgeListPos(int const edge_list_pos) { edge_list_pos_ = edge_list_pos; }
    Vertex* ends(VertexType vtype) const { return (vtype == vrow) ? row_ : col_; };
private:
    void add();
    void remove();
    Vertex * const row_, * const col_;
    std::vector<Edge*>* const edge_list_;
    double weight_, p_, lambda_, tolerance_;
    int col_pos_, row_pos_, edge_list_pos_;
    int const fixed_;
};

#include <Rcpp.h>

#endif