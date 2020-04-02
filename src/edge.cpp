#include "edge.h"
#include <algorithm>

using namespace std;

Edge::Edge(Vertex* const row, Vertex* const col, vector<Edge*>* edge_list, double const weight, int const fixed, double const p, double const lambda, double const tol): 
    row_(row), 
    col_(col),
    edge_list_(edge_list),
    weight_(weight), 
    p_(p),
    lambda_(lambda),
    tolerance_(tol),
    fixed_(fixed)
{
    // add to vertex structure if positive weight and free
    if(weight > tol && !fixed)  
        add();
}

void Edge::add() 
{
    // add to col vertex
    setPos(col_->edges.size(), vcol);
    col_->edges.push_back(this);
    // add to row vertex
    setPos(row_->edges.size(), vrow);
    row_->edges.push_back(this);
    // add to edge list
    setEdgeListPos(edge_list_->size());
    edge_list_->push_back(this);
}

void Edge::remove() 
{
    // remove from col datastructure
    col_->edges.back()->setPos(col_pos_, vcol);
    swap(col_->edges[col_pos_], col_->edges.back());
    col_->edges.pop_back();
    // remove from row datastructure
    row_->edges.back()->setPos(row_pos_, vrow);
    swap(row_->edges[row_pos_], row_->edges.back());
    row_->edges.pop_back();
    // remove from edge list
    edge_list_->back()->setEdgeListPos(edge_list_pos_);
    swap((*edge_list_)[edge_list_pos_], edge_list_->back());
    edge_list_->pop_back();
}

void Edge::weight(double w)
{
  if (weight_ < tolerance_ && w > tolerance_) add();
  if (weight_ > tolerance_ && w < tolerance_) remove();
  weight_ = w;
}
