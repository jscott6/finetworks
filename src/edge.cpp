
#include "edge.h"
#include <algorithm>

using namespace std;

double const EPS = 1e-9;

Edge::Edge(Vertex* const tail, Vertex* const head, vector<Edge*>* edge_list, double const weight, int const fixed, double const p, double const lambda): 
    tail_(tail), 
    head_(head),
    edge_list_(edge_list),
    weight_(weight), 
    p_(p),
    lambda_(lambda),
    fixed_(fixed)
{
    // add to vertex structure if positive weight and free
    if(weight > EPS && !fixed)  
        add();
}

void Edge::add() 
{
    // add to head vertex
    setPos(head_->edges.size(), vcol);
    head_->edges.push_back(this);
    // add to tail vertex
    setPos(tail_->edges.size(), vrow);
    tail_->edges.push_back(this);
    // add to edge list
    setEdgeListPos(edge_list_->size());
    edge_list_->push_back(this);
}

void Edge::remove() 
{
    // remove from head datastructure
    head_->edges.back()->setPos(head_pos_, vcol);
    swap(head_->edges[head_pos_], head_->edges.back());
    head_->edges.pop_back();
    // remove from tail datastructure
    tail_->edges.back()->setPos(tail_pos_, vrow);
    swap(tail_->edges[tail_pos_], tail_->edges.back());
    tail_->edges.pop_back();
    // remove from edge list
    edge_list_->back()->setEdgeListPos(edge_list_pos_);
    swap((*edge_list_)[edge_list_pos_], edge_list_->back());
    edge_list_->pop_back();
}

void Edge::weight(double w)
{
  if (weight_ < EPS && w > EPS) add();
  if (weight_ > EPS && w < EPS) remove();
  weight_ = w;
}
