
#include "edge.h"
#include <algorithm>

using namespace std;

double const EPS = 1e-9;

Edge::Edge(Vertex* const head, Vertex* const tail, double const weight, int const fixed, double const p, double const lambda):
    head_(head), 
    tail_(tail), 
    weight_(weight), 
    p_(p),
    lambda_(lambda),
    fixed_(fixed)
{
    // add to vertex structure if positive weight and free
    if(weight > EPS && !fixed)
    {
        add();
        head_->pos++;
        tail_->pos++;
    }
}

void Edge::add() 
{
    // add to head vertex
    setHeadPos(head_->edges.size());
    head_->edges.push_back(this);
    // add to tail vertex
    setTailPos(tail_->edges.size());
    tail_->edges.push_back(this);
}

void Edge::remove() 
{
    // remove from head datastructure
    head_->edges.back()->head_pos_ = head_pos_;
    swap(head_->edges[head_pos_], head_->edges.back());
    head_->edges.pop_back();
    // remove from tail datastructure
    tail_->edges.back()->tail_pos_ = tail_pos_;
    swap(tail_->edges[tail_pos_], tail_->edges.back());
    tail_->edges.pop_back();
}

void Edge::weight(double w)
{
  if (weight_ < EPS && w > EPS) add();
  if (weight_ > EPS && w < EPS) remove();
  weight_ = w;
}

/*
method moves the edge to final position
in each of the vertex structures
*/

void Edge::moveBack() 
{
    // move to back of head ds
    head_->edges[head_->pos]->setHeadPos(head_pos_); 
    swap(head_->edges[head_pos_], head_->edges[head_->pos]);
    setHeadPos(head_->pos);
    head_->pos--; 
    // move to back of tail ds
    tail_->edges[tail_->pos]->setTailPos(tail_pos_); 
    swap(tail_->edges[tail_pos_], tail_->edges[tail_->pos]);
    setTailPos(tail_->pos);
    tail_->pos--; 
}


