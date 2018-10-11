
#include "edge.h"
#include <algorithm>

using namespace std;

double const EPS = 1e-9;

Edge::Edge(Vertex* const head, Vertex* const tail, double const weight, int const fixed):
    head_(head), 
    tail_(tail), 
    weight_(weight), 
    fixed_(fixed) 
{
    // add to vertex structure if positive weight and free
    if(weight > EPS && !fixed)
        add();
    return;
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


/*
method moves the edge to final position
in each of the vertex structures
*/

void Edge::moveBack() 
{
    // move to back of head ds
    head_->edges[head_->pos]->setHeadPos(head_pos_); 
    setHeadPos(head_->pos);
    swap(head_->edges[head_pos_], head_->edges[head_->pos]);
    head_->pos--; 
    // move to back of tail ds
    tail_->edges[tail_->pos]->setTailPos(tail_pos_); 
    setTailPos(tail_->pos);
    swap(tail_->edges[tail_pos_], tail_->edges[tail_->pos]);
    tail_->pos--; 
}

