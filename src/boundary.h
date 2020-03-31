#ifndef GUARD_boundary
#define GUARD_boundary

#include "edge.h"

struct Boundary
{
    double dlow, dup;
    unsigned int nlow, nup;
    Edge *elow, *eup;
    Boundary(Edge* e, Edge* f)
      : dlow(-e->weight()), 
        dup(f->weight()),
        nlow(1), 
        nup(1),
        elow(e),
        eup(f)
  {
  }
};

#endif