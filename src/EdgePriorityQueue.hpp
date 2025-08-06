#ifndef _EDGE_PRIORITY_QUEUE_HPP
#define _EDGE_PRIORITY_QUEUE_HPP

#include "bucket_queue.hpp"

class EdgePriorityQueue : public BucketQueue<Edge *>
{
public:
    EdgePriorityQueue();
    void push_back(Edge *);
};

#endif
