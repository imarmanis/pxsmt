#include "EdgePriorityQueue.hpp"
#include "graph.hpp"

EdgePriorityQueue::EdgePriorityQueue() : BucketQueue<Edge *>(MAX_EDGE_COST + 1) {}

void EdgePriorityQueue::push_back(Edge *e) { insert(e->getCost(), e); }
