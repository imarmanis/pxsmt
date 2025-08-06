#ifndef _BUCKET_QUEUE_HPP
#define _BUCKET_QUEUE_HPP

#include "utility"
#include "vector"
#include "graph.hpp"

template <typename T>
class BucketQueue
{
private:
    std::vector<std::vector<T>> buckets;
    unsigned size;
    unsigned cost;
public:
    BucketQueue(unsigned d) : buckets(d), size(0), cost(0) { }
    bool empty() { return size == 0; }
    void insert(unsigned b, T elem)
    { 
        buckets[(cost + b) % buckets.size()].push_back(elem);
        size++;
    }
    T back()
    {
        unsigned d = buckets.size();
        while(buckets[cost % d].empty()) cost++;
        return buckets[cost % d].back();
    }
    void pop_back()
    {
        unsigned d = buckets.size();
        while(buckets[cost % d].empty()) cost++;
        buckets[cost % d].pop_back();
        size--;
    }
    unsigned getCost()
    {
        return cost;
    }
};


#endif
