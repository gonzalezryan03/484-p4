/*
 * TODO: Implement partition, probe functions in Join.cpp
 */
#ifndef _JOIN_HPP_
#define _JOIN_HPP_

#include "Bucket.hpp"
#include "Mem.hpp"
#include <cstdint> // For uint32_t, if needed
#include <utility> // For std::pair
#include <vector>

// Define 'uint' if it's not defined elsewhere
typedef unsigned int uint;

/*
 * partition function
 *
 * Input:
 * disk: pointer of Disk object
 * mem: pointer of Memory object
 * left_rel: [left_rel.first, left_rel.second) range of page ids for the left relation
 * right_rel: [right_rel.first, right_rel.second) range of page ids for the right relation
 *
 * Output:
 * A vector of buckets of size (MEM_SIZE_IN_PAGE - 1).
 * Each bucket represents a partition of both relations.
 * See Bucket class for more information.
 */
std::vector<Bucket> partition(Disk* disk, Mem* mem, std::pair<uint, uint> left_rel, std::pair<uint, uint> right_rel);

/*
 * probe function
 * Input:
 * disk: pointer of Disk object
 * mem: pointer of Memory object
 * partitions: a reference to a vector of buckets from the partition function
 *
 * Output:
 * A vector of page ids that contains the join result.
 */
std::vector<uint> probe(Disk* disk, Mem* mem, std::vector<Bucket>& partitions);

#endif // _JOIN_HPP_
