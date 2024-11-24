/*
 * TODO: Implement partition, probe functions in Join.cpp
 */
#ifndef _JOIN_HPP_
#define _JOIN_HPP_

#include "Bucket.hpp"
#include "Mem.hpp"

/*
 * partition function
 *
 * Input:
 * disk: pointer of Disk object
 * mem: pointer of Memory object
 * left_rel: [left_rel.first, left_rel.second) will be the range of page ids of left relation to join
 * right_rel: [right_rel.first, right_rel.second) will be the range of page ids of right relation to join
 *
 * Output:
 * A vector of buckets of size (MEM_SIZE_IN_PAGE - 1).
 * Each bucket represent a partition of both relation.
 * See Bucket class for more information.
*/
std::vector<Bucket> partition(Disk* disk, Mem* mem,
                              std::pair<uint, uint> left_rel,
                              std::pair<uint, uint> right_rel){
                                // partitions init
                                vector<Bucket> partitions(MEM_SIZE_IN_PAGE - 1, Bucket(disk));

                                // process left rel
                                for (uint page_id = left_rel.first; page_id < left_rel.second; ++page_id){
                                    // load into mem page 0
                                    mem->loadFromDisk(disk, page_id, 0);
                                    Page* input_page = mem->mem_page(0);
                                    for (uint i = 0; i < input_page->size(); ++i){
                                        Record record = input_page->get_record(i);
                                         // Compute h1
                                        uint h1 = record.partition_hash() % (MEM_SIZE_IN_PAGE - 1);

                                        // Output record into mem->mem_page(h1 + 1)
                                        Page* output_page = mem->mem_page(h1 + 1);
                                        output_page->loadRecord(record);

                                        // If output page is full, flush to disk
                                        if (output_page->full()) {
                                            uint new_disk_page_id = mem->flushToDisk(disk, h1 + 1);
                                            partitions[h1].add_left_rel_page(new_disk_page_id);
                                        }
                                    }
                                }
                                  // Flush any remaining records in output buffers for left relation
                                for (uint h1 = 0; h1 < MEM_SIZE_IN_PAGE - 1; ++h1) {
                                    Page* output_page = mem->mem_page(h1 + 1);
                                    if (!output_page->empty()) {
                                        uint new_disk_page_id = mem->flushToDisk(disk, h1 + 1);
                                        partitions[h1].add_left_rel_page(new_disk_page_id);
                                    }
                                }

                                // Similarly process right relation
                                for (uint page_id = right_rel.first; page_id < right_rel.second; ++page_id) {
                                    // Load disk page into mem page 0
                                    mem->loadFromDisk(disk, page_id, 0);

                                    // Process each record in mem->mem_page(0)
                                    Page* input_page = mem->mem_page(0);
                                    for (uint i = 0; i < input_page->size(); ++i) {
                                        Record record = input_page->get_record(i);

                                        // Compute h1
                                        uint h1 = record.partition_hash() % (MEM_SIZE_IN_PAGE - 1);

                                        // Output record into mem->mem_page(h1 + 1)
                                        Page* output_page = mem->mem_page(h1 + 1);
                                        output_page->loadRecord(record);

                                        // If output page is full, flush to disk
                                        if (output_page->full()) {
                                            uint new_disk_page_id = mem->flushToDisk(disk, h1 + 1);
                                            partitions[h1].add_right_rel_page(new_disk_page_id);
                                        }
                                    }
                                }

                                // Flush any remaining records in output buffers for right relation
                                for (uint h1 = 0; h1 < MEM_SIZE_IN_PAGE - 1; ++h1) {
                                    Page* output_page = mem->mem_page(h1 + 1);
                                    if (!output_page->empty()) {
                                        uint new_disk_page_id = mem->flushToDisk(disk, h1 + 1);
                                        partitions[h1].add_right_rel_page(new_disk_page_id);
                                    }
                                }

                                return partitions;
                                }

/*
 * probe function
 * Input:
 * disk: pointer of Disk object
 * mem: pointer of Memory object
 * partition: a reference to a vector of buckets from partition function
 *
 * Output:
 * A vector of page ids that contains the join result.
*/
std::vector<uint> probe(Disk* disk, Mem* mem, std::vector<Bucket>& partitions){
    vector<uint> disk_pages;
    for (uint i = 0; i < partitions.size(); ++i) {
        Bucket& partition = partitions[i];

        // Skip empty ones
        if (partition.num_left_rel_record == 0 || partition.num_right_rel_record == 0) {
            continue;
        }

        // find the smaller relation
        bool left_smaller = (partition.num_left_rel_record <= partition.num_right_rel_record);

        uint hash_table_size = MEM_SIZE_IN_PAGE - 2; // pages 1 to MEM_SIZE_IN_PAGE - 2 are for hash table
        uint hash_table_capacity = hash_table_size * RECORDS_PER_PAGE;

        if ((left_smaller && partition.num_left_rel_record > hash_table_capacity) ||
            (!left_smaller && partition.num_right_rel_record > hash_table_capacity)) {
            continue;
        }

        // Reset memory
        mem->reset();

        if (left_smaller) {
            // Build hash table from left relation
            // mem->mem_page(1) to mem->mem_page(MEM_SIZE_IN_PAGE - 2) are for hash table
            // mem->mem_page(0) is input buffer
            // mem->mem_page(MEM_SIZE_IN_PAGE - 1) is output buffer

            // Build hash table
            vector<Page*> hash_table(hash_table_size, nullptr);
            for (uint k = 0; k < hash_table_size; ++k) {
                hash_table[k] = mem->mem_page(k + 1);
                hash_table[k]->reset();
            }

            // For each page in left relation
            vector<uint> left_pages = partition.get_left_rel();
            for (uint page_id : left_pages) {
                // Load page into mem->mem_page(0)
                mem->loadFromDisk(disk, page_id, 0);
                Page* input_page = mem->mem_page(0);

                // For each record
                for (uint j = 0; j < input_page->size(); ++j) {
                    Record record = input_page->get_record(j);

                    // Compute h2
                    uint h2 = record.probe_hash() % hash_table_size;

                    // Insert record into hash table at mem->mem_page(h2 + 1)
                    Page* hash_page = hash_table[h2];
                    hash_page->loadRecord(record);
                }
            }

            // Now probe right relation
            vector<uint> right_pages = partition.get_right_rel();
            for (uint page_id : right_pages) {
                // Load page into mem->mem_page(0)
                mem->loadFromDisk(disk, page_id, 0);
                Page* input_page = mem->mem_page(0);

                // For each record
                for (uint j = 0; j < input_page->size(); ++j) {
                    Record record = input_page->get_record(j);

                    // Compute h2
                    uint h2 = record.probe_hash() % hash_table_size;

                    // Get matching records from hash table
                    Page* hash_page = hash_table[h2];

                    // For each record in hash_page
                    for (uint m = 0; m < hash_page->size(); ++m) {
                        Record left_record = hash_page->get_record(m);

                        // Check if records match
                        if (left_record == record) { // overloaded == operator

                            // Output the pair
                            Page* output_page = mem->mem_page(MEM_SIZE_IN_PAGE - 1);
                            output_page->loadPair(left_record, record);

                            // If output page is full, flush to disk
                            if (output_page->full()) {
                                uint new_disk_page_id = mem->flushToDisk(disk, MEM_SIZE_IN_PAGE - 1);
                                disk_pages.push_back(new_disk_page_id);
                            }
                        }
                    }
                }
            }

            // Flush any remaining output records
            Page* output_page = mem->mem_page(MEM_SIZE_IN_PAGE - 1);
            if (!output_page->empty()) {
                uint new_disk_page_id = mem->flushToDisk(disk, MEM_SIZE_IN_PAGE - 1);
                disk_pages.push_back(new_disk_page_id);
            }
        } else {
            // Build hash table from right relation
            // Similar to above, but switch roles

            vector<Page*> hash_table(hash_table_size, nullptr);
            for (uint k = 0; k < hash_table_size; ++k) {
                hash_table[k] = mem->mem_page(k + 1);
                hash_table[k]->reset();
            }

            // For each page in right relation
            vector<uint> right_pages = partition.get_right_rel();
            for (uint page_id : right_pages) {
                // Load page into mem->mem_page(0)
                mem->loadFromDisk(disk, page_id, 0);
                Page* input_page = mem->mem_page(0);

                // For each record
                for (uint j = 0; j < input_page->size(); ++j) {
                    Record record = input_page->get_record(j);

                    // Compute h2
                    uint h2 = record.probe_hash() % hash_table_size;

                    // Insert record into hash table at mem->mem_page(h2 + 1)
                    Page* hash_page = hash_table[h2];
                    hash_page->loadRecord(record);
                }
            }

            // Now probe left relation
            vector<uint> left_pages = partition.get_left_rel();
            for (uint page_id : left_pages) {
                // Load page into mem->mem_page(0)
                mem->loadFromDisk(disk, page_id, 0);
                Page* input_page = mem->mem_page(0);

                // For each record
                for (uint j = 0; j < input_page->size(); ++j) {
                    Record record = input_page->get_record(j);

                    // Compute h2
                    uint h2 = record.probe_hash() % hash_table_size;

                    // Get matching records from hash table
                    Page* hash_page = hash_table[h2];

                    // For each record in hash_page
                    for (uint m = 0; m < hash_page->size(); ++m) {
                        Record right_record = hash_page->get_record(m);

                        // Check if records match
                        if (right_record == record) { // overloaded == operator

                            // Output the pair
                            Page* output_page = mem->mem_page(MEM_SIZE_IN_PAGE - 1);
                            output_page->loadPair(record, right_record);

                            // If output page is full, flush to disk
                            if (output_page->full()) {
                                uint new_disk_page_id = mem->flushToDisk(disk, MEM_SIZE_IN_PAGE - 1);
                                disk_pages.push_back(new_disk_page_id);
                            }
                        }
                    }
                }
            }

            // Flush  remaining output records
            Page* output_page = mem->mem_page(MEM_SIZE_IN_PAGE - 1);
            if (!output_page->empty()) {
                uint new_disk_page_id = mem->flushToDisk(disk, MEM_SIZE_IN_PAGE - 1);
                disk_pages.push_back(new_disk_page_id);
            }
        }
    }

    return disk_pages;

}

#endif
