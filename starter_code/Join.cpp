#include "Join.hpp"
#include "Page.hpp"   // Include if Page is used
#include "Record.hpp" // Include if Record is used
#include <utility>    // For std::pair
#include <vector>

using namespace std;

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
std::vector<Bucket> partition(Disk* disk, Mem* mem, pair<uint, uint> left_rel, pair<uint, uint> right_rel) {
	// Initialize partitions
	vector<Bucket> partitions(MEM_SIZE_IN_PAGE - 1, Bucket(disk));

	// Process left relation
	for (uint page_id = left_rel.first; page_id < left_rel.second; ++page_id) {
		// Load into mem page 0
		mem->loadFromDisk(disk, page_id, 0);
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
				partitions[h1].add_left_rel_page(new_disk_page_id);
				output_page->reset();
			}
		}
	}
	// Flush any remaining records in output buffers for left relation
	for (uint h1 = 0; h1 < MEM_SIZE_IN_PAGE - 1; ++h1) {
		Page* output_page = mem->mem_page(h1 + 1);
		if (!output_page->empty()) {
			uint new_disk_page_id = mem->flushToDisk(disk, h1 + 1);
			partitions[h1].add_left_rel_page(new_disk_page_id);
			output_page->reset();
		}
	}

	// Process right relation
	for (uint page_id = right_rel.first; page_id < right_rel.second; ++page_id) {
		// Load into mem page 0
		mem->loadFromDisk(disk, page_id, 0);
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
				output_page->reset();
			}
		}
	}
	// Flush any remaining records in output buffers for right relation
	for (uint h1 = 0; h1 < MEM_SIZE_IN_PAGE - 1; ++h1) {
		Page* output_page = mem->mem_page(h1 + 1);
		if (!output_page->empty()) {
			uint new_disk_page_id = mem->flushToDisk(disk, h1 + 1);
			partitions[h1].add_right_rel_page(new_disk_page_id);
			output_page->reset();
		}
	}

	return partitions;
}

/*
 * probe function
 *
 * Input:
 * disk: pointer of Disk object
 * mem: pointer of Memory object
 * partitions: a reference to a vector of buckets from the partition function
 *
 * Output:
 * A vector of page ids that contains the join result.
 */
std::vector<uint> probe(Disk* disk, Mem* mem, vector<Bucket>& partitions) {
	vector<uint> disk_pages;

	// Initialize the output page outside the partition loop
	Page* output_page = mem->mem_page(MEM_SIZE_IN_PAGE - 1);
	output_page->reset();

	for (uint i = 0; i < partitions.size(); ++i) {
		Bucket& partition = partitions[i];

		// Skip empty partitions
		if (partition.num_left_rel_record == 0 || partition.num_right_rel_record == 0) {
			continue;
		}

		// Find the smaller relation
		bool left_smaller = (partition.num_left_rel_record <= partition.num_right_rel_record);

		uint hash_table_size = MEM_SIZE_IN_PAGE - 2; // Pages 1 to MEM_SIZE_IN_PAGE - 2 are for hash table
		uint hash_table_capacity = hash_table_size * RECORDS_PER_PAGE;

		if ((left_smaller && partition.num_left_rel_record > hash_table_capacity)
		    || (!left_smaller && partition.num_right_rel_record > hash_table_capacity)) {
			continue;
		}

		// Reset only the necessary pages (input and hash table pages)
		mem->mem_page(0)->reset(); // Input page
		for (uint k = 1; k <= hash_table_size; ++k) {
			mem->mem_page(k)->reset(); // Hash table pages
		}
		// Do NOT reset the output page here

		if (left_smaller) {
			// Build hash table from left relation
			vector<Page*> hash_table(hash_table_size, nullptr);
			for (uint k = 0; k < hash_table_size; ++k) {
				hash_table[k] = mem->mem_page(k + 1);
				// hash_table[k]->reset(); // Already reset above
			}

			// Build phase
			vector<uint> left_pages = partition.get_left_rel();
			for (uint page_id : left_pages) {
				mem->loadFromDisk(disk, page_id, 0);
				Page* input_page = mem->mem_page(0);

				for (uint j = 0; j < input_page->size(); ++j) {
					Record record = input_page->get_record(j);

					// Compute h2
					uint h2 = record.probe_hash() % hash_table_size;

					Page* hash_page = hash_table[h2];
					hash_page->loadRecord(record);
				}
			}

			// Probe phase
			vector<uint> right_pages = partition.get_right_rel();
			for (uint page_id : right_pages) {
				mem->loadFromDisk(disk, page_id, 0);
				Page* input_page = mem->mem_page(0);

				for (uint j = 0; j < input_page->size(); ++j) {
					Record record = input_page->get_record(j);

					// Compute h2
					uint h2 = record.probe_hash() % hash_table_size;

					Page* hash_page = hash_table[h2];

					for (uint m = 0; m < hash_page->size(); ++m) {
						Record left_record = hash_page->get_record(m);

						// Check if records match (overloaded '==' operator)
						if (left_record == record) {
							// Output the pair
							output_page->loadPair(left_record, record);

							// If output page is full, flush to disk
							if (output_page->full()) {
								uint new_disk_page_id = mem->flushToDisk(disk, MEM_SIZE_IN_PAGE - 1);
								disk_pages.push_back(new_disk_page_id);
								output_page->reset();
							}
						}
					}
				}
			}

			// Do NOT flush output page here
		} else {
			// Build hash table from right relation
			vector<Page*> hash_table(hash_table_size, nullptr);
			for (uint k = 0; k < hash_table_size; ++k) {
				hash_table[k] = mem->mem_page(k + 1);
				// hash_table[k]->reset(); // Already reset above
			}

			// Build phase
			vector<uint> right_pages = partition.get_right_rel();
			for (uint page_id : right_pages) {
				mem->loadFromDisk(disk, page_id, 0);
				Page* input_page = mem->mem_page(0);

				for (uint j = 0; j < input_page->size(); ++j) {
					Record record = input_page->get_record(j);

					// Compute h2
					uint h2 = record.probe_hash() % hash_table_size;

					Page* hash_page = hash_table[h2];
					hash_page->loadRecord(record);
				}
			}

			// Probe phase
			vector<uint> left_pages = partition.get_left_rel();
			for (uint page_id : left_pages) {
				mem->loadFromDisk(disk, page_id, 0);
				Page* input_page = mem->mem_page(0);

				for (uint j = 0; j < input_page->size(); ++j) {
					Record record = input_page->get_record(j);

					// Compute h2
					uint h2 = record.probe_hash() % hash_table_size;

					Page* hash_page = hash_table[h2];

					for (uint m = 0; m < hash_page->size(); ++m) {
						Record right_record = hash_page->get_record(m);

						// Check if records match (overloaded '==' operator)
						if (right_record == record) {
							// Output the pair
							output_page->loadPair(record, right_record);

							// If output page is full, flush to disk
							if (output_page->full()) {
								uint new_disk_page_id = mem->flushToDisk(disk, MEM_SIZE_IN_PAGE - 1);
								disk_pages.push_back(new_disk_page_id);
								output_page->reset();
							}
						}
					}
				}
			}

			// Do NOT flush output page here
		}
	}

	// After processing all partitions, flush any remaining output records
	if (!output_page->empty()) {
		uint new_disk_page_id = mem->flushToDisk(disk, MEM_SIZE_IN_PAGE - 1);
		disk_pages.push_back(new_disk_page_id);
		output_page->reset();
	}

	return disk_pages;
}