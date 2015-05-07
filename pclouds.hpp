#pragma once

#include "../threadpool/ThreadPool.h"
#include "../sequence/sequence.hpp"

#include <json.hpp>

#include <vector>
#include <map>
#include <algorithm>

namespace grizzly {
namespace pclouds {

struct parameters {
	unsigned int core{ };
	unsigned int shell{ };
	unsigned int kmer_length{ };
};

class cloud {
	sequence::kmer root;
	std::vector <sequence::kmer> kmers{ };

	public:
	cloud( const sequence::kmer &root );
	sequence::kmer get_root_kmer( ) const;
	void add_kmer( sequence::kmer );

	void json( nlohmann::json );
	nlohmann::json json( );
};

std::vector <cloud> build( const parameters param, const std::string &sequence );
// get starting clouds
std::vector <cloud> generate_clouds( const unsigned int &core_threshold, const std::vector <sequence::kmer> &kmers );
// find any new clouds and edges
std::vector <cloud> expand_clouds( const unsigned int &shell_threshold, const std::vector <cloud> &clouds, const std::vector <sequence::kmer> &kmers );
void recurse_edges(const unsigned int &threshold, cloud &cloud, ThreadPool &pool,
		std::vector <sequence::kmer>::const_iterator start,
		std::vector <sequence::kmer>::const_iterator end);

std::vector <cloud> build( const parameters param, const std::string &sequence ) {
	std::vector <sequence::kmer> kmers = sequence::count( sequence, param.kmer_length );

	std::vector <cloud> clouds = generate_clouds( param.core, kmers );
	clouds = expand_clouds( param.shell, clouds, kmers );

	return clouds;
}


std::vector <cloud> generate_clouds( const unsigned int &core_threshold, const std::vector <sequence::kmer> &kmers ) {
	auto sorted = sequence::sort( kmers );

	std::vector <cloud> clouds;
	for (const auto &kmer: sorted) {
		if (kmer.get_count() > core_threshold)
			clouds.push_back( cloud( kmer ));
	};
	return clouds;
}

std::vector <cloud> expand_clouds( const unsigned int &shell_threshold, const std::vector <cloud> &clouds, const std::vector <sequence::kmer> &kmers ) {
	unsigned int threads = std::thread::hardware_concurrency( );
	if ( threads > 1 )
		threads = threads - 1;
	ThreadPool pool( threads );

	// after setting up our thread pool, we then want to recursively search cloud space
	// for each cloud. So we search and if we find another kmer that should be a cloud,
	// we call it a cloud, initialize it, and add it to our thread pool of clouds to search
	std::vector <cloud> new_clouds = clouds;
	for (auto &cloud: new_clouds) {
		// this is the function that is recursive
		recurse_edges( shell_threshold, cloud, pool, kmers.cbegin(), kmers.cend() );
	}
	// right now, it's possible (and will happen) that a core kmer will become part of another cloud
	// we may want to change this behavior later, but right now it doesn't matter since in the end
	// all we care about is if the kmer is there or not.

	return new_clouds;
}

void recurse_edges( const unsigned int &threshold, cloud &cloud, ThreadPool &pool,
		std::vector <sequence::kmer>::const_iterator start,
		std::vector <sequence::kmer>::const_iterator end ) {
	// algorithm: if <= 3 away and count > threshold, add to this->kmers and start new recurse at position
	// finish: once next kmer count is below threshold

	std::mutex cloud_lock;
	sequence::kmer root_kmer = cloud.get_root_kmer();
	for ( auto kmer_iter = start; kmer_iter < end; ++kmer_iter ) {
		if ( sequence::distance (kmer_iter->get_sequence(), root_kmer.get_sequence()) <= 3 && kmer_iter->get_count() > threshold ) {
			// lock the cloud so we don't get a race condition in any of our branching recurse_edges-es
			cloud_lock.lock();
			cloud.add_kmer( *kmer_iter );
			cloud_lock.unlock();

			auto status = pool.enqueue( recurse_edges, threshold, cloud, pool, ++kmer_iter, end );
			status.get(); // wait for current recursive branch to finish
			// I don't like my logic here. I would prefer to move this up into expand_clouds
			// Plus, shouldn't the pool automatically finish and start new jobs on its own?
			// It might, in which case, status.get() just means wait if it hasn't already finished
		}
	}
}

cloud::cloud( const sequence::kmer &root ) {
	this->root = root;
}

}
}
