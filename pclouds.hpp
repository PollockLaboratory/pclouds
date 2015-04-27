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
	unsigned int core;
	unsigned int shell;
};

class cloud {
	const sequence::kmer root;
	std::vector <sequence::kmer> kmers;

	std::vector <sequence::kmer>
		recurse_edges( const std::vector <sequence::kmer>::iterator &begin,
				const std::vector <sequence::kmer>::iterator &end,
				const unsigned int threshold );

	public:
	cloud( const sequence::kmer &root );
	sequence::kmer get_root_kmer( ) const;
	void add_kmer( sequence::kmer );

	void json( nlohmann::json );
	nlohmann::json json( );
};

std::vector <cloud> build( const parameters param, const std::string &sequence );
std::vector <cloud> generate_clouds( const unsigned int &core_threshold, const std::vector <sequence::kmer> &kmers );
std::vector <cloud> expand_clouds( const unsigned int &shell_threshold, const std::vector <cloud> &clouds, const std::vector <sequence::kmer> &kmers );
void expand_cloud(const unsigned int &threshold, cloud &cloud, ThreadPool &pool,
		std::vector <sequence::kmer>::const_iterator start,
		std::vector <sequence::kmer>::const_iterator end);

std::vector <cloud> build( const parameters param, const std::string &sequence ) {
	std::vector <sequence::kmer> kmers = sequence::count( sequence );

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

	std::vector <cloud> new_clouds = clouds;
	for (auto &cloud: new_clouds) {
		expand_cloud( shell_threshold, cloud, pool, kmers.cbegin(), kmers.cend() );
	}
	return new_clouds;
}

void expand_cloud( const unsigned int &threshold, cloud &cloud, ThreadPool &pool,
		std::vector <sequence::kmer>::const_iterator start,
		std::vector <sequence::kmer>::const_iterator end ) {
	// algorithm: if <= 3 away and count > threshold, add to this->kmers and start new recurse at position
	// finish: once next kmer count is below threshold

	std::mutex cloud_lock;
	sequence::kmer root_kmer = cloud.get_root_kmer();
	for ( auto kmer = start; kmer < end; ++kmer ) {
		if ( sequence::distance (kmer->get_sequence(), root_kmer.get_sequence()) <= 3 && kmer->get_count() > threshold ) {
			cloud_lock.lock(); // move to add_kmer ?
			cloud.add_kmer( *kmer );
			cloud_lock.unlock();

			pool.enqueue( expand_cloud, threshold, cloud, pool, ++kmer, end );
		}
	}
}

}
}
