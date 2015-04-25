#pragma once

#include "../threadpool/ThreadPool.h"
#include "../sequence/sequence.hpp"

#include <json.hpp>

#include <set>
#include <map>
#include <algorithm>

namespace grizzly {
namespace pclouds {

struct parameters {
	unsigned int core;
	unsigned int shell;
};

class cloud : public sequence::kmer {
	sequence::kmer root;
	std::set <sequence::kmer> kmers;

	public:
	cloud( const sequence::kmer &root );
	std::set <sequence::kmer> find_edge( const unsigned int shell, const std::set <sequence::kmer> &kmers ) const;

	void json( nlohmann::json );
	nlohmann::json json( );
};

std::set <cloud> get_core( const unsigned int &core_threshold, const std::set <sequence::kmer> &kmers );

std::set <cloud> build( const parameters param, const std::string &sequence ) {
	std::set <sequence::kmer> kmers = sequence::count( sequence );

	std::set <cloud> clouds = get_core( param.core, kmers );
	for (auto &cloud: clouds) {
		kmers = cloud.find_edge( param.shell, kmers );
	}

	return clouds;
}

std::set <sequence::kmer> cloud::find_edge( const unsigned int shell, const std::set <kmer> &kmers ) const {
	unsigned int threads = std::thread::hardware_concurrency( );
	if ( threads > 1 )
		threads = threads - 1;
	ThreadPool pool( threads );

	std::set <sequence::kmer> leftover = kmers;
	// recurse function
	// algorithm: pursue any that are at least 3 nucleotides away
	// branch: if < 3 away, start new recurse at position
	// finish: once count is below threshold, stop


	leftover.erase( kmers.find( this->root) );
	std::remove_if( leftover.begin(), leftover.end(),
			[&](const std::set<sequence::kmer>::iterator kmer_iter){
				return leftover.count( *kmer_iter );
			});
	return leftover;
}


}
}
