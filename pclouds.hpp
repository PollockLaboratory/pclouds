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
	const sequence::kmer root;
	mutable std::set <sequence::kmer> kmers;

	public:
	cloud( const sequence::kmer &root );
	std::set <sequence::kmer> populate_shell( const unsigned int shell, const std::set <sequence::kmer> &kmers ) const;

	void json( nlohmann::json );
	nlohmann::json json( );
};

std::set <cloud> get_core( const unsigned int &core_threshold, const std::set <sequence::kmer> &kmers );

std::set <cloud> build( const parameters param, const std::string &sequence ) {
	std::set <sequence::kmer> kmers = sequence::count( sequence );

	std::set <cloud> clouds = get_core( param.core, kmers );
	for (auto &cloud: clouds) {
		kmers = cloud.populate_shell( param.shell, kmers );
	}

	return clouds;
}

std::set <sequence::kmer>
recurse_edges( const std::set <sequence::kmer>::iterator &begin, const std::set <sequence::kmer>::iterator &end, const unsigned int threshold );

std::set <sequence::kmer> cloud::populate_shell( const unsigned int shell_threshold, const std::set <kmer> &kmers ) const {
	this->kmers = recurse_edges( kmers.cbegin(), kmers.cend(), shell_threshold );

	std::set <sequence::kmer> leftover = kmers;
	leftover.erase( kmers.find( this->root) );
	std::remove_if( leftover.begin(), leftover.end(),
			[&](const std::set<sequence::kmer>::iterator kmer_iter){
				return this->kmers.count( *kmer_iter );
			});
	return leftover;
}

std::set <sequence::kmer>
recurse_edges( const std::set <sequence::kmer> &begin, const std::set <sequence::kmer> &end, const unsigned int threshold ) {
	// algorithm: pursue any that are at least 3 nucleotides away
	// branch: if < 3 away, add to this->kmers and start new recurse at position
	// finish: once count is below threshold, stop
	unsigned int threads = std::thread::hardware_concurrency( );
	if ( threads > 1 )
		threads = threads - 1;
	ThreadPool pool( threads );
	// for (  iter = kmer_inc; iter < kmer_end; iter++ )
		// if ( seq::dist (iter, iter inc)j < max_dist && iter.count > shell_thresh )
			// lock mutex
			// add kmer to iter_th_safe_container
			// unlock mutex
			// add (recurse (iter, kmer_begin, kmer_end) )
		// else
			// thread.join
}

}
}
