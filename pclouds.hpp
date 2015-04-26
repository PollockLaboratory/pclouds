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

	void json( nlohmann::json );
	nlohmann::json json( );
};

std::vector <cloud> build( const parameters param, const std::string &sequence );
std::vector <cloud> generate_clouds( const unsigned int &core_threshold, const std::vector <sequence::kmer> &kmers );
std::vector <cloud> expand_clouds( const unsigned int &shell_threshold, const std::vector <cloud> &clouds, const std::vector <sequence::kmer> &kmers );

std::vector <cloud> build( const parameters param, const std::string &sequence ) {
	std::vector <sequence::kmer> kmers = sequence::count( sequence );

	std::vector <cloud> clouds = generate_clouds( param.core, kmers );
	clouds = expand_clouds( param.shell, clouds, kmers );

	return clouds;
}


std::vector <cloud> generate_clouds( const unsigned int &core_threshold, const std::vector <sequence::kmer> &kmers ) {
}

std::vector <cloud> expand_clouds( const unsigned int &shell_threshold, const std::vector <cloud> &clouds, const std::vector <sequence::kmer> &kmers ) {
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

	std::vector <sequence::kmer> leftover = kmers;
	leftover.erase( kmers.find( this->root) );
	std::remove_if( leftover.begin(), leftover.end(),
			[&](const std::vector<sequence::kmer>::iterator kmer_iter){
				return this->kmers.count( *kmer_iter );
			});
	return leftover;
}

}
}
