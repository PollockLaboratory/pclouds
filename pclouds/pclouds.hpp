#pragma once

#include "../../threadpool/ThreadPool.h"
#include "../../sequence/sequence/sequence.hpp"

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
	sequence::kmer root{ };
	std::vector <sequence::kmer> kmers{ };

	public:
	cloud( const sequence::kmer &root );
	sequence::kmer get_root_kmer( ) const;
	void add_kmer( const sequence::kmer& );

	void json( nlohmann::json );
	nlohmann::json json( );
};

std::vector <cloud> build( const parameters param, const std::string &sequence );
std::vector <cloud> generate( const std::vector <sequence::kmer> &kmers, const unsigned int &core_threshold );
void expand( std::vector <cloud> &clouds, const std::vector <sequence::kmer> &kmers, const unsigned int &shell_threshold );
void recurse_edges(const unsigned int &threshold, cloud &cloud, ThreadPool &pool,
		std::vector <sequence::kmer>::const_iterator start,
		std::vector <sequence::kmer>::const_iterator end);


}
}
