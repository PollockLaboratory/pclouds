#pragma once

#include "../sequence/sequence.hpp"

#include <json.hpp>

#include <vector>
#include <map>
#include <set>

namespace grizzly {
namespace pclouds {

struct parameters {
	unsigned int core;
	unsigned int outer;
};

class cloud {
	unsigned int id{ };
	std::vector <sequence::kmer> kmers;

	public:
	cloud( const std::string id );
	void find_edge( const unsigned int outer );
	std::vector <sequence::kmer> get_kmers( );
	unsigned int get_id( );

	void json( nlohmann::json );
	nlohmann::json json( );
};

std::vector <cloud> get_core( const unsigned int &core_threshold, const std::vector <sequence::kmer> &kmers );
void outline( std::vector <cloud> &clouds );

std::vector <cloud> build( const parameters param, const std::string &sequence ) {
	vector <sequence::kmer> kmers = sequence::count( sequence );

	std::vector <cloud> clouds = get_core( param.core, kmers );

	std::set <unsigned int> unprocessed;
	for (auto &cloud: clouds) {
		unprocessed.insert( cloud.get_id() );
	}

	for (auto &cloud: clouds) {
		if ( unprocessed.count( cloud.get_id() ) == 0 ) {
			cloud.find_edge( param.outer );
			//erase any in kmer list that were found as an edge
			unprocessed.erase( cloud.get_id( ) );
		}
	}

	return clouds;
}


}
}
