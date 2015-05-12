#include "pclouds.hpp"

#include "../../threadpool/ThreadPool.h"
#include "../../sequence/sequence/sequence.hpp"

#include <json.hpp>

#include <vector>
#include <map>
#include <algorithm>

namespace grizzly {
namespace pclouds {

cloud::cloud( const sequence::kmer &root ) {
	this->root = root;
}

sequence::kmer cloud::get_root_kmer( ) const {
	return this->root;
}

void cloud::add_kmer( const sequence::kmer &kmer) {
	this->kmers.push_back( kmer );
}

void cloud::json( nlohmann::json json ) {
	const auto jroot = *json.find("sequence");
	sequence::kmer root;
	root.json( jroot );
	this->root = root;
	const auto jkmers = json.find("kmers")->get<std::vector <nlohmann::json>>();
	for (const auto &jkmer: jkmers) {
		sequence::kmer kmer;
		kmer.json( jkmer );
		this->kmers.push_back( kmer );
	}
}

nlohmann::json cloud::json( ) {
	nlohmann::json json;
	json["root"] = this->root.json();
	std::vector <nlohmann::json> jkmers;
	for (const auto &kmer: this->kmers) {
		jkmers.push_back( kmer.json() );
	}
	json["kmers"] = jkmers;
	return json;
}

void recurse_edges(const unsigned int &threshold, cloud &cloud, ThreadPool &pool,
		std::vector <sequence::kmer>::const_iterator start,
		std::vector <sequence::kmer>::const_iterator end);

std::vector <cloud> build( const parameters param, const std::string &sequence ) {
	std::vector <sequence::kmer> kmers = sequence::count( sequence, param.kmer_length );

	std::vector <cloud> clouds = generate( kmers, param.core );
	expand( clouds, kmers, param.shell );

	return clouds;
}

std::vector <cloud> generate( const std::vector <sequence::kmer> &kmers, const unsigned int &core_threshold ) {
	auto sorted = kmers;
	std::sort( sorted.begin(), sorted.end() );

	std::vector <cloud> clouds;
	for (const auto &kmer: sorted) {
		if (kmer.get_count() > core_threshold)
			clouds.push_back( cloud( kmer ));
	};
	return clouds;
}

void expand(  std::vector <cloud> &clouds, const std::vector <sequence::kmer> &kmers, const unsigned int &shell_threshold ) {
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
			// I don't like my logic here. I would prefer to move this up into expand
			// Plus, shouldn't the pool automatically finish and start new jobs on its own?
			// It might, in which case, status.get() just means wait if it hasn't already finished
		}
	}
}

}
}
