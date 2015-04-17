#include <string>
#include <unordered_map>

namespace grizzly {
namespace pclouds {

// options:
// - file location
// - cloud location
// - reverse complement
// - merge
// - print
struct kmer_map {
	std::unordered_map <std::string, std::string> all;
	std::unordered_map <std::string, std::string> rev;
	int length;
};

std::string reverse_compliment ( const std::string &sequence );


}
}
