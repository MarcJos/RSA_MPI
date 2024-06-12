//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once
#include <random>
#include <iostream>

namespace law {
template<class TYPE = double>
struct uniform {
	static_assert(std::is_same<TYPE, double>::value or std::is_same<TYPE, int>::value);
	using unif_distrib =
		typename std::conditional<std::is_same<TYPE, double>::value,
		std::uniform_real_distribution<>,
		std::uniform_int_distribution<>>::type;
public:
	//! @brief constructor : uniform law in [a_min, a_max]
	uniform(TYPE a_min, TYPE a_max, std::mt19937& a_gen) : gen{ &a_gen } {
		dis = unif_distrib(a_min, a_max); // takes boundaries
	}
	//! @brief default constructor : uniform law in (0, 1)
	uniform(std::mt19937& a_gen) : uniform<TYPE>(0, 1, a_gen) {}
	//! @brief constructor
	//! @param a_gen
	//! @param a_dis
	uniform(std::mt19937& a_gen, unif_distrib a_dis) : gen{ &a_gen }, dis(a_dis) {}
	//! @return evaluation of a random value
	TYPE operator()() {
		return dis(*gen);
	}

private:
	// should be initialized
	std::mt19937* gen;
	unif_distrib dis;
};

inline std::mt19937 create_random_generator(size_t i) {
	std::mt19937 random_generator;
	std::seed_seq seq{ static_cast<size_t>(rsa_mpi::get_my_rank()), i };
	random_generator.seed(seq);
	return random_generator;
}

} // namespace law
