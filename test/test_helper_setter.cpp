//! Copyright : Apache 2.0, see LICENSE 
//! 
#include <iostream>
#include <array>
#include <vector>
#include <helper.hxx>

template<int DIM>
bool check_setter() {
	constexpr int N = 1000;
	constexpr int id = 346;
	static_assert(id < N);
	constexpr double init = 5.434;
	constexpr double val = 54645.434;
	arr_vec_double<DIM> pos(N);
	for (int i = 0; i < N; i++) {
		for (int d = 0; d < DIM; d++) {
			pos[i][d] = init;
		}
	}

	rsa_helper::helper_setter<DIM> setter;

	if constexpr (DIM == 1)	setter(id, pos, val);
	if constexpr (DIM == 2)	setter(id, pos, val, val);
	if constexpr (DIM == 3)	setter(id, pos, val, val, val);
	if constexpr (DIM == 10)	setter(id, pos, val, val, val, val, val, val, val, val, val, val);

	bool ret = true;
	for (int d = 0; d < DIM; d++) {
		if (pos[d][id] != val) {
			std::cout << " d: " << d << " DIM: " << DIM << " val: " << pos[d][id] << " instead of " << val << std::endl;
			std::cout << " helper_setter doesn't set correctly the value for dim = " << DIM << std::endl;
			return false;
		}
	}

	for (int d = 0; d < DIM; d++) {
		for (int i = 0; i < N; i++) {
			if (i != id) {
				if (pos[d][i] != init) {
					std::cout << " helper_setter doesn't set the right id for dim = " << DIM << std::endl;
					return false;
				}
			}
		}
	}
	std::cout << " test success for DIM = " << DIM << std::endl;
	return true;
}

int main(int argc, char** argv) {
	if (argc != 2) {
		std::cout << "error: wrong number of elements" << std::endl;
	}
	long test = strtol(argv[1], NULL, 10);
	bool check = false;
	switch (test) {
	case 1: check = check_setter<1>(); break;
	case 2: check = check_setter<2>(); break;
	case 3: check = check_setter<3>(); break;
	case 10: check = check_setter<10>(); break;
	default: check = false;
	}
	if (check) return EXIT_SUCCESS;
	else return check;
}
