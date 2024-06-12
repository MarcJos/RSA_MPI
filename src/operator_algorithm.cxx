//! Copyright : Apache 2.0, see LICENSE 
//! 
#include<operator_algorithm.hxx>

namespace algorithm {

void fill_priority(uint64_t* a_in, size_t a_size) {
    const double inf = 0;
    const double sup = uint64_t(-1) - 1;
    std::random_device rd;
    std::uniform_real_distribution<> tmp_dis(inf, sup); // takes boundaries
    std::mt19937 tmp_gen(rd());
    for (int idx = 0; idx < a_size; idx++) {
        a_in[idx] = uint64_t(tmp_dis(tmp_gen));
    }
}

} // namespace  algorithm
