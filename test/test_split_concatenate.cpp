//! Copyright : Apache 2.0, see LICENSE 
//! 
#include <iostream>
#include <array>
#include <vector>
#include "../include/auxiliary/AuxiFunctions.hxx"

int main() {
    std::vector<double> test_vec = { 0, 2, 24, 90, 2, 5, 49, 38, 34, 75 };
    auto res = sac_de_billes::auxi_function::split(3, test_vec);
    std::vector<double> res1_ = { 0, 2, 24, 90 };
    std::vector<double> res2_ = { 2, 5, 49 };
    std::vector<double> res3_ = { 38, 34, 75 };
    std::vector<std::vector<double>> res_ref = { res1_, res2_, res3_ };
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < res_ref[i].size(); j++) {
            if (res[i][j] != res_ref[i][j]) {
                sac_de_billes::auxi_function::writeVectorToString(res[i], std::cerr, ",");
                std::cerr << " : ";
                sac_de_billes::auxi_function::writeVectorToString(res_ref[i], std::cerr, ",");
                std::cerr << "\n";
                throw std::runtime_error("Impossible! --");
            }
        }
    }

    auto res2 = sac_de_billes::auxi_function::concatenate(res);
    if (res2.size() != test_vec.size()) {
        throw std::runtime_error("Impossible!");
    }
    for (size_t i = 0; i < test_vec.size(); i++) {
        if (test_vec[i] != res2[i]) {
            throw std::runtime_error("Impossible_2 !");
        }
    }
    return EXIT_SUCCESS;
}
