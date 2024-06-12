//! Copyright : Apache 2.0, see LICENSE 
//! 
#include "AuxiFunctions.hxx"
//! Copyright : see license.txt
//!
//! \brief Implementation of AuxiFunctions.hxx
//
#ifndef AUXIFUNCTIONS_IXX_
#define AUXIFUNCTIONS_IXX_

namespace sac_de_billes {
namespace auxi_function {

template<class NB, class C>
NB productOf(const C& table) {
    NB res = 1;
    for (const auto& nb : table) {
        res *= nb;
    }
    return res;
}

template<class T1, class T2>
T1 fast_modulo(T1 discCoord, T2 s) {
    while (discCoord < 0) {
        discCoord += s;
    }
    while (discCoord >= s) {
        discCoord -= s;
    }
    return discCoord;
}

template<unsigned short EXPONENT, class C>
inline constexpr C puissance(const C& x) {
    if constexpr (EXPONENT == 0) {
        return 1;
    } else if constexpr (EXPONENT == 1) {
        return x;
    } else {
        constexpr unsigned short expo1 = EXPONENT / 2;
        return puissance<expo1>(x) * puissance<EXPONENT - expo1>(x);
    }
}

template<class C>
C puissance(C x, size_t exponent) {
    if (exponent == 0) {
        return 1;
    } else if (exponent == 1) {
        return x;
    } else {
        size_t expo1 = exponent / 2;
        return puissance(x, expo1) * puissance(x, exponent - expo1);
    }
}

template<typename _InputIterator1, typename _InputIterator2,
    typename _OutputIterator, typename _Function1, typename _Function2>
_OutputIterator extract_list(_InputIterator1 __first1, _InputIterator1 __last1,
    _InputIterator2 __first2, _InputIterator2 __last2,
    _OutputIterator __result, _Function1 __fun1, _Function2 __fun2) {
    sort(__first1, __last1, [__fun1](const auto& i1, const auto& i2) {
        return __fun1(i1) < __fun1(i2);
        });
    sort(__first2, __last2, [__fun2](const auto& i1, const auto& i2) {
        return __fun2(i1) < __fun2(i2);
        });
    while (__first1 != __last1 and __first2 != __last2) {
        auto f1 = __fun1(*__first1);
        auto f2 = __fun2(*__first2);
        if (f1 == f2) {
            *__result = *__first1;
            __result++;
            __first1++;
        } else if (f1 < f2) {
            __first1++;
        } else {
            __first2++;
        }
    }
    return __result;
}

template<class VECTOR>
inline void writeVectorToString(const VECTOR& vect, std::ostream& f, std::string separator) {
    for (auto it = vect.begin(); it != vect.end(); it++) {
        if (it != vect.begin()) f << separator;
        f << *it;
    }
}

template<class MAP, typename PREDICATE>
size_t erase_if(MAP& c, PREDICATE pred) {
    auto old_size = c.size();
    for (auto i = c.begin(), last = c.end(); i != last; ) {
        if (pred(*i)) {
            i = c.erase(i);
        } else {
            ++i;
        }
    }
    return old_size - c.size();
}

template<size_t N, class A, class B>
array<B, N> convertArray(array<A, N> array_) {
    array<B, N> result{};
    for (size_t i = 0; i < N; i++) {
        result[i] = static_cast<B>(array_[i]);
    }
    return result;
}

template<int64_t  Begin, uint64_t... I, typename Func>
void for_constexpr(Func f, std::index_sequence<I...>) {
    (f(I + Begin), ...);
}

template<int64_t Begin, int64_t End, typename Func>
void for_constexpr(Func f) {
    for_constexpr<Begin>(f, std::make_index_sequence<End - Begin>{});
}


template<class C>
vector<vector<C>> split(size_t nb_of_parts, const vector<C>& vector_to_be_splitted) {
    std::cerr << "begin split" << std::endl;
    size_t total_size = vector_to_be_splitted.size();
    size_t average_size = total_size / nb_of_parts;
    vector<size_t> local_size(nb_of_parts, average_size);
    int64_t reste = total_size - nb_of_parts * average_size;
    for (int64_t i = 0; i < reste; i++) {
        local_size[i] += 1;
    }

    vector<vector<C>> result(nb_of_parts, vector<C>{});

    size_t current_index = 0;
    for (size_t i = 0; i < nb_of_parts; i++) {
        result[i].resize(local_size[i]);
        std::copy(vector_to_be_splitted.begin() + current_index, vector_to_be_splitted.begin() + current_index + local_size[i], result[i].begin());
        current_index += local_size[i];
    }
    if (current_index != total_size) {
        throw runtime_error("Incorrect splitting");
    }
    std::cerr << "end split" << std::endl;
    return result;
}

template<class C>
vector<C> concatenate(const vector<vector<C>>& vectors_to_be_gathered) {
    vector<C> result{};
    for (size_t i = 0; i < vectors_to_be_gathered.size(); i++) {
        size_t size_0 = result.size();
        const auto& to_be_copied = vectors_to_be_gathered[i];
        result.resize(size_0 + to_be_copied.size());
        std::copy(to_be_copied.begin(), to_be_copied.end(), result.begin() + size_0);
    }
    return result;
}

} // namespace auxi_function
} // namespace sac_de_billes


#endif /* AUXIFUNCTIONS_IXX_ */
