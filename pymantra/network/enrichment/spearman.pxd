# distutils: language = c++

from libcpp cimport bool as bool_t


cdef extern from "statsutils.hpp":
    tuple spearmans_rank_1by1(list, list) except +
    tuple spearmans_rank_2by2(list, list, int) except +
    bool_t HAS_PARALLEL()
