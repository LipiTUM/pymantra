# distutils: language = c++
from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.string cimport string


cdef extern from "LocalSearch.hpp":

    cdef cppclass LocalSearch:
        LocalSearch(
            dict, dict, int, double, double,
            size_t, size_t, size_t, size_t, float
        ) except +

        double score_solution(vector[string]&)
        #
        void run_local_search(vector[string]&, int) except +
        #
        void set_seed_py(string&, size_t) except +
        #
        void set_seed_py(list) except +
        #
        void precompute_objectives_py(vector[string]&) except +

        bool get_converged() except +
        double get_best_score() except +
        set get_best_solution() except +
        vector[double] get_score_progression() except +

        # getters
        double get_temp()
        size_t get_lmax()
        size_t get_lmin()
        size_t get_min_reactions()
        size_t get_maxiter()
        double get_deltamin()
        float get_lp_norm()
        dict get_objective_values() except +

        # setters
        void set_temp(double)
        void set_lmax(size_t)
        void set_lmin(size_t)
        void set_min_reactions(size_t)
        void set_maxiter(size_t)
        void set_deltamin(double)
        void set_lp_norm(float)


cdef extern from "reaction_graph.hpp":
    set extract_reaction_graph(dict, dict, bool) except +
