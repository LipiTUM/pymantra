# What requirements do we need here?
cmake_minimum_required(VERSION 3.20.1)
project(enrichment)

set(CMAKE_CXX_STANDARD 17)

# TODO: are these fully distribution independent?

# OpenMP
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    add_compile_definitions(USE_OPENMP)
endif()

# python header
# set(PYTHON_HEADERS $(shell python -c 'from distutils.sysconfig import get_python_inc; print(get_python_inc())'))
# set(PYTHON_STATIC $(shell python -c 'from distutils import sysconfig; print(sysconfig.get_config_var("LIBDIR"))'))
find_package ( PythonLibs COMPONENTS Development Interpreter REQUIRED )
include_directories( ${PYTHON_INCLUDE_DIRS} )

# boost
find_package( Boost 1.70 COMPONENTS graph REQUIRED )
include_directories( ${Boost_INCLUDE_DIR} )


### library
add_library(
    enrichment SHARED
        Exceptions.hpp Exceptions.cpp
        pyutils.hpp pyutils.cpp
        utils.hpp utils.cpp
        statsutils.hpp statsutils.cpp
        RFA/ReferenceFree.hpp RFA/ReferenceFree.cpp
        RBA/ReferenceBased.hpp RBA/ReferenceBased.cpp
        RFA/ReferenceNetwork.hpp RFA/ReferenceNetwork.cpp
        LSO/local_dependence.hpp LSO/local_dependence.cpp
        LSO/LocalSearch.hpp LSO/LocalSearch.tpp LSO/LocalSearch.cpp
        LSO/reaction_graph.hpp LSO/reaction_graph.cpp
        LSO/objective_functions.hpp LSO/objective_functions.cpp LSO/lso_utils.cpp)
target_link_libraries( enrichment PRIVATE OpenMP::OpenMP_CXX )
target_link_libraries( enrichment PRIVATE ${PYTHON_LIBRARIES} )
target_link_libraries( enrichment PRIVATE ${Boost_LIBRARIES} )


### Doctest testing
# ReferenceNetwork
add_executable(
    ReferenceNetworkTest
        tests/test_ReferenceNetwork.cpp
)
add_test( NAME ReferenceNetworkTest COMMAND ReferenceNetworkTest )
target_link_libraries( ReferenceNetworkTest PRIVATE OpenMP::OpenMP_CXX )
target_link_libraries( ReferenceNetworkTest PRIVATE ${PYTHON_LIBRARIES} )


# Stats Functions
add_executable(
        StatsTests
        Exceptions.hpp Exceptions.cpp
        pyutils.hpp pyutils.cpp
        LSO/objective_functions.hpp LSO/objective_functions.cpp
        statsutils.hpp statsutils.cpp
        tests/test_stats.cpp
)
add_test( NAME StatsTest COMMAND StatsTest )
target_link_libraries( StatsTests PRIVATE ${Boost_LIBRARIES} )
target_link_libraries( StatsTests PRIVATE ${PYTHON_LIBRARIES} )

# ReferenceFree
add_executable(
        ReferenceFreeTest
        tests/test_ReferenceFree.cpp
)
add_test( NAME ReferenceFreeTest COMMAND ReferenceFreeTest )
target_link_libraries( ReferenceFreeTest PRIVATE OpenMP::OpenMP_CXX )
target_link_libraries( ReferenceFreeTest PRIVATE ${PYTHON_LIBRARIES} )

### Local Search related tests
# objective functions
add_executable(
        objective_functions_test
        statsutils.hpp statsutils.cpp
        LSO/reaction_graph.hpp LSO/reaction_graph.cpp
        pyutils.hpp pyutils.cpp
        Exceptions.hpp Exceptions.cpp
        LSO/objective_functions.hpp LSO/objective_functions.cpp
        tests/test_objective_functions.cpp
)
add_test( NAME objective_functions_test COMMAND objective_functions_tests )
target_link_libraries( objective_functions_test PRIVATE ${PYTHON_LIBRARIES} )
target_link_libraries( objective_functions_test PRIVATE ${Boost_LIBRARIES} )

# reaction graphs
add_executable(
        reaction_graph_test
        LSO/reaction_graph.hpp LSO/reaction_graph.cpp
        pyutils.hpp pyutils.cpp
        Exceptions.hpp Exceptions.cpp
        tests/test_reaction_graph.cpp
)
add_test( NAME reaction_graph_test COMMAND reaction_graph_test )
target_link_libraries( reaction_graph_test PRIVATE ${PYTHON_LIBRARIES} )
target_link_libraries( reaction_graph_test PRIVATE ${Boost_LIBRARIES} )

# sorted merging
add_executable(
        sorted_merging_test
        statsutils.hpp statsutils.cpp
        pyutils.hpp pyutils.cpp
        Exceptions.hpp Exceptions.cpp
        LSO/LocalSearch.hpp LSO/LocalSearch.tpp LSO/LocalSearch.cpp
        LSO/reaction_graph.hpp LSO/reaction_graph.cpp
        LSO/objective_functions.hpp LSO/objective_functions.cpp
        tests/test_sorted_merging.cpp
)
add_test( NAME sorted_merging_test COMMAND sorted_merging_test )
target_link_libraries( sorted_merging_test PRIVATE OpenMP::OpenMP_CXX )
target_link_libraries( sorted_merging_test PRIVATE ${PYTHON_LIBRARIES} )
target_link_libraries( sorted_merging_test PRIVATE ${Boost_LIBRARIES} )

# string formatting
add_executable(
        StrFormatTest
        Exceptions.hpp Exceptions.cpp
        tests/test_str_formatting.cpp
)
add_test( NAME StrFormatTest COMMAND StrFormatTest )

# LocalSearch
add_executable(
        LocalSearchTest
        LSO/LocalSearch.hpp LSO/LocalSearch.tpp LSO/LocalSearch.cpp
        LSO/reaction_graph.hpp LSO/reaction_graph.cpp
        LSO/objective_functions.hpp LSO/objective_functions.cpp
        pyutils.hpp pyutils.cpp
        statsutils.hpp statsutils.cpp
        Exceptions.hpp Exceptions.cpp
        tests/test_LocalSearch.cpp
        LSO/lso_utils.cpp)
add_test( NAME LocalSearchTest COMMAND LocalSearchTest )
target_link_libraries( LocalSearchTest PRIVATE OpenMP::OpenMP_CXX )
target_link_libraries( LocalSearchTest PRIVATE ${PYTHON_LIBRARIES} )
target_link_libraries( LocalSearchTest PRIVATE ${Boost_LIBRARIES} )


### package export
export(PACKAGE enrichment)
