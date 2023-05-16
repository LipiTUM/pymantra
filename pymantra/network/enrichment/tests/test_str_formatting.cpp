#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN


#include "../doctest.h"
#include "../Exceptions.hpp"
#include <iostream>
#include <string>
#include <boost/format.hpp>


using std::string;
using boost::format;


TEST_SUITE_BEGIN("string formatting error");

TEST_CASE("") {
    string DATA_KEY = "data";
    string node_name = "test_node";
    throw MissingAttribute(
            boost::str(format("Group data key '%1%' not found for node '%2%'") % DATA_KEY % node_name)
    );
}
