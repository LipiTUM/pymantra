#ifndef ENRICHMENT_OBJECTIVE_FUNCTIONS_HPP
#define ENRICHMENT_OBJECTIVE_FUNCTIONS_HPP


#include <string>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include "reaction_graph.hpp"
#include "../Exceptions.hpp"
#include "../statsutils.hpp"


using std::string;
using std::vector;


// subgraph
typedef boost::unordered_set<size_t> vertex_subgraph;

// utility functions
double dist(double x, double y);
double lp_norm(vector<double>& vec, float p);

// metabolomics
void set_reaction_kcl_objectives(vector<string>&, Graph&);


void set_reaction_ld_objectives(vector<string>&, Graph&);
double reaction_dysregulation(vector<string>&, Graph&, vertex_subgraph&, float p);

void precomputed_objectives(vector<string>&, Graph&);

// multi-omics
vector<double> residual_zscores(vector<double>&, vector<double>&);
void set_multiomics_reaction_dysregulation_objectives(vector<string> &groups, Graph &graph,
                                                      const string& edge_type);
double multiomics_reaction_dysregulation(vector<string> &groups, Graph &graph, vertex_subgraph &subnetwork,
                                         const string& node_type, float p);


void set_microbiome_reaction_dysregulation_objectives(vector<string>&, Graph&);
double microbiome_reaction_dysregulation(vector<string>&, Graph&, vertex_subgraph&, float);


void set_transcriptome_reaction_dysregulation_objectives(vector<string>&, Graph&);
double transcriptome_reaction_dysregulation(vector<string>&, Graph&, vertex_subgraph&, float);


#endif //ENRICHMENT_OBJECTIVE_FUNCTIONS_HPP
