#include <cmath>
#include <iostream>

#include <boost/format.hpp>
#include <boost/math/statistics/univariate_statistics.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include "objective_functions.hpp"

using std::pow;
using std::abs;
using boost::str;
using boost::format;


typedef boost::graph_traits<Graph>::edge_iterator EdgeIter;
typedef std::pair<EdgeIter, EdgeIter> Edge;


/* ===================
 * KCL-based Objective
 * ===================
 */
double dist(double x, double y) {
   return abs(x - y) / abs(x);
}


void set_reaction_kcl_objectives(vector<string>& groups, Graph& graph) {
    Graph::vertex_iterator vertex, vend;
    for (boost::tie(vertex, vend) = boost::vertices(graph); vertex != vend; ++vertex) {
        if (graph[*vertex].type == "reaction") {
            if (graph[*vertex].value.find(groups[0]) == graph[*vertex].value.end()) {
                throw std::runtime_error(
                    str(format("Group '%1%' data not found for reaction node %2%") % groups[0] % graph[*vertex].name)
                );
            }
            if (graph[*vertex].value.find(groups[1]) == graph[*vertex].value.end()) {
                throw std::runtime_error(
                    str(format("Group '%1%' data not found for reaction node %2%") % groups[1] % graph[*vertex].name)
                );
            }
            graph[*vertex].objective = dist(
                graph[*vertex].value[groups[0]], graph[*vertex].value[groups[1]]);
        }
    }
}


double lp_norm(vector<double>& vec, float p) {
    if (std::numeric_limits<float>::infinity() == p) {
        double max = 0.;
        for (double x : vec) {
            if (abs(x) > max) max = x;
        }
        return max;
    }
    double sum = 0.;
    for (double x : vec) sum += pow(abs(x), p);
    return pow(sum, 1. / p);
}


/*
 * NOTE: these objective functions are meant for group-wise parameter estimation
 */
/**
 * Compute the pure metabolite reaction dysregulation defined as the mean the
 * of the relative differences of the reaction values between groups
 *
 * \f$ \Gamma(a, b, R) = \frac{1}{N}\sum_{r \in R} \frac{r_a - r_b}{r_a}\f$
 *
 * The input graph is assumed to be a reaction-reaction network
 *
 * @param groups Sample groups to compare, represented as an array of strings of
 *               size 2. The first element will be treated as the control group.
 * @param graph Reaction to use as a basis for objective function. For the required
 *              vertex and edge properties see ::objective_. Please note that no
 *              checks for the correctness of vertex types is done!
 * @param subnetwork Vertex subgraph specifying which subnetwork to investigate.
 *
 * @return dysregulation value between the two groups based on the input subnetwork
 */
double reaction_dysregulation(
    vector<string>& groups, Graph& graph,
    vertex_subgraph& subnetwork, float p
) {
    if (subnetwork.empty()) return 0.;
    boost::graph_traits<Graph>::vertex_descriptor vertex;
    vector<double> objectives(subnetwork.size());
    int i = 0;
    for (auto v_idx : subnetwork) {
        vertex = boost::vertex(vertex_index_t(v_idx), graph);
        objectives[i] = graph[vertex].objective;
        i++;
    }
    double dysreg = lp_norm(objectives, p);
    return dysreg / double(subnetwork.size());
}


/* ==================================
 * Linear Dependence-based Objective
 * ==================================
 */
void set_reaction_ld_objectives(vector<string>& groups, Graph& graph) {
    Graph::vertex_iterator vertex, vend;
    double pval;
    vector<string> invalid_nodes;
    for (boost::tie(vertex, vend) = boost::vertices(graph); vertex != vend; ++vertex) {
        string node_name = graph[*vertex].name;
        if (graph[*vertex].type == "reaction") {
            if (graph[*vertex].vec_value.find(groups[0]) == graph[*vertex].vec_value.end()) {
                throw std::runtime_error(
                    str(format("Group '%1%' data not found for reaction node %2%") % groups[0] % graph[*vertex].name)
                );
            }
            if (graph[*vertex].vec_value.find(groups[1]) == graph[*vertex].vec_value.end()) {
                throw std::runtime_error(
                    str(format("Group '%1%' not found for reaction node %2%") % groups[1] % graph[*vertex].name)
                );
            }
            // objective is 1 - the probability that the residual distributions are different
            try {
                pval = welch_test(
                    graph[*vertex].vec_value[groups[0]], graph[*vertex].vec_value[groups[1]]);
            } catch (const std::invalid_argument& err) {
                throw std::invalid_argument(
                        str(format("Welch test comparing residual distributions failed for reaction node %1%. %2%")
                        % graph[*vertex].name % err.what()
                    )
                );
            } catch (const std::domain_error& err) {
                // usually this means there was an error in the Welch's test computation from a
                // NaN statistics value, so we check the potential (common) sources for this
                // => too little samples, 0 variance in *both* groups
                size_t v0, v1;
                double var0 = boost::math::statistics::sample_variance(
                        graph[*vertex].vec_value[groups[0]]);
                double var1 = boost::math::statistics::sample_variance(
                        graph[*vertex].vec_value[groups[1]]);

                // if variances are both zero or any of the variance is nan or
                // inf this explains why the test failed => we can "safely"
                // return NaN
                if (var0 + var1 <= std::numeric_limits<double>::epsilon() ||
                        boost::math::isnan(var0) || boost::math::isnan(var1) ||
                        boost::math::isinf(var0) || boost::math::isinf(var1))
                {
                    pval = std::numeric_limits<double>::quiet_NaN();
                    invalid_nodes.push_back(graph[*vertex].name);
                } else {
                    v0 = graph[*vertex].vec_value[groups[0]].size();
                    v1 = graph[*vertex].vec_value[groups[1]].size();
                    if (v0 < 3 || v1 < 3) {
                            throw std::domain_error(
                            str(format("At least 3 samples per group required but found only %1% and %2% for node")
                                % v0 % v1 % graph[*vertex].name
                            )
                        );
                    }
                    throw std::domain_error(
                        str(format("Invalid value in Welch's test for node %1%. This usually means "
                                   "there are invalid values (+- NaN or += inf) in the group vectors or their variances. "
                                   "%2%")
                            % graph[*vertex].name % err.what()
                        )
                    );
                }
            }
            if (isnan(pval)) graph[*vertex].objective = 0.;
            else graph[*vertex].objective = 1 - pval;
        }
    }
    if (!invalid_nodes.empty()) {
        // logging if nodes for which Welch's test crashed with zero variance
        size_t i;
        std::ostringstream stream;
        stream << "zero variances found for the following nodes: ";
        for (i = 0; i < invalid_nodes.size() - 1; ++i)
            stream << invalid_nodes[i] << ", ";
        stream << invalid_nodes[i];
        // we could use cerr as well, since we don't need the buffering of clog,
        // I personally find it more consistent to use clog for logging
        std::clog << std::endl << stream.str() << std::endl;
    }
}


void precomputed_objectives(vector<string>& groups, Graph& graph) {
    // NOTE: this is just a placeholder for implementation simplicity
}


/* ==================================
 * Multi-omics Association Objectives
 * ==================================
 */
vector<double> residual_zscores(vector<double>& control, vector<double>& case_) {
    double mu = mean(control);
    double sigma = sqrt(var(control, mu));

    vector<double> zcase(case_.size());
    unsigned i;
    for (i = 0; i < case_.size(); i++) {
        zcase[i] = (case_[i] - mu) / sigma;
    }
    return zcase;
}


void set_multiomics_reaction_dysregulation_objectives(
        vector<string>& groups, Graph& graph, const string& edge_type) {
    // set node objectives => TODO: keep automatic option detection?
    Graph::vertex_iterator vertex, vend;
    bool data_empty, vec_data_empty;
    for (boost::tie(vertex, vend) = vertices(graph); vertex != vend; ++vertex) {
        if (graph[*vertex].type == "reaction") {
            data_empty = graph[*vertex].value.empty();
            vec_data_empty = graph[*vertex].vec_value.empty();
            if (vec_data_empty) {
                if (data_empty) {
                    throw std::runtime_error(
                        str(format("No node values found for at least one node (%1%)") % graph[*vertex].name)
                    );
                }
                set_reaction_kcl_objectives(groups, graph);
                break;
            }
            set_reaction_ld_objectives(groups, graph);
            break;
        }
    }
    // set edge objectives
    double objective;
    for (auto edge : boost::make_iterator_range(boost::edges(graph))) {
       if (graph[edge].type == edge_type) {
           // add difference in correlation between groups as objective
           // since we work with correlations ([-1, 1]) diff / 2 scales to range [0, 1]
           // TODO: should we try to take the direction of change into account?
           objective = abs(graph[edge].value[groups[0]] - graph[edge].value[groups[1]]) / 2;
           if (isnan_(objective)) graph[edge].objective = 0.;
           else graph[edge].objective = objective;
       }
    }
}


double multiomics_reaction_dysregulation(
        vector<string>& groups, Graph& graph, vertex_subgraph& subnetwork, const string& node_type, float p
) {
    if (subnetwork.empty()) return 0.;

    double reaction_dysreg, association_change;
    boost::graph_traits<Graph>::vertex_descriptor vertex;

    string edge_type;
    if (node_type == "organism") edge_type = "REACTION_ORGANISM";
    else edge_type = "REACTION_GENE";

    vector<double> reaction_objectives;
    vector<double> association_objectives;
    int i = 0;
    int n_reactions = 0;
    int n_associations = 0;
    for (auto v_idx : subnetwork) {
        vertex = boost::vertex(vertex_index_t(v_idx), graph);
        if (graph[vertex].type == "reaction") {
            reaction_objectives.push_back(graph[vertex].objective);
            n_reactions++;
            // TODO: check if type names are correct
        } else if (graph[vertex].type == node_type) {
           // go over all adjacent edges and add objectives if reaction - organism edge
           for (auto edge : boost::make_iterator_range(boost::out_edges(vertex, graph))) {
               if (graph[edge].type == edge_type) {
                   association_objectives.push_back(graph[edge].objective);
                   n_associations++;
               }
           }
        }
        i++;
    }
    // both should be in the range of [0, 1] since each individual element is in this range
    reaction_dysreg = lp_norm(reaction_objectives, p);
    association_change = lp_norm(association_objectives, p);
    if (n_reactions == 0) {
        if (n_associations == 0) return -1;
        return association_change / n_associations;
    }
    if (n_associations == 0) {
        return reaction_dysreg / n_reactions;
    }
    return reaction_dysreg / n_reactions + association_change / n_associations;
}


/* ===========================
 * Organism-Reaction Objective
 * ===========================
 */
void set_microbiome_reaction_dysregulation_objectives(vector<string>& groups, Graph& graph) {
    set_multiomics_reaction_dysregulation_objectives(groups, graph, "REACTION_ORGANISM");
}


double microbiome_reaction_dysregulation(vector<string>& groups, Graph& graph, vertex_subgraph& subnetwork, float p) {
    return multiomics_reaction_dysregulation(groups, graph, subnetwork, "organism", p);
}


/* ==================================
 * Gene/Transcript-Reaction Objective
 * ==================================
 */
void set_transcriptome_reaction_dysregulation_objectives(vector<string>& groups, Graph& graph) {
    set_multiomics_reaction_dysregulation_objectives(groups, graph, "REACTION_GENE");
}


double transcriptome_reaction_dysregulation(
    vector<string>& groups, Graph& graph,
    vertex_subgraph& subnetwork, float p
) {
    return multiomics_reaction_dysregulation(groups, graph, subnetwork, "gene", p);
}
