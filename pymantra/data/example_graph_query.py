"""
Generating the mantra example graph. The main is not meant to be used and is
only included for reproducibility of the example graph.
"""
if __name__ == "__main__":
    import os
    import pathlib
    from neo4j.graph import Node
    import networkx as nx

    from pymantra.database import NetworkGenerator, read_env

    def add_graph_data(node_set, edge_set, results):
        for conn in results:
            for elem in conn:
                if isinstance(elem, Node):
                    node_set.add((elem.get("nodeLabel"), list(elem.labels)[0]))
                else:
                    edge_set.add((
                        elem.nodes[0].get("nodeLabel"),
                        elem.nodes[1].get("nodeLabel"),
                        elem.type
                    ))

    read_env(pathlib.Path(__file__).parent.parent.parent / ".env")
    AUTH = (os.getenv("NEO4J_USER"), os.getenv("NEO4J_PASSWORD"))
    gen = NetworkGenerator("bolt://127.0.0.1:7687", AUTH)

    res = gen.session.run(
        'match (ms:metabolite)-[a]-(s)-[b]-(m:metabolite {KeggID: "C01236"})'
        '-[c]-(p)-[d]-(mp:metabolite) return *'
    ).values()

    nodes = set()
    edges = set()
    add_graph_data(nodes, edges, res)

    reaction_nodes = {node[0] for node in nodes if node[1] == "reaction"}

    join_str = "','"
    joined_r = f"'{join_str.join(reaction_nodes)}'"
    o_res = gen.session.run(
        "match (r:reaction)-[e]-(o:organism) where r.nodeLabel in "
        f"[{joined_r}] return o, e, r limit 10"
    ).values()
    add_graph_data(nodes, edges, o_res)

    g_res = gen.session.run(
        f"match (r:reaction)-[e]-(g:gene) where r.nodeLabel in [{joined_r}] "
        "return g, e, r limit 10"
    ).values()
    add_graph_data(nodes, edges, g_res)

    g = nx.DiGraph()
    for node in nodes:
        g.add_node(node[0], node_type=node[1])
    for edge in edges:
        g.add_edge(edge[0], edge[1], edge_type=edge[2])

    nx.write_graphml(
        g, pathlib.Path(__file__).parent.absolute() / "example_graph.graphml")
