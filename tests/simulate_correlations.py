import numpy as np
from scipy.stats import random_correlation


np.random.seed(42)


def simulate_correlated_data(n_features: int, max_: int = 25) -> np.ndarray:
    eigs = np.random.uniform(1, max_, n_features) / \
            np.abs(np.random.randn(n_features))
    eigs /= (np.sum(eigs) / n_features)
    return random_correlation.rvs(eigs)


if __name__ == "__main__":
    import pandas as pd
    import networkx as nx
    import matplotlib.pyplot as plt

    m = 25
    g = nx.barabasi_albert_graph(m, 15)
    corr = simulate_correlated_data(m)
    corr_adj = corr * nx.to_numpy_array(g)

    fig, ax = plt.subplots(figsize=(16, 9), ncols=2)
    ax[0].hist(np.triu(corr, k=-1))
    ax[1].hist(np.triu(corr_adj, k=-1))
    plt.show()
    plt.close(fig)

    # adapting to test graph
    metabos = pd.read_csv("test_data/test_metabolites.csv").iloc[:, 0]
    sim_corr = simulate_correlated_data(metabos.size)
    pd.DataFrame(
        sim_corr, index=metabos.values, columns=metabos.values
    ).to_csv("test_data/simulated_metabolite_correlations.csv")
