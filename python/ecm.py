'''Entropy c-Means Clustering

Paper Source:
Gupta A., Datta S. and Das S., "Fuzzy Clustering to Identify Clusters at Different Levels of Fuzziness: An Evolutionary Multiobjective Optimization Approach," in IEEE Transactions on Cybernetics, 2019.
'''

# Author: Avisek Gupta


import numpy as np
import pygmo as pg
from scipy.spatial.distance import cdist


class ecm:
    '''Definition of the multi-objective problem for Entropy c-Means

    Usage
    -----

    Set up the problem as:
    spobj = ecm(k * data.shape[1])
    spobj.set_data(data)
    prob = pg.problem(spobj)

    where k is the number of clusters,
    and data.shape[1] is the number of features
    data is the data array.
    '''

    def __init__(self, dim):
        self.dim = dim

    def fitness(self, x):
        centers = np.reshape(x, (x.shape[0] // self.data.shape[1], self.data.shape[1]))
        if centers.ndim == 1:
            dist = np.fmax(
                cdist(centers[:,None], self.data, metric='sqeuclidean'), 1e-15
            )
        else:
            dist = np.fmax(
                cdist(centers, self.data, metric='sqeuclidean'), 1e-15
            )
        u = np.exp(-dist / (2 * (self.sigma ** 2)))
        u = u / u.sum(axis=0)

        f1 = np.sum(u * dist)
        f2 = np.sum(- u * np.log(np.fmax(u, 1e-15)))

        return [f1, f2]

    def get_nobj(self):
        return 2

    def get_bounds(self):
        return ([-1] * self.dim, [+1] * self.dim)

    def get_name(self):
        return "Entropy c-Means"

    def set_data(self, data):
        self.data = data
        self.sigma = (((data - data.mean(axis=0))**2).sum(axis=1)).std()


def ecm_nsga2(
    data, k=2, popsize=40, gen=500, cr=0.95, eta_c=10, m=0.01, eta_m=10
):
    '''Entropy c-means - NSGA-II Clustering

    Parameters
    ----------

    data : array, shape (n_data_points, n_features)
        The data array.

    k : int, default: 2
    The number of clusters.

    popsize : int, default: 40
        The population size.

    gen : int, default: 500
        The number of generations to be iterated.

    cr : float, default: 0.95
        NSGA-II parameter for crossover probability.

    eta_c : float, default: 10
        NSGA-II parameter for crossover distribution index.

    m : float, default: 0.01
        NSGA-II parameter for mutation probability.

    eta_m : float, default: 10
        NSGA-II parameter for mutation distribution index.

    Returns
    -------

    vectors: array, shape (popsize, k * n_features)
        The resulting cluster centers, flattened to 1 dimensional arrays.

    pareto_front: array, shape (popsize, 2)
        The pareto front of mapped solutions

    '''

    # set up the problem
    spobj = ecm(k * data.shape[1])
    spobj.set_data(data)
    prob = pg.problem(spobj)

    # create population
    pop = pg.population(prob, popsize)
    # select the MO algorithm
    algo = pg.algorithm(pg.nsga2(
        gen=gen, cr=cr, eta_c=eta_c, m=m, eta_m=eta_m,
    ))
    # run optimization
    pop = algo.evolve(pop)

    # extract results
    pareto_front, vectors = pop.get_f(), pop.get_x()
    # sort the Pareto front and solutions
    sorted_idxs = np.argsort(pareto_front[:, 0])
    pareto_front = pareto_front[sorted_idxs, :]
    vectors = vectors[sorted_idxs, :]

    return vectors, pareto_front


def ecm_moead(
    data, k=2, popsize=40, gen=500, weight_generation="grid",
    decomposition="tchebycheff", neighbours=20, CR=1, F=0.5, eta_m=20,
    realb=0.9, limit=2, preserve_diversity=True
):
    '''Entropy c-means - MOEA/D Clustering

    Parameters
    ----------

    data : array, shape (n_data_points, n_features)
        The data array.

    k : int, default: 2
    The number of clusters.

    popsize : int, default: 40
        The population size.

    gen : int, default: 500
        The number of generations to be iterated.

    weight_generation : str, default: 'grid'
        Method used for weight generation.

    decomposition : str, default: 'tchebycheff'
        Method for object decomposition.

    neighbours : int, default: 20
        Size of weight neighbourhood.

    CR : float, default: 1
        MOEA/D parameter for Differential Evolution crossover parameter.

    F : float, default: 0.5
        MOEA/D parameter for Differential Evolution operator.

    eta_m : float, default: 20
        MOEA/D parameter for polynomial mutation distribution index.

    realb : float, default: 0.9
        Chance that the neighbourhood is considered at each generation,
        rather than the whole population (only if preserve_diversity is true).

    limit : int, default: 2
        Maximum number of copies reinserted in the population
        (only if m_preserve_diversity is true).

    preserve_diversity : bool, default: True
        When true activates diversity preservation mechanisms.

    Returns
    -------

    vectors: array, shape (popsize, k * n_features)
        The resulting cluster centers, flattened to 1 dimensional arrays.

    pareto_front: array, shape (popsize, 2)
        The pareto front of mapped solutions

    '''
    # set up the problem
    spobj = ecm(k * data.shape[1])
    spobj.set_data(data)
    prob = pg.problem(spobj)

    # create population
    pop = pg.population(prob, popsize)
    # select the MO algorithm
    algo = pg.algorithm(pg.moead(
        gen=gen, weight_generation=weight_generation,
        decomposition=decomposition, neighbours=neighbours,
        CR=CR, F=F, eta_m=eta_m, realb=realb, limit=limit,
        preserve_diversity=preserve_diversity
    ))
    # run optimization
    pop = algo.evolve(pop)

    # extract results
    pareto_front, vectors = pop.get_f(), pop.get_x()
    # sort the Pareto front and solutions
    sorted_idxs = np.argsort(pareto_front[:, 0])
    pareto_front = pareto_front[sorted_idxs, :]
    vectors = vectors[sorted_idxs, :]

    return vectors, pareto_front


if __name__ == '__main__':
    # DEBUG
    '''
    from sklearn.datasets import load_iris
    X = load_iris().data
    X = (
        (X - X.min(axis=0)) / np.fmax(X.max(axis=0) - X.min(axis=0), 1e-15)
    ) * 2 - 1

    k = 3
    '''

    '''
    vectors, pareto_front = ecm_nsga2(
        X, k=k, popsize=60, gen=1000, cr=0.9, eta_c=20, m=0.1, eta_m=10
    )
    '''

    '''
    vectors, pareto_front = ecm_moead(
        X, k=3, popsize=60, gen=1000, weight_generation="grid",
        decomposition="tchebycheff", neighbours=40, CR=0.6, F=0.3, eta_m=10,
        realb=0.9, limit=2, preserve_diversity=True
    )
    '''

    '''
    from sklearn.metrics import adjusted_rand_score as ARI
    ari_max = -1
    for i in range(vectors.shape[0]):
        ecm_centers1 = vectors[i, :].reshape(
            vectors[i, :].shape[0] // X.shape[1], X.shape[1]
        )
        if ecm_centers1.ndim == 1:
            dist = cdist(ecm_centers1[:, None], X, metric='sqeuclidean')
        else:
            dist = cdist(ecm_centers1, X, metric='sqeuclidean')
        sigma = (((X - X.mean(axis=0))**2).sum(axis=1)).std()
        u = np.exp(-dist / (2 * (sigma ** 2)))
        u = u / u.sum(axis=0)

        #import matplotlib.pyplot as plt
        #plt.scatter(X[:,2], X[:,3], marker='x', c='gray')
        #plt.scatter(ecm_centers1[:,2], ecm_centers1[:,3], marker='o', c='r')
        #plt.show()

        ari = ARI(load_iris().target, u.argmax(axis=0))
        if ari_max < ari:
            ari_max = ari
        print('soln',i,'ARI:', ari)
    print('Max ARI =', ari_max)
    '''
