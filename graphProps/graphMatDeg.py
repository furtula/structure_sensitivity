# -*- coding: utf-8 -*-
"""
===================
Topological indices
===================

Napisati dalje!!!!
"""

import numpy as np
import networkx as nx


def __numpymatrix(func):
    def to__numpymatrix(arg):
        if isinstance(arg, list):
            for row in arg:
                if not isinstance(row, list):
                    raise TypeError("You need to provide list of lists!")
                elif len(row) != len(arg):
                    raise TypeError("Matrix must be square one!")
            arg = np.matrix(arg)
        elif isinstance(arg, (str, int, float, dict)):
            raise TypeError('This is not a matrix!')
        return func(arg)
    return to__numpymatrix


@__numpymatrix
def distance_matrix(A):
    """Returns distance matrix of a graph G.

    Requrements
    -----------

    **numpy**

    Parameters
    ----------

    A : numpy matrix or list of lists
        This is n x n adjacency matrix of a graph.

    Returns
    -------

    dist : numpy matrix
        This is n x n distance matrix of a graph.

    Examples
    --------

    >>> a = [[0,1,0,0,0],
             [1,0,1,0,0],
             [0,1,0,1,0],
             [0,0,1,0,1],
             [0,0,0,1,0]]
    >>> print(distance_matrix(a))
    [[ 0.  1.  2.  3.  4.]
     [ 1.  0.  1.  2.  3.]
     [ 2.  1.  0.  1.  2.]
     [ 3.  2.  1.  0.  1.]
     [ 4.  3.  2.  1.  0.]]

    >>> import networkx as nx
    >>> g = nx.Graph()
    >>> g.add_edges_from([(0,1), (1,2), (2,3), (3,4), (4,5), (5,0), (5,6),
    (6,7), (7,8), (8,9), (9,0), (1,10), (6,11)])
    >>> a = nx.adj_matrix(g)
    >>> print(distance_matrix(a.toarray()))
    [[ 0.  1.  2.  3.  2.  1.  2.  3.  2.  1.  2.  3.]
     [ 1.  0.  1.  2.  3.  2.  3.  4.  3.  2.  1.  4.]
     [ 2.  1.  0.  1.  2.  3.  4.  5.  4.  3.  2.  5.]
     [ 3.  2.  1.  0.  1.  2.  3.  4.  5.  4.  3.  4.]
     [ 2.  3.  2.  1.  0.  1.  2.  3.  4.  3.  4.  3.]
     [ 1.  2.  3.  2.  1.  0.  1.  2.  3.  2.  3.  2.]
     [ 2.  3.  4.  3.  2.  1.  0.  1.  2.  3.  4.  1.]
     [ 3.  4.  5.  4.  3.  2.  1.  0.  1.  2.  5.  2.]
     [ 2.  3.  4.  5.  4.  3.  2.  1.  0.  1.  4.  3.]
     [ 1.  2.  3.  4.  3.  2.  3.  2.  1.  0.  3.  4.]
     [ 2.  1.  2.  3.  4.  3.  4.  5.  4.  3.  0.  5.]
     [ 3.  4.  5.  4.  3.  2.  1.  2.  3.  4.  5.  0.]]

    Notes
    -----

    Generating distance matrix from an adjacency matrix of a graph using
    Floyd--Warshall algorithm [1]_.

    References
    ----------

    .. [1] \O. Ivanciuc, Representing two--dimensional (2D) chemical structures
    with molecular graphs, in: J. L. Faulon, A. Bender, *Handbook of
    Chemoinformatics Algorithms*, CRC Press, Boca Raton, 2010, pp. 1--36.
    """
    infinity = np.inf
    dist = np.zeros((len(A), len(A)))

    # Preparing matrix dist for Floyd-Warshall algorithm
    for i in range(len(A) - 1):
        for j in range(i + 1, len(A)):
            if A[i, j] == 0:
                dist[i, j] = dist[j, i] = infinity
            else:
                dist[i, j] = dist[j, i] = A[i, j]

    # Floyd-Warshall algorithm
    for k in range(len(A)):
        for i in range(len(A)):
            for j in range(len(A)):
                dist[i, j] = min(dist[i, j],
                                 dist[i, k] + dist[k, j])
    return dist


@__numpymatrix
def harary_matrix(A):
    """
    """
    dist = distance_matrix(A).tolist()
    hmatrix = [[0] * len(dist)  for _ in range(len(dist))]
    for i in range(len(dist) - 1):
        for j in range(i + 1, len(dist)):
            hmatrix[i][j] = hmatrix[j][i] = 1 / dist[i][j]
    return np.matrix(hmatrix)


@__numpymatrix
def degree(A):
    """Returns list of vertex degrees of a graph.

    Requrements
    -----------

    **numpy**

    Parameters
    ----------

    A : numpy matrix or list of lists
        This is n x n adjacency matrix of a graph.

    Returns
    -------

    deg : list
        This is list of vertex degrees of a graph.

    Examples
    --------

    >>> a = [
                [0, 1, 0, 0, 1],
                [1, 0, 1, 0, 0],
                [0, 1, 0, 1, 0],
                [0, 0, 1, 0, 1],
                [1, 0, 0, 1, 0]
            ]
    >>> print(degree(a))
    [2, 2, 2, 2, 2]

    >>> import networkx as nx
    >>> g = nx.Graph()
    >>> g.add_edges_from([(0,1), (1,2), (2,3), (3,4), (4,5), (5,0),
    (5,6), (6,7), (7,8), (8,9), (9,0), (1,10), (6,11)])
    >>> a = nx.adj_matrix(g)
    >>> print(degree(a.toarray()))
    [3, 3, 2, 2, 2, 3, 3, 2, 2, 2, 1, 1]
    """
    deg = list()
    for row in A.tolist():
        deg.append(row.count(1))
    return deg


@__numpymatrix
def laplacian_matrix(A):
    """Returns Laplacian matrix of a graph G.

    Requrements
    -----------

    **numpy**

    Parameters
    ----------

    A : numpy matrix or list of lists
        This is n x n adjacency matrix of a graph.

    Returns
    -------

    lapl : numpy matrix
        This is n x n Laplacian matrix of a graph.

    Examples
    --------

    >>> a = [
                [0, 1, 0, 0, 1],
                [1, 0, 1, 0, 0],
                [0, 1, 0, 1, 0],
                [0, 0, 1, 0, 1],
                [1, 0, 0, 1, 0]
            ]
    >>> print(laplacian_matrix(a))
    [[ 2. -1.  0.  0. -1.]
     [-1.  2. -1.  0.  0.]
     [ 0. -1.  2. -1.  0.]
     [ 0.  0. -1.  2. -1.]
     [-1.  0.  0. -1.  2.]

    >>> import networkx as nx
    >>> g = nx.Graph()
    >>> g.add_edges_from([(0,1), (1,2), (2,3), (3,4), (4,5), (5,0),
    (5,6), (6,7), (7,8), (8,9), (9,0), (1,10), (6,11)])
    >>> a = nx.adj_matrix(g)
    >>> print(laplacian_matrix(a.toarray()))
    [[ 3. -1.  0.  0.  0. -1.  0.  0.  0. -1.  0.  0.]
     [-1.  3. -1.  0.  0.  0.  0.  0.  0.  0. -1.  0.]
     [ 0. -1.  2. -1.  0.  0.  0.  0.  0.  0.  0.  0.]
     [ 0.  0. -1.  2. -1.  0.  0.  0.  0.  0.  0.  0.]
     [ 0.  0.  0. -1.  2. -1.  0.  0.  0.  0.  0.  0.]
     [-1.  0.  0.  0. -1.  3. -1.  0.  0.  0.  0.  0.]
     [ 0.  0.  0.  0.  0. -1.  3. -1.  0.  0.  0. -1.]
     [ 0.  0.  0.  0.  0.  0. -1.  2. -1.  0.  0.  0.]
     [ 0.  0.  0.  0.  0.  0.  0. -1.  2. -1.  0.  0.]
     [-1.  0.  0.  0.  0.  0.  0.  0. -1.  2.  0.  0.]
     [ 0. -1.  0.  0.  0.  0.  0.  0.  0.  0.  1.  0.]
     [ 0.  0.  0.  0.  0.  0. -1.  0.  0.  0.  0.  1.]]
    """
    deg = np.diag(degree(A))
    return deg - A


@__numpymatrix
def signless_laplacian_matrix(A):
    """Returns signless Laplacian matrix of a graph G.

    Requrements
    -----------

    **numpy**

    Parameters
    ----------

    A : numpy matrix or list of lists
        This is n x n adjacency matrix of a graph.

    Returns
    -------

    lapl : numpy matrix
        This is n x n signless Laplacian matrix of a graph.
    """
    deg = np.diag(degree(A))
    return deg + A


@__numpymatrix
def randic_matrix(A):
    """
    Requirements
    ------------

    **numpy**

    Parameters
    ----------
    A : numpy matrix or list of lists
        This is n x n adjacency matrix of a graph.

    Returns
    -------
    RM : numpy matrix
        This is n x n Randic matrix of a graph.

    Examples
    --------
    >>> a = [[0,1,0,0,0], [1,0,1,1,0], [0,1,0,0,0], [0,1,0,0,1], [0,0,0,1,0]]
    >>> print(randic_matrix(a))
    [[ 0.          0.57735027  0.          0.          0.        ]
     [ 0.57735027  0.          0.57735027  0.40824829  0.        ]
     [ 0.          0.57735027  0.          0.          0.        ]
     [ 0.          0.40824829  0.          0.          0.70710678]
     [ 0.          0.          0.          0.70710678  0.        ]]
    """
    deg = degree(A)
    degstep = [i ** -.5   for i in deg]
    degstep = np.diag(degstep)
    RM = (degstep @ A) @ degstep
    return RM

@__numpymatrix
def degree_extended_matrix(A):
    """

    """
    degs = degree(A)
    xm = np.zeros((len(A), len(A)))
    for i in range(len(A)-1):
        for j in range(i+1, len(A)):
            if A[i, j] > 0:
                xm[i, j] = xm[j, i] = (
                    degs[i] / degs[j] + degs[j] / degs[i]) / 2
    return xm

@__numpymatrix
def detour_matrix(A):
    """Detour matrix - a matrix of largest distances between pairs of vertices.

    Parameters
    ==========

    A: numpy.matrix
        Adjacency matrix as numpy matrix.

    Returns
    =======

    detour: numpy.matrix
        Matrix of maximal lengths between pairs of vertices.
    """
    nxgraph = nx.from_numpy_matrix(A)
    n = nxgraph.number_of_nodes()
    detour = [[0] * n  for _ in range(n)]
    for i in range(n - 1):
        for j in range(i+1, n):
            pp = nx.all_simple_paths(nxgraph, i, j)
            detour[i][j] = detour[j][i] = max([len(l) - 1  for l in pp])
    return np.matrix(detour)


@__numpymatrix
def resistance_distance_matrix(A):
    """

    :param A:
    :return:
    """
    def entry(mat, x, y):
        return mat[x, x] + mat[y, y] - 2 * mat[x, y]

    lmatrix = laplacian_matrix(A)
    invlmat = np.linalg.pinv(lmatrix)

    resmat = np.zeros((len(invlmat), len(invlmat)))

    for i in range(len(invlmat)):
        for j in range(len(invlmat)):
            resmat[i, j] = entry(invlmat, i, j)
    return resmat


def matrix_spectrum(mat):
    """
    Requirements
    ------------

    **numpy**

    Parameters
    ----------

    mat : numpy matrix or list of lists
        Square symmetric matrix n x n with real entries.

    Returns
    -------

    matspectrum : numpy array
        Eigenvalue of square symmetric matrix.

    Examples
    --------

    >>> a = [[0,1,0,0,0], [1,0,1,1,0], [0,1,0,0,0], [0,1,0,0,1], [0,0,0,1,0]]
    >>> print(matrix_spectrum(a))
    [ -1.84775907e+00  -7.65366865e-01   1.32825242e-16   7.65366865e-01
   1.84775907e+00]
    """
    matspectrum = np.linalg.eigvalsh(mat)
    return matspectrum

@__numpymatrix
def number_of_vertices(A):
    return len(A)

@__numpymatrix
def number_of_edges(A):
    m = 0
    for i in range(len(A)-1):
        for j in range(i+1, len(A)):
            if A[i, j] > 0:
                m += 1
    return m
