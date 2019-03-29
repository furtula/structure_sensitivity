"""Informacioni indeksi zasnovani na entropijama."""

from graphProps import distance_matrix
from moldes.distanceBased import w, funcd
import numpy as np


def distance_partition(adjmatrix):
    """Distance partition.

    Parametri
    ---------
        adjmatrix : numpy.matrix
            Matrica susedstva grafa G.

    Rezultati
    ---------
        k : int
            Broj javljanja rastojanja k u matrici rastojanja.

    """
    dm = distance_matrix(adjmatrix)
    k = dict()
    for i in range(len(dm) - 1):
        for j in range(i+1, len(dm)):
            k[int(dm[i, j])] = k.get(int(dm[i, j]), 0) + 1
    return k


def bt_entropy(adjmatrix):
    """Bonchev Trinajstic entropija.

    Parametri
    ---------
        adjmatrix : numpy.matrix
            Matrica susedstva grafa G.

    Rezultati
    ---------
        idg, idgm : float
            Entropijska mera zasnovana na metrickim osobinama grafa.
            Bonchev i Trinajstic.

    """
    k = distance_partition(adjmatrix)
    v = len(adjmatrix)
    idg = 0
    idgm = 0
    for i in k.keys():
        idg += v * v * np.log2(v * v) - v * np.log2(v) \
            - 2 * k[i] * np.log2(2 * k[i])
        idgm += (-1) * np.log2(1 / v) / v \
            - 2 * k[i] / (v * v) * np.log2(2 * k[i] / (v * v))
    return idg, idgm


def br_entropy(adjmatrix):
    """Bonchev-Rovray entropija.

    Parametri
    ---------
        adjmatrix : numpy.matrix
            Matrica susedstva grafa G.

    Rezultati
    ---------,
        iwg, iwgm : float
            Bonchev i Rouvray indeks.

    """
    k = distance_partition(adjmatrix)
    iwg = 0
    iwgm = 0
    wiener = funcd(adjmatrix, w)
    for i in k.keys():
        iwg += wiener * np.log2(wiener) - i * k[i] * np.log2(i)
        iwgm += (-1) * i * k[i] * np.log2(i / wiener) / wiener
    return iwg, iwgm
