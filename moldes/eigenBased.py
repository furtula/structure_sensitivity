"""Eigenvalue-based graph invariants."""
import graphProps as gp
import numpy as np


def energy(M, mat='adjacency', c=0):
    """Some kind of graph energy.

    Parametri
    ---------
        M : numpy.matrix
            Matrica iz koje se racuna spektar grafa i potom
            energija.

        mat : str
            Nacin na koji biramo matricu.
            Podrazumevana vrednost je 'adjacency'.

        c : float
            Konstanta od koje se oduzima sopstvena vrednost grafa.
            Podrazumavana vrednost je 0.
            Npr. kod laplasove energije c = 2*m/n .

    Rezultati
    ---------
        e : float
            Neka energija grafa.
    """
    if mat == 'laplacian':
        M = gp.laplacian_matrix(M)
    elif mat == "distance":
        M = gp.distance_matrix(M)
    elif mat == 'randic':
        M = gp.randic_matrix(M)
    elif mat == 'xu':
        M = gp.degree_extended_matrix(M)
    elif mat == 'harary':
        M = gp.harary_matrix(M)
    else:
        pass

    e = 0
    spec = gp.matrix_spectrum(M)
    for l in spec:
        e += abs(l - c)
    return e


def lel(adjmatrix):
    """Laplacian energy like.

    Parametri
    ---------
        laplacian_matrix : numpy.matrix
            Laplasova matrica.

    Rezultati
    ---------
        lel :


    """
    lel = 0
    lspec = gp.matrix_spectrum(gp.laplacian_matrix(adjmatrix))
    for l in lspec:
        if l < 1e-8:
            continue
        else:
            lel += l ** .5
    return lel


def estrada(M, mat='adjacency', c=0):
    """Some kind of the Estrada index.

    Parametri
    ---------
        M : numpy.matrix
            Matrica iz koje se racuna spektar grafa i potom
            energija.

        mat : str
            Nacin na koji biramo matricu.
            Podrazumevana vrednost je 'adjacency'.

        c : float
            Konstanta od koje se oduzima sopstvena vrednost grafa.
            Podrazumavana vrednost je 0.
            Npr. kod laplasove energije c = 2*m/n .

    Rezultati
    ---------
        ee : float
            Neki Estradin indeks.
    """
    if mat == 'laplacian':
        M = gp.laplacian_matrix(M)
    elif mat == "distance":
        M = gp.distance_matrix(M)
    elif mat == 'randic':
        M = gp.randic_matrix(M)
    elif mat == 'xu':
        M = gp.degree_extended_matrix(M)
    elif mat == 'harary':
        M = gp.harary_matrix(M)
    else:
        pass

    ee = 0
    spec = gp.matrix_spectrum(M)
    for l in spec:
        ee += np.exp(l - c)
    return ee


def resolvent_energy(adjmatrix):
    """Resolvent energy of a graph.

    Parametri
    ---------
        adjmatrix : numpy.matrix
            Matrica susedstva.

    Rezultati
    ---------
        re : float
            Resolventna energija.
    """
    re = 0
    n = len(adjmatrix)
    aspec = gp.matrix_spectrum(adjmatrix)
    for l in aspec:
        re += 1 / (n - l)
    return re


def kirchhoff(adjmatrix):
    """

    :param adjmatrix:
    :return:
    """
    lmatrix = gp.laplacian_matrix(adjmatrix)
    lspec = gp.matrix_spectrum(lmatrix)
    kf = 0
    for l in lspec:
        if l <= 1e-8:
            continue
        else:
            kf += len(lspec) / l
    return kf


def homolumo(adjmatrix):
    """

    :param adjmatrix:
    :return:
    """
    spec = gp.matrix_spectrum(adjmatrix)
    h = None
    l = None
    for i in spec:
        if i > 0:
            if h == None or h > i:
                h = i
        elif i < 0:
            if l == None or l < i:
                l = i
    return h - l

def incidence_energy(adjmatrix):
    """

    :param adjmatrix:
    :return:
    """
    qmatrix = gp.signless_laplacian_matrix(adjmatrix)
    qspec = gp.matrix_spectrum(qmatrix)
    ie = 0
    for i in qspec:
        if abs(i) >= 1e-8:
            ie += np.sqrt(i)
    return ie
