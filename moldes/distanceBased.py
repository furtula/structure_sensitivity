"""Distance-based graph invariants."""
import numpy as np

from graphProps import (degree, detour_matrix, distance_matrix,
                        resistance_distance_matrix)

# distance_matrix = gm.distance_matrix


def n1n2(adjmatrix):
    """Pravi listu sa vrednostima za n1 i n2.

    Parametri
    ---------
        adjmatrix : numpy.matrix
            Matrica susedstva grafa G.

    Rezultati
    ---------
        n1n2 : int
            Proizvod broja cvorova sa jedne i druge strane
            zadate grane.
    """
    dmatrix = distance_matrix(adjmatrix)
    n1n2 = []
    for u in range(len(dmatrix) - 1):
        for v in range(u + 1, len(dmatrix)):
            if dmatrix[u, v] == 1:
                n1 = 0
                n2 = 0
                for p in range(len(dmatrix)):
                    if dmatrix[u, p] < dmatrix[v, p]:
                        n1 += 1
                    elif dmatrix[u, p] > dmatrix[v, p]:
                        n2 += 1
                n1n2.append((n1, n2))
    return n1n2


def diameter(adjmatrix):
    """Izracunava diametar grafa (maksimalno rastojanje).

    Parametri
    ---------

        adjmatrix : numpy.matrix
            Matrica susedstva grafa G.

    Rezultati
    ---------

        dia : int
            Maksimalno rastojanje u matrici rastojanja.
    """
    dist = distance_matrix(adjmatrix)
    dia = None
    for row in dist:
        if dia is None or dia < max(row):
            dia = max(row)
    return dia


def w(distance, x=1):
    """Vraca rastojanje na neki stepen.

    Parametri
    ---------

        distance : int
            Rastojanje izmedju dva cvora.

        x : float
            Stepen na koji se dize rastojanje.

    Rezultati
    ---------

        distance ** x : float
            Rastojanje na x-ti stepen.
    """
    return distance**x


def harary(distance, x=1):
    return 1 / distance


def ww(distance, x=1):
    """Vraca elemenat za izracunavanje hiper-Viner indeksa.

    Parametri
    ---------

        distance : int
            Rastojanje izmedju dva cvora.

        x : float
            Stepen na koji se dize rastojanje.

    Rezultati
    ---------

        ((distance + distance * distance) / 2) ** x : float
            HiperWienner-ov indeks.
    """
    return ((distance + distance * distance) / 2)**x


def tsz(distance, x=1):
    r"""Vraca elemenat za izracunavanje Trac-Stankivec-Zefirov indeksa.

    Parametri
    ---------

        distance : int
            Rastojanje izmedju dva cvora.

        x : float
            Stepen na koji se dize rastojanje.

    Rezultati
    ---------

        ((2 * distance + 3 * distance * distance + distance ** 3) \
        / 6) ** x : float
            Trac-Stankevic-Zefir-ov indeks.
    """
    return ((2 * distance + 3 * distance * distance + distance**3) / 6)**x


class DistanceTI(object):
    """Distance based topological index."""

    def __init__(self, funkcija, ime, x=1, choose='r'):
        self.func = funkcija
        self.x = x
        self.choose = choose
        self._ime = ime

    def ti(self, adjmatrix):
        dist = distance_matrix(adjmatrix)
        deg = degree(adjmatrix)
        y = 0
        if self.choose == 't':
            for i in range(len(dist) - 1):
                for j in range(i + 1, len(dist)):
                    if deg[i] == 1 and deg[j] == 1:
                        y += self.func(dist[i, j], self.x)
            return y
        elif self.choose == 's':
            for i in range(len(dist) - 1):
                for j in range(i + 1, len(dist)):
                    if deg[i] == 1 or deg[j] == 1:
                        y += self.func(dist[i, j], self.x)
            return y
        else:
            for i in range(len(dist) - 1):
                for j in range(i + 1, len(dist)):
                    y += self.func(dist[i, j], self.x)
        return y

    def __str__(self):
        return self._ime

# def funcd(adjmatrix, funkcija, x=1, choose='r'):
#     """Funkcija koja racuna indekse pomocu gornjih funkija.
#
#     Parametri
#     ---------
#
#         adjmatrix : numpy.matrix
#             Matrica susedstva grafa G.
#
#         funkcija : function
#             Funkcija koja racuna terminalni, semi-terminalni
#             ili obican(regularan) indeks.
#
#         x : float
#             Varijabilni parametar.
#
#         choose : string
#             Bira se terminalni, semi-terminalni
#             ili obican(regularan) indeks.
#
#     Rezultati
#     ---------
#
#         y : float
#             Regularan, terminalni ili semi-terminalni indeks.
#     """
#     dist = distance_matrix(adjmatrix)
#     deg = degree(adjmatrix)
#     y = 0
#
#     if choose == 't':
#         for i in range(len(dist) - 1):
#             for j in range(i + 1, len(dist)):
#                 if deg[i] == 1 and deg[j] == 1:
#                     y += funkcija(dist[i, j], x)
#         return y
#     elif choose == 's':
#         for i in range(len(dist) - 1):
#             for j in range(i + 1, len(dist)):
#                 if deg[i] == 1 or deg[j] == 1:
#                     y += funkcija(dist[i, j], x)
#         return y
#     else:
#         for i in range(len(dist) - 1):
#             for j in range(i + 1, len(dist)):
#                 y += funkcija(dist[i, j], x)
#     return y


def szeged(adjmatrix, x=1):
    """
    Parametri
    ---------
        adjmatrix : numpy.matrix
            Matrica susedstva grafa G.

        x : float
            Varijabilni parametar.

    Rezultati
    ---------
        sz : float
            Szegedinski indeks.
    """
    parovi = n1n2(adjmatrix)
    sz = 0
    for n1, n2 in parovi:
        sz += (n1 * n2)**x
    return sz


def pi(adjmatrix, x=1):
    """
    Parametri
    ---------

        adjmatrix : numpy.matrix
            Matrica susedstva grafa G.

        x : float
            Varijabilni parametar.

    Rezultati
    ---------
        pi : float
            Padmakar-Ivan-ov indeks.


    """
    parovi = n1n2(adjmatrix)
    pi = 0
    for n1, n2 in parovi:
        pi += (n1 + n2)**x

    return pi


def second_geometric_arithmetic(adjmatrix):
    """
    Parametri
    ---------
        adjmatrix : numpy.matrix
            Matrica susedstva grafa G.

    Rezultati
    ---------
        ga : float
            Second Geometric-Arithmetic-indeks.

    """
    parovi = n1n2(adjmatrix)
    ga2 = 0
    for n1, n2 in parovi:
        ga2 += 2 * (n1 * n2)**.5 / (n1 + n2)
    return ga2


def graovac_ghorbani(adjmatrix):
    """
    Parametri
    ---------

        adjmatrix : numpy.matrix
            Matrica susedstva grafa G.

    Rezultati
    ---------

        gg : float
            Graovac-Ghorbani-indeks.
    """
    parovi = n1n2(adjmatrix)
    gg = 0
    for n1, n2 in parovi:
        gg += ((n1 + n2 - 2) / (n1 * n2))**.5

    return gg


def reciprocal_complementary_w(adjmatrix, x=1):
    """
    Parametri
    ---------

        adjmatrix : numpy.matrix
            Matrica susedstva grafa G.

        x : float
            Varijabilni parametar.

    Rezultati
    ---------

        rcw : float
            Reciprocal-Complementary-Wienner indeks.

    """
    dist = distance_matrix(adjmatrix)
    dia = diameter(adjmatrix)
    rcw = 0
    for i in range(len(dist) - 1):
        for j in range(i + 1, len(dist)):
            rcw += 1 / (dia + 1 - dist[i, j])**x
    return rcw


def balaban(adjmatrix, sum_balaban=False):
    '''
    Parametri
    ---------

        adjmatrix : numpy matrix ili list of lists
            Matrica susedstva grafa G.

        sum_balaban : bool
            Uslov koji ukoliko je ispunjen, funkcija racuna
            sum-Balabanov indeks. Ukoliko uslov nije ispunjen
            racuna se originalni Balabanov indeks.

    Rezultati
    ---------

        jj : float
            Balabanov ili sum-Balabanov indeks u zavisnosti od
            ispunjenosti uslova sum_balaban.
    '''
    n = len(adjmatrix)
    m = 0
    for row in adjmatrix.tolist():
        m += sum(row) / 2

    dm = distance_matrix(adjmatrix)

    c = m / (m - n + 2)

    jj = 0

    rowsum = []

    for row in dm:
        rowsum.append(sum(row))

    for i in range(len(dm) - 1):
        for j in range(i + 1, len(dm)):
            if dm[i, j] == 1:
                if not sum_balaban:
                    jj += 1 / np.sqrt(rowsum[i] * rowsum[j])
                else:
                    jj += 1 / np.sqrt(rowsum[i] + rowsum[j])
    return c * jj


def detour(adjmatrix):
    """
    :param adjmatrix:
    :return:
    """
    detmatrix = detour_matrix(adjmatrix)
    dd = 0
    for row in detmatrix.tolist():
        dd += sum(row) / 2
    return dd


def hyper_detour(adjmatrix):
    """

    :param adjmatrix:
    :return:
    """
    detmatrix = detour_matrix(adjmatrix)
    hdd = 0
    for i in range(len(detmatrix) - 1):
        for j in range(i + 1, len(detmatrix)):
            hdd += detmatrix[i, j] * (detmatrix[i, j] + 1)
    return hdd


def sdegree_kirchhoff(adjmatrix):
    """

    :param adjmatrix:
    :return:
    """
    rdmatrix = resistance_distance_matrix(adjmatrix).tolist()
    rdmrowsums = []
    for row in rdmatrix:
        rdmrowsums.append(sum(row))
    sdkf = 0
    for i in range(len(rdmrowsums) - 1):
        for j in range(i + 1, len(rdmrowsums)):
            if adjmatrix[i, j] > 0:
                sdkf += rdmrowsums[i] + rdmrowsums[j]
    return sdkf
