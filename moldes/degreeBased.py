"""Degree-based graph invariants."""
import numpy as np

from graphProps import degree


def zajednicka_promenljiva(adjmatrix, ind=True):
    """Auxillary function.

    Parametri
    ----------
    adjmatrix : numpy.matrix
        Matrica susedstva grafa G.

    ind : Bool
        Ukoliko je ind tacno, onda se racuna degree-based indeks.
        U suprotnom racuna se koindeks.

    Rezultati
    ---------
    lista_stepena : list
        Lista stepena parova susednih cvorova ukoliko je ind tacno.
        U suprotnom to je lista parova nesusednih cvorova.

    """
    deg = degree(adjmatrix)
    lista_stepena = []
    if ind:
        for i in range(len(adjmatrix) - 1):
            for j in range(i + 1, len(adjmatrix)):
                if adjmatrix[i, j] == 1:
                    lista_stepena.append((deg[i], deg[j]))
    else:
        for i in range(len(adjmatrix) - 1):
            for j in range(i + 1, len(adjmatrix)):
                if adjmatrix[i, j] != 1:
                    lista_stepena.append((deg[i], deg[j]))
    return 0, lista_stepena


def first_zagreb(adjmatrix, x=1, indeks=True):
    """Izracunava prvi zagrebacki indeks sa varijablinim parametrom x.

    Parametri
    ---------

    adjmatrix : numpy.matrix
        Matrica susedstva grafa G.

    x : float
        Varijabilni parametar.

    Rezultati
    ---------

    m1 : float
        Prvi varijabilni zagrebacki indeks.

    Primer
    ------

    >>> import numpy as np
    >>> a = np.matrix([[0,1,0], [1,0,1], [0,1,0]])
    >>> print(firstZagreb(a, 2))
    6
    """
    m1, dlist = zajednicka_promenljiva(adjmatrix, indeks)
    for d_u, d_v in dlist:
        m1 += d_u**x + d_v**x
    return m1


def sigma(adjmatrix, x=1, indeks=True):
    """Izracunava sigma indeks sa varijablinim parametrom x.

    Parametri
    ---------

    adjmatrix : numpy.matrix
        Matrica susedstva grafa G.

    x : float
        Varijabilni parametar.

    Rezultati
    ---------

    s : float
        Prvi varijabilni zagrebacki indeks.

    Primer
    ------

    >>> import numpy as np
    >>> a = np.matrix([[0,1,0], [1,0,1], [0,1,0]])
    >>> print(firstZagreb(a, 2))
    6
    """
    s, dlist = zajednicka_promenljiva(adjmatrix, indeks)
    for d_u, d_v in dlist:
        s += (d_u**x - d_v**x) ** 2
    return s


def first_product_zagreb(adjmatrix, x=1, indeks=True):
    """Izracunava prvi multiplikativni zagrebacki indeks sa 
    varijablinim parametrom x.

    Parametri
    ---------

    adjmatrix : numpy.matrix
        Matrica susedstva grafa G.

    x : float
        Varijabilni parametar.

    Rezultati
    ---------

    pm1 : float
        Prvi varijabilni multiplikativni zagrebacki indeks.
    """
    degs = []
    for row in adjmatrix:
        degs.append(sum(row))
    pm1 = 1
    for d_u in degs:
        pm1 *= d_u**2
    return pm1


def first_modified_product_zagreb(adjmatrix, x=1, indeks=True):
    """Izracunava prvi modifikovani multiplikativni zagrebacki 
    indeks sa varijablinim parametrom x.

    Parametri
    ---------

    adjmatrix : numpy.matrix
        Matrica susedstva grafa G.

    x : float
        Varijabilni parametar.

    Rezultati
    ---------

    mpm1 : float
        Prvi modifikovani multiplikativni varijabilni zagrebacki indeks.
    """
    mpm1, dlist = zajednicka_promenljiva(adjmatrix, indeks)
    mpm1 += 1
    for d_u, d_v in dlist:
        mpm1 *= d_u**x + d_v**x
    return mpm1


def second_zagreb(adjmatrix, x=1, indeks=True):
    """Izracunava drugi zagrebacki indeks sa varijablinim parametrom x.

    Parametri
    ---------

        adjmatrix : numpy.matrix
            Matrica susedstva grafa G.

        x : float
            Varijabilni parametar.

    Rezultati
    ---------

        m_2 : float
             Drugi varijabilni zagrebacki indeks.

    """
    m_2, dlist = zajednicka_promenljiva(adjmatrix, indeks)

    for d_u, d_v in dlist:
        m_2 += (d_u * d_v) ** x
    return m_2


def second_product_zagreb(adjmatrix, x=1, indeks=True):
    """Izracunava drugi multiplikativni zagrebacki indeks sa varijablinim 
    parametrom x.

    Parametri
    ---------

        adjmatrix : numpy.matrix
            Matrica susedstva grafa G.

        x : float
            Varijabilni parametar.

    Rezultati
    ---------

        mm2 : float
             Drugi varijabilni multiplikativni zagrebacki indeks.

    """
    mm2, dlist = zajednicka_promenljiva(adjmatrix, indeks)

    mm2 += 1

    for d_u, d_v in dlist:
        mm2 *= (d_u * d_v) ** x
    return mm2


def sum_connectivity(adjmatrix, x=1, indeks=True):
    """Izracunava sum connectivity indeks sa varijablinim parametrom x.

    Parametri
    ---------

        adjmatrix : numpy.matrix
            Matrica susedstva grafa G.

        x : float
            Varijabilni parametar.

    Rezultati
    ---------

        sci : float
             Drugi varijabilni zagrebacki indeks.

    """
    sci, dlist = zajednicka_promenljiva(adjmatrix, indeks)

    for d_u, d_v in dlist:
        sci += (d_u + d_v) ** x
    return sci


def narumi_katayama(adjmatrix, x=1, indeks=True):
    """Izracunava Narumi-Katayama indeks sa varijabilnim parametrom x.

    Parametri
    ---------

        adjmatrix : numpy.matrix
            Matrica susedstva grafa G.

        x : float
            Varijabilni parametar.

    Rezultati
    ---------

        lnk : float
              Narumi-Katajama indeks.
    """
    lnk, dlist = zajednicka_promenljiva(adjmatrix, indeks)

    for d_u, d_v in dlist:
        lnk += (np.log(d_u) / d_u)**x + (np.log(d_v) / d_v)**x
    return np.exp(lnk)


def atom_bond_connectivity(adjmatrix, x=.5, indeks=True):
    """Izracunava atom-bond connectivity indeks sa varijablinim parametrom x.

    Parametri
    ---------

        adjmatrix : numpy.matrix
            Matrica susedstva grafa G.

        x : float
            varijabilni parametar.

    Rezultati
    ---------

        abc : float
              Atom Bond Connectivity-indeks.
    """
    abc, dlist = zajednicka_promenljiva(adjmatrix, indeks)

    for d_u, d_v in dlist:
        abc += ((d_u + d_v - 2) / (d_u * d_v)) ** x
    return abc


def geometric_arithmetic(adjmatrix, indeks=True):
    """Izracunava geometric-arithmetic indeks.

    Parametri
    ---------

        adjmatrix : numpy.matrix
            Matrica susedstva grafa G.

    Rezultati
    ---------

        ga : float
            Geometric-Arithmetic-indeks.

    """
    ga, dlist = zajednicka_promenljiva(adjmatrix, indeks)

    for d_u, d_v in dlist:
        ga += 2 * np.sqrt(d_u * d_v) / (d_u + d_v)
    return ga


def harmonic(adjmatrix, x=1, indeks=True):
    """Izracunava harmonic indeks sa varijablinim parametrom x.

    Parametri
    ---------

        adjmatrix : numpy.matrix
            Matrica susedstva grafa G.

        x : float
            Varijabilni parametar.

    Rezultati
    ---------

        harmonic : float
            Harmonic-indeks.
    """
    harmonic, dlist = zajednicka_promenljiva(adjmatrix, indeks)

    for d_u, d_v in dlist:
        harmonic += 1 / (d_u + d_v) ** x
    return harmonic


def albertson(adjmatrix, indeks=True):
    """Izracunava albertson indeks sa varijablinim parametrom x.

    Parametri
    ---------

        adjmatrix : numpy.matrix
            Matrica susedstva grafa G.

    Rezultati
    ---------

        alb : float
            Albertson-indeks.
    """
    alb, dlist = zajednicka_promenljiva(adjmatrix, indeks)

    for d_u, d_v in dlist:
        alb += abs(d_u - d_v)
    return alb


def reduced_second_zagreb(adjmatrix, x=1, indeks=True):
    """Izracunava redukovani drugi zagrebacki indeks sa varijablinim 
    parametrom x.

    Parametri
    ---------

        adjmatrix : numpy.matrix
            Matrica susedstva grafa G.

        x : float
            Varijabilni parametar.

    Rezultati
    ---------

        rm2 : float
            Reduced-Second-Zagreb-indeks.
    """
    rm2, dlist = zajednicka_promenljiva(adjmatrix, indeks)

    for d_u, d_v in dlist:
        rm2 += ((d_u - 1) * (d_v - 1)) ** x
    return rm2


def sdd(adjmatrix, x=1, indeks=True):
    """Izracunava sdd indeks sa varijablinim parametrom x.

    Parametri
    ---------

        adjmatrix : numpy.matrix
            Matrica susedstva grafa G.

        x : float
            varijabilni parametar.

    Rezultati
    ---------

        y : float
            Yang-ov indeks (J. Chem. Comput. Sci. 34 (1994) 1140-1145).
    """
    y, dlist = zajednicka_promenljiva(adjmatrix, indeks)

    for d_u, d_v in dlist:
        y += (d_u / d_v + d_v / d_u) ** x
    return y


def inverse_sum_indeg(adjmatrix, x=1, indeks=True):
    """Izracunava inverse sum indeg indeks sa varijablinim parametrom x.

    Parametri
    ---------

        adjmatrix : numpy.matrix
            Matrica susedstva grafa G.

        x : float
            varijabilni parametar.

        indeks : Bool
            Odlucuje da li je indeks ili koindeks.

    Rezultati
    ---------

        isi : float
            Inverse sum indeg indeks.
            (Croat. Chem. Acta 83 (2010) 243-260.).
    """
    isi, dlist = zajednicka_promenljiva(adjmatrix, indeks)

    for d_u, d_v in dlist:
        isi += (d_u * d_v / (d_u + d_v)) ** x
    return isi
