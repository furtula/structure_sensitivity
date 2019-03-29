import graphProps as gp


def eccentricity(adjmatrix):
    """

    :param adjmatrix:
    :return:
    """
    dmatrix = gp.distance_matrix(adjmatrix).tolist()
    ecc = list()
    for row in dmatrix:
        ecc.append(max(row))
    return ecc


def vertex_distance(adjmatrix):
    """

    :param adjmatrix:
    :return:
    """
    dmatrix = gp.distance_matrix(adjmatrix).tolist()
    vdist = list()
    for row in dmatrix:
        vdist.append(sum(row))
    return vdist


def eccentric_connectivity(adjmatrix):
    """

    :param adjmatrix:
    :return:
    """
    degs = gp.degree(adjmatrix)
    ecc = eccentricity(adjmatrix)
    eccin = 0
    for i in range(len(degs)):
        eccin += degs[i] * ecc[i]
    return eccin


def degree_distance(adjmatrix):
    """

    :param adjmatrix:
    :return:
    """
    vdist = vertex_distance(adjmatrix)
    dd = 0
    for i in range(len(vdist) - 1):
        for j in range(i + 1, len(vdist)):
            if adjmatrix[i, j] > 0:
                dd += vdist[i] + vdist[j]
    return dd


def gutman(adjmatrix):
    """

    :param adjmatrix:
    :return:
    """
    dmatrix = gp.distance_matrix(adjmatrix)
    degs = gp.degree(adjmatrix)
    gut = 0
    for i in range(len(degs) - 1):
        for j in range(i + 1, len(degs)):
            gut += degs[i] * degs[j] * dmatrix[i, j]
    return gut


def pdegree_kirchhoff(adjmatrix):
    """

    Parameters
    ----------
    adjmatrix

    Returns
    -------

    """
    rdmatrix = gp.resistance_distance_matrix(adjmatrix)
    degs = gp.degree(adjmatrix)
    pdkf = 0
    for i in range(len(degs) - 1):
        for j in range(i + 1, len(degs)):
            pdkf += degs[i] * degs[j] * rdmatrix[i, j]
    return pdkf
