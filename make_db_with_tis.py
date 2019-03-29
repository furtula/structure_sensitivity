from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import click
import sqlite3 as sql
from moldes import degreeBased as degb
from moldes import distanceBased as distb
from moldes import distancedegreeBased as distdegb
from moldes import eigenBased as eigb


def is_alkane(molecule):
    """
    Checking if the molecule is alkane or not.
        :param molecule: The rdkit class of a molecule.
    """
    for bond in molecule.GetBonds():
        if bond.GetBondType().name != 'SINGLE':
            return False
    return True

def choose_alkanes(container):
    """
    Filtering a sdf-container for eliminating non-alkane-molecules.
        :param container: Container of molecules in sdf format.
    """
    molecules = []
    for molecule in container:
        if is_alkane(molecule):
            molecules.append(molecule)
    return molecules

def tanimoto(ref_MorganFP, test_MorganFP):
    """
    Calculating Tanimoto similarity index betwen the Morgan-type fingerprints
    between ref_ and test_ molecules. It is rounded at upto 2 decimal places.
        :param ref_MorganFP: Morgan-type fingerprint of the ref_ molecule.
        :param test_MorganFP: Morgan-type fingerprint of the test_ molecule.
    """
    return round(DataStructs.TanimotoSimilarity(ref_MorganFP, test_MorganFP), 2)

@click.command()
@click.argument('infile')
@click.argument('dbase')
def main(infile, dbase):
    molecules = Chem.SDMolSupplier(infile)
    alkanes = choose_alkanes(molecules)
    conn = sql.connect(dbase)
    cur = conn.cursor()
    cur.execute(
        """CREATE TABLE IF NOT EXISTS molecules(moleculeID INTEGER PRIMARY KEY 
        UNIQUE, molecule TEXT UNIQUE);
        """
    )


    cur.execute("""
    CREATE TABLE IF NOT EXISTS degree (
        degreeID INTEGER PRIMARY KEY UNIQUE,
        firstZagreb INTEGER,
        secondZagreb INTEGER,
        randic REAL,
        abc REAL,
        ga REAL,
        sum_connectivity REAL,
        azi REAL,
        harmonic REAL,
        sigma INTEGER,
        indeg REAL, 
        sdd REAL,
        inverseDegree REAL,
        forgotten INTEGER,
        reducedReciprocalRandic REAL,
        reciprocalRandic REAL,
        reducedSecondZagreb INTEGER
    )
    """)


    cur.execute("""
    CREATE TABLE IF NOT EXISTS distance (
        distanceID INTEGER PRIMARY KEY UNIQUE,
        wiener INTEGER,
        modiefied_wiener2 INTEGER,
        modified_wiener3 INTEGER,
        modified_wienerMINUS2 REAL,
        modified_wienerMINUS3 REAL,
        harary REAL,
        hyper_wiener INTEGER,
        tratch_stankievich_zefirov INTEGER,
        szeged INTEGER,
        padmakar_ivan INTEGER,
        ga2 REAL,
        graovac_ghorbani REAL,
        reciprocal_complementary_w REAL,
        balaban REAL,
        sum_balaban REAL,
        detour INTEGER,
        kirchhoff REAL
    )
    """)

    wiener = distb.DistanceTI(distb.w, 'wiener', x=1, choose='r')
    hyper_wiener = distb.DistanceTI(
        distb.ww, 'hyper_wiener', x=1, choose='r')
    tsz = distb.DistanceTI(
        distb.tsz, 'tratch_stankievich_zefirov', x=1, choose='r')
    mwiener2 = distb.DistanceTI(
        distb.w, 'modified_wiener(2)', x=2, choose='r')
    mwiener3 = distb.DistanceTI(
        distb.w, 'modified_wiener(3)', x=3, choose='r')
    harary = distb.DistanceTI(distb.harary, 'harary', x=1, choose='r')
    mwienerM2 = distb.DistanceTI(
        distb.w, 'modified_wiener(-2)', x=-2, choose='r')
    mwienerM3 = distb.DistanceTI(
        distb.w, 'modified_wiener(-3)', x=-3, choose='r')

    cur.execute("""
    CREATE TABLE IF NOT EXISTS distance_degree (
        distance_degreeID INTEGER PRIMARY KEY UNIQUE,
        eccentric_connectivity INTEGER,
        degree_distance INTEGER,
        gutman INTEGER,
        pdegree_kirchhoff REAL,
        sdegree_kirchhoff REAL
    )
    """)

    cur.execute("""
    CREATE TABLE IF NOT EXISTS eigenvalue (
        eigenvalueID INTEGER PRIMARY KEY UNIQUE,
        energy REAL,
        laplacian_energy REAL,
        distance_energy REAL,
        ex_adjacency_energy REAL,
        randic_energy REAL,
        harary_energy REAL,
        resolvent_energy REAL,
        incidence_energy REAL,
        laplacian_energy_like REAL,
        estrada REAL,
        laplacian_estrada REAL,
        distance_estrada REAL,
        randic_estrada REAL,
        harary_estrada REAL,
        ex_adjacency_estrada REAL,
        homo_lumo REAL
    )
    """)


    for i, alkane in enumerate(alkanes):
        cur.execute(
            """INSERT OR IGNORE INTO molecules (moleculeID, molecule)
            VALUES (?,?)""", (i+1, Chem.MolToSmiles(alkane)))

        a = Chem.GetAdjacencyMatrix(alkane)

        cur.execute("""
            INSERT OR IGNORE INTO degree (
                degreeID, firstZagreb, secondZagreb, randic, abc, ga, 
                sum_connectivity, azi, harmonic, sigma, indeg, sdd,
                inverseDegree, forgotten, reducedReciprocalRandic, 
                reciprocalRandic, reducedSecondZagreb)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """,(i+1, degb.first_zagreb(a), degb.second_zagreb(a, x=1),
        degb.second_zagreb(a, x=-.5), degb.atom_bond_connectivity(a),
        degb.geometric_arithmetic(a), degb.sum_connectivity(a, x=-.5),
        degb.atom_bond_connectivity(a, x=-3), 2*degb.harmonic(a, x=1),
        degb.sigma(a), degb.inverse_sum_indeg(a),
        degb.sdd(a), degb.first_zagreb(a, x=-2), degb.first_zagreb(a, x=2),
        degb.reduced_second_zagreb(a, x=.5),degb.second_zagreb(a, x=.5), 
        degb.reduced_second_zagreb(a)))

        cur.execute("""
            INSERT OR IGNORE INTO distance (
                distanceID, wiener, modiefied_wiener2, modified_wiener3,
                modified_wienerMINUS2, modified_wienerMINUS3, harary,
                hyper_wiener, tratch_stankievich_zefirov, szeged, padmakar_ivan,
                ga2, graovac_ghorbani, reciprocal_complementary_w, balaban,
                sum_balaban, detour, kirchhoff)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """, (i+1, wiener.ti(a),  mwiener2.ti(a),  mwiener3.ti(a),
              mwienerM2.ti(a),  mwienerM3.ti(a), harary.ti(a),
              hyper_wiener.ti(a), tsz.ti(a), distb.szeged(a), distb.pi(a),
              distb.second_geometric_arithmetic(a), distb.graovac_ghorbani(a),
              distb.reciprocal_complementary_w(a), distb.balaban(a),
              distb.balaban(a, sum_balaban=True), distb.detour(a),
              eigb.kirchhoff(a))
        )

        cur.execute("""
            INSERT OR IGNORE INTO distance_degree (
                distance_degreeID, eccentric_connectivity, degree_distance,
                gutman, pdegree_kirchhoff, sdegree_kirchhoff)
            VALUES (?,?,?,?,?,?)
        """, (i+1, distdegb.eccentric_connectivity(a),
             distdegb.degree_distance(a), distdegb.gutman(a),
             distdegb.pdegree_kirchhoff(a), distb.sdegree_kirchhoff(a))
        )

        numAtoms = alkane.GetNumAtoms()
        numBonds = alkane.GetNumBonds()

        cur.execute("""
            INSERT OR IGNORE INTO eigenvalue (
                 eigenvalueID, energy, laplacian_energy, distance_energy,
                 ex_adjacency_energy, randic_energy, harary_energy,
                 resolvent_energy, incidence_energy, laplacian_energy_like,
                 estrada, laplacian_estrada, distance_estrada, randic_estrada,
                 harary_estrada, ex_adjacency_estrada, homo_lumo)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """, (i + 1, eigb.energy(a, mat='adjacency', c=0), 
             eigb.energy(a, mat='laplacian', c=2 * numBonds / numAtoms),
             eigb.energy(a, mat='distance', c=0), eigb.energy(a, mat='xu', c=0),
             eigb.energy(a, mat='randic', c=0), eigb.energy(a, mat='harary', c=0),
             eigb.resolvent_energy(a), eigb.incidence_energy(a), eigb.lel(a),
             eigb.estrada(a, mat='adjacency', c=0),
             eigb.estrada(a, mat='laplacian', c=2 * numBonds / numBonds),
             eigb.estrada(a, mat='distance', c=0), 
             eigb.estrada(a, mat='randic', c=0),
             eigb.estrada(a, mat='harary', c=0),
            eigb.estrada(a, mat='xu', c=0),
            eigb.homolumo(a))
        )





    cur.execute("""
    CREATE TABLE IF NOT EXISTS tanimoto (
        tanimotoID INTEGER PRIMARY KEY UNIQUE,
        rootMol INTEGER,
        destMol INTEGER,
        value REAL,
        FOREIGN KEY(rootMol) REFERENCES molecules(moleculeID),
        FOREIGN KEY(destMol) REFERENCES molecules(moleculeID)
    )
    """)

    rb = 0

    for i, ref_alkane in enumerate(alkanes):
        fp_ref = AllChem.GetMorganFingerprint(ref_alkane, 2)
        for j, test_alkane in enumerate(alkanes):
            if j <= i:
                continue
            else:
                rb += 1
                fp_test = AllChem.GetMorganFingerprint(test_alkane, 2)
                cur.execute(
                    """INSERT OR IGNORE INTO tanimoto (tanimotoID, rootMol,
                    destMol, value) VALUES (?, ?, ?, ?)
                    """, (rb, i+1, j+1, tanimoto(fp_ref, fp_test)))

    conn.commit()
    cur.close()
    conn.close()



if __name__ == "__main__":
    main()
