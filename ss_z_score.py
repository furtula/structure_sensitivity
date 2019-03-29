"""This program calculates the structure sensitivity and abruptness
of topological indices from a given sqlite-database and a value of
the tanimoto similarity index."""

import sqlite3 as sql
import click
import pprint
import progressbar as pbar
import csv
import os

def db_cursor(name):
    """
    Makes cursor object for sqlite database if it exists or raise an error.

    Parameters
    ----------
    name : str
           The name of the sqlite database.

    Returns
    -------
    cur : sqlite3.cursor
          The cursor of the sqlite-database with given 'name'.

    Raises
    ------
    NameError
        If database 'name' does not exist.
    """
    if name in os.listdir():
        con = sql.connect(name)
        return con.cursor()
    else:
        raise NameError('The database "{}" does not exist!'.format(name))

def table_column_names(dbName, tabName):
    """
    Makes a list with names of columns in a table of a given sqlite-database.
        :param dbName: The sqlite-database's name.
        :param tabName: The table's name.
    """
    cur = db_cursor(dbName)
    cur.execute('PRAGMA table_info({})'.format(tabName))
    names = []
    for row in cur.fetchall():
        names.append(row[1])
    cur.close()
    return names

def make_dict_with_names_and_values(dbName, tabName, tabID, index):
    """
    Makes a dictionary with a names of topological indices as keys and
    their values for a molecule of a given index.
        :param dbName: The sqlite-database's name.
        :param tabName: The table's name.
        :param tabID: The table ID number.
        :param index: An index of a molecule in a molecules' table in sqlite
                      database.
    """
    cur = db_cursor(dbName)
    names = table_column_names(dbName, tabName)
    cur.execute('SELECT * FROM {} WHERE {}={}'.format(tabName,
                                                      tabID,
                                                      index))
    row = cur.fetchall()
    cur.close()
    return dict(zip(names, row[0]))

def select_similar(dbName, refID, z_score=.7):
    """
    Choosing ID's of a molecules similar to the molecule with refID.
        :param dbName: The sqlite-database's name.
        :param refID: The ID of a referent molecule.
        :param z_score=.7: The value of the z_score calculated on Tanimoto index similarity index.
            Defalut value is .7 because p for that value is greater than 75% .
    """
    cur = db_cursor(dbName)
    cur.execute("""
    SELECT * FROM tanimoto 
    WHERE (rootMol = {0} OR destMol = {0}) AND z_score >= {1}
    """.format(refID, z_score))
    rows = cur.fetchall()
    similar_ids = []
    for row in rows:
        if row[1] == int(refID):
            similar_ids.append(row[2])
        else:
            similar_ids.append(row[1])
    return similar_ids

def make_ss_dict(db_name):
    """
    Constructing a python dictionary of dictionaries with table names as
    keys and values are dictionaries containing names of topological indices
    as keys, and 0's as values.
        :param db_name: The sqlite-database's name.
    """
    not_allowed_table_names = ['molecules', 'tanimoto']

    ss = dict()

    con = sql.connect(db_name)
    cur = con.cursor()
    cur.execute("SELECT name FROM sqlite_master WHERE type='table'")
    for (i,) in cur.fetchall():
        if i not in not_allowed_table_names:
            tabColNames = table_column_names(db_name, i)
            ss.setdefault(i, {})
            for name in tabColNames:
                if name.endswith('ID'):
                    continue
                else:
                    ss[i].setdefault(name, 0)
    cur.close()
    con.close()
    return ss

def max_min(db_name, tab_name, ti_name):
    """
    Finding maximal and minimal value of topological index in a given
    sqlite-database.
        :param db_name: The sqlite-database's name.
        :param tab_name: The table's name.
        :param ti_name: The name of molecular descriptor.
    """
    cur = db_cursor(db_name)
    cur.execute('SELECT max({}) FROM {}'.format(ti_name, tab_name))
    maxTI = cur.fetchall()[0][0]
    cur.execute('SELECT min({}) FROM {}'.format(ti_name, tab_name))
    minTI = cur.fetchall()[0][0]
    cur.close()
    return maxTI, minTI


@click.command()
@click.argument('dbname')
@click.option('--zscore',
              '-z',
              default=.7,
              type=float,
              help='z-score value.')
def main(dbname, zscore):
    """
    Main function.
        :param dbname: The sqlite-database's name.
        :param z_score=.7: The value of the z_score calculated on Tanimoto index similarity index.
            Defalut value is .7 because p for that value is greater than 75% .
    """
    ss = make_ss_dict(dbname)
    abr = make_ss_dict(dbname)
    cur = db_cursor(dbname)
    cur.execute('SELECT COUNT(*) FROM molecules')
    n = cur.fetchall()[0][0] * len(ss)
    f = lambda x, y: (x - y) * (x - y)
    nn = 0
    with pbar.ProgressBar(max_value=n) as bar:
        for key in ss.keys():
            tabid = key + "ID"
            cur.execute('SELECT moleculeID FROM molecules')
            ukupan_broj = 0
            for (i,) in cur.fetchall():
                nn += 1
                ukupan_broj += 1
                bar.update(nn)
                refDict = make_dict_with_names_and_values(dbname, key, tabid, i)
                similarIds = select_similar(dbname, i, z_score=zscore)
                broj_slicnih = 0
                privremeniSS = make_ss_dict(dbname)
                privremeniAbr = make_ss_dict(dbname)
                if not similarIds:
                    continue
                for simID in similarIds:
                    broj_slicnih += 1
                    similarDicts = make_dict_with_names_and_values(dbname, 
                                                            key,
                                                            tabid,
                                                            simID)
                    for tikey in refDict.keys():
                        if tikey.endswith('ID'):
                            continue
                        privremeniSS[key][tikey] += \
                        f(refDict[tikey], similarDicts[tikey])
                        privremeniAbr[key][tikey] = \
                        max(privremeniAbr[key].get(tikey, 0),\
                        f(refDict[tikey], similarDicts[tikey]))

                if broj_slicnih == 0:
                    continue
                for k in ss[key].keys():
                    ss[key][k] += (privremeniSS[key][k] / broj_slicnih) ** .5
                    abr[key][k] += (privremeniAbr[key][k]) ** .5
            for k in ss[key].keys():
                maxti, minti = max_min(dbname, key, k)
                if maxti - minti == 0:
                    ss[key][k] = 0
                    abr[key][k] = 0
                else:
                    ss[key][k] = ss[key][k] / (ukupan_broj * (maxti - minti))
                    abr[key][k] = abr[key][k] / (ukupan_broj * (maxti - minti))
            
        

        for key in ss.keys():
            name = dbname[:-3] + "_" + key + '_Z' + str(zscore) + '.csv'
            with open(name, 'w') as f:
                f.write('topological_index,ss,abr\n')
                for k in ss[key].keys():
                    f.write('{},{:.3f},{:.3f}\n'.format(k,
                                                        ss[key][k],
                                                        abr[key][k]))

        pprint.pprint(ss)
        print('\n'*3)
        pprint.pprint(abr)
        cur.close()
    
if __name__ == "__main__":
    main()
