"""Program for alternating the original SQLite-database, adding
a column with z-score.
"""
import sqlite3 as sql
import click


def mean(X):
    return sum(X) / len(X)


def std(X):
    ex = mean(X)
    sum_of_suqares = [(x - ex) ** 2 for x in X]
    return (sum(sum_of_suqares) / (len(X) - 1)) ** .5


@click.command()
@click.argument('dbase')
def main(dbase):
    conn = sql.connect(dbase)
    cur = conn.cursor()
    cur.execute('ALTER TABLE tanimoto ADD COLUMN z_score float')
    cur.execute('SELECT * FROM tanimoto')
    tid = []
    tanimoto = []
    for i, _, _, t, _ in cur.fetchall():
        tid.append(i)
        tanimoto.append(t)

    ex = mean(tanimoto)
    sdev = std(tanimoto)

    zz = []

    for i, t in zip(tid, tanimoto):
        z = (t - ex) / sdev
        zz.append(round(z, 2))
        l = 'UPDATE tanimoto set z_score={:.2f} WHERE tanimotoID={}'.format(
            z, i)
        print(l)
        cur.execute(l)

    conn.commit()
    cur.close()
    conn.close()


if __name__ == '__main__':
    main()
