# Structure sensitivity of topological molecular descriptors

The structure sensitivity of topological molecular descriptors is
an important feature that gives insight into their quality and
prediction potential. These Python scripts are used for assessing
the structure sensitivity of a bunch of topological descriptors in
the paper M. Rakić, B. Furtula, A novel method for measuring the
structure sensitivity of molecular descriptors, J. Chemom., submitted.

In order to obtain results, the set of molecules needs to be provided
as sdf-file (hydrocarbons in this particular case). As output the
SQLite base will be produced containing table of molecules in SMILES
format and the table of Tanimoto similarity measure calculated between
all pairs of molecules in database. Tanimoto similarity measure is
calculating based on the so-called Morgan circular fingerprints.
Furthermore, this DB contains table with a number of degree-based
topological indices, table with a number of distance-based topological
indices, and the table of a number of eigenvalue-based topological
indices. For obtaining this SQLite database, the following command needs
to be executed:

```Console
foo@bar:~$ python make_db_with_tis.py sdfFile outfile.db
```

The Tanimoto values from DB are standardized by the following
command:

```Console
foo@bar:~$ python z_score.py outfile.db
```
Finally, the structure sensitivity and abruptness of topological
indices are obtained by the following command:

```Console
foo@bar:~$ python ss_z_score.py outfile.db --zscore .7
```

The default value of `--zscore` is 0.7 and you may use this
option only if you want to change its default value.
