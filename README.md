# Preprocessing tool for protein graph representation

This tool turns pattern of pdb files to featurized graph of atoms. Atoms are graph nodes. They're connected to each other based on distance within cutoff (default is 3.5A). Each atom also is featurized with few numerical data points:
* Wan der vaal radius
* Atomic number
* Valence electrons
* Electronegativity

For each pdb file, corresponding `<filename>_graph.json` is created with serialized graph representation.

## Usage

```
graphein --pdb-glob dir_with_pdbs/**/*.pdb --cutoff 3.5
```