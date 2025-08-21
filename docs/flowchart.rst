STICkit workflow
================

.. mermaid::

   flowchart TD
     A[Start] --> B[Load config (stickit.yaml)]
     B --> C[Load inputs (SMILES/SDF/CSV)]
     C --> D[Sanitize & standardize
            • remove salts/solvents
            • neutralize if needed
            • rdkit.Mol validation]
     D --> E{Enumerate\nStereoisomers}
     E --> F{Enumerate\nTautomers}
     F --> G{Enumerate\nIonization states}
     G --> H[Generate conformers (ETKDG)]
     H --> I[Optimize + rank (MMFF/UFF)]
     I --> J[Prune near-duplicates\n(RMSD/energy)]
     J --> K[Compose STIC records
            (S,T,I,C indices + metadata)]
     K --> L[Deduplicate microstates
            (canonical SMILES/InChIKey)]
     L --> Z[End]
