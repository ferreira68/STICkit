
.. _stic_flowcharts:

STICkit Developer Flowcharts
============================

This page provides an architectural overview of the STICkit pipeline and focused subflows for developers.

Overview
--------

.. mermaid::
   :zoom:
   :caption: End to end STIC pipeline
   :align: center

   flowchart TD
     A[Start Inputs SMILES names] --> B[Load and Validate RDKit sanitize standardize]
     B --> C{Enumeration}
     C --> C1[Stereochemistry assign and enumerate]
     C --> C2[Tautomerization ruleset based]
     C --> C3[Ionization states pH windows]
     C1 --> D[Conformer generation ETKDG n]
     C2 --> D
     C3 --> D
     D --> E[Geometry optimization MMFF or UFF]
     E --> E0[Create simulation context AmberTools invoked here]
     E0 --> F[Property calculation and filtering energy RMSD rules]
     F --> G[Deduplicate and index STIC ID mol st tt ion cf]
     G --> H[Export SDF with 3D coords and props]
     H --> I[Downstream OpenMM Docking QSAR]

Entrypoint and Orchestration
----------------------------

.. mermaid::
   :zoom:
   :caption: Python API orchestration no CLI yet
   :align: center

   flowchart LR
     P1[Python API pipeline py stic_generation] --> P2[Parse config or args]
     P2 --> P3{Parallel backend}
     P3 --> P4[Local multiprocessing]
     P4 --> P5[Chunk inputs]
     P5 --> P6[Map tasks to workers]
     P6 --> P7[Collect results and retries]
     P7 --> P8[Write SDF output]
     P7 --> P9[Emit logs and metrics]

Input and Validation
--------------------

.. mermaid::
   :zoom:
   :caption: Input and preparation SMILES with names
   :align: center

   flowchart LR
     A1[Read SMILES and names] --> A2[Assign IDs source id]
     A2 --> A3[Salt or solvent strip optional]
     A3 --> A4[Neutralize and canonicalize]
     A4 --> A5[RDKit sanitize]
     A5 --> A6[Standardize aromaticity or kekulize]
     A6 --> A7[Optional initial 3D embed]

Conformers and Optimization
---------------------------

.. mermaid::
   :zoom:
   :caption: Conformer workflow with AmberTools at context creation
   :align: center

   flowchart LR
     E1[Enumerated microstate] --> F1[Conformer generation ETKDG n]
     F1 --> F2[Minimize MMFF or UFF]
     F2 --> F3[Create simulation context]
     F3 --> F4[AmberTools invoked antechamber GAFF]
     F4 --> G1[Filter by energy window and RMSD]
     G1 --> G2[Select N best per microstate]

Deduplication and Indexing
--------------------------

.. mermaid::
   :zoom:
   :caption: Uniqueness and STIC IDs
   :align: center

   flowchart TD
     H1[Conformers] --> H2[Canonicalization InChIKey charge aware]
     H2 --> H3{Duplicate}
     H3 -->|Yes| H4a[Drop or merge keep best]
     H3 -->|No| H5a[Unique]
     H4a --> H5[Assign STIC ID mol st tt ion cf]
     H5a --> H5
     H5 --> H6[Write records with 3D coords to SDF]


Enumeration Strategy
--------------------

.. mermaid::
   :zoom:
   :caption: Stereo then Tautomer then Ionization
   :align: center

   flowchart TD
     B1[Input molecule] --> B2{Has stereochemistry}
     B2 -->|No| B3[Pass through]
     B2 -->|Yes| B4[Enumerate chiral centers]
     B4 --> B5[Enumerate E Z and atropisomers]
     B3 --> C1{Tautomer needed}
     B5 --> C1
     C1 -->|Yes| C2[Generate tautomers]
     C1 -->|No| C3[Keep as is]
     C2 --> D1{Ionization}
     C3 --> D1
     D1 -->|Yes| D2[Enumerate protonation and charge microstates]
     D1 -->|No| E1[Label stereo id taut id ion id]
     D2 --> E1[Label stereo id taut id ion id]


Appendix Data Model
-------------------

.. mermaid::
   :zoom:
   :caption: STIC record schema
   :align: center

   classDiagram
     class STIC {
       string stic_id
       string parent_mol_id
       int stereo_id
       int taut_id
       int ion_id
       int conf_id
       string smiles
       string inchikey
       float energy
       float pH
       dict props
       bytes molblock3D
     }


Call Flow: stic_generation
--------------------------

.. mermaid::
   :zoom:
   :caption: Sequence of the main pipeline call
   :align: center

   sequenceDiagram
     participant Caller as Caller
     participant Pipeline as stic_generation
     participant Preprocess as Preprocess
     participant Enumerate as Enumerate
     participant Conformers as Conformers
     participant AmberTools as AmberTools Context
     participant FilterDedup as Filter Dedup
     participant Writer as SDF Writer

     Caller->>Pipeline: call with SMILES names config
     Pipeline->>Preprocess: load validate standardize
     Preprocess-->>Pipeline: molecules
     Pipeline->>Enumerate: run enumeration
     alt Stereochemistry present
       Enumerate-->>Pipeline: generated stereoisomers
     else No stereochemistry
       Enumerate-->>Pipeline: pass through
     end
     alt Tautomer rules apply
       Enumerate-->>Pipeline: generated tautomers
     else No tautomers
       Enumerate-->>Pipeline: keep as is
     end
     alt Ionization enumeration
       Enumerate-->>Pipeline: microstates
     else No ionization changes
       Enumerate-->>Pipeline: unchanged states
     end
     Pipeline->>Conformers: embed minimize select
     Conformers-->>Pipeline: conformer set
     Pipeline->>AmberTools: create simulation context
     AmberTools-->>Pipeline: charges parameters
     Pipeline->>FilterDedup: score filter deduplicate
     FilterDedup-->>Pipeline: final set
     Pipeline->>Writer: write SDF with 3D coords props
     Writer-->>Pipeline: output path
     Pipeline-->>Caller: summary counts ids
