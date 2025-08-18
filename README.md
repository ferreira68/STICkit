# stickit

STIC generation pipeline: Stereoisomer, Tautomer, Ionization state, Conformer.



a **STIC** is essentially the fully resolved microstate of a molecule in 3D, incorporating:

1. **Stereoisomer** – full stereochemical assignment (absolute configuration at chiral centers, double bond E/Z, atropisomerism if relevant).
2. **Tautomer** – fixed tautomeric form (proton location and formal bond arrangement locked).
3. **Ionization state** – fixed protonation/deprotonation and formal charges (specific pH-dependent microstate).
4. **Conformer** – a single geometry from the conformational ensemble.

That means a STIC is one node in a **hierarchical chemical state space**:

```
Molecule → set of Stereoisomers
    → set of Tautomers for each stereoisomer
        → set of Ionization states for each tautomer
            → set of Conformers for each ionization state
```

This makes STICs the most granular, physically-realizable “version” of a molecule you can feed into a physics-based engine like OpenMM or for which you can compute a single quantum chemistry result.

------

### Why the concept is useful

- **PhysChem property calculation** – Many properties (logP, pKa shifts, dipole moments, binding affinities) depend on the exact STIC.
- **Virtual screening** – Docking accuracy can change drastically across stereoisomers, tautomers, protonation states, and conformers.
- **Reproducibility** – Saying “molecule X” is ambiguous; saying “STIC #42 of molecule X” is unambiguous.
- **Microstate enumeration** – This is essentially the formalization of the microstate space for QSAR, docking, FEP, etc.



## Install

### Dependencies
STICkit requires AmberTools 23.6 or newer.
This can be install easily with conda or at the system level.
See the [AmberTools website](https://ambermd.org/Installation.php) for more details.

```bash
poetry install --with dev        # add --extras ray if you want Ray parallelism

