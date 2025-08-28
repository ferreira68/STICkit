# Paste into a Jupyter notebook cell
import nglview as nv
import MDAnalysis as mda

u = mda.Universe("topology.psf", "trajectory.dcd")  # replace with your files
nv.show_mdanalysis(u)
