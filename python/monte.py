import numpy as np
import time
from scipy.spatial import Voronoi
import MDAnalysis as mda
import module

def load_system(pdb, xtc):
    u = mda.Universe(pdb, xtc)
    u.guess_TopologyAttrs(to_guess=["elements"])
    protein = u.select_atoms("protein")
    waters = u.select_atoms("resname SOL and name OW")
    return u, protein, waters

def monte_volume(protein, waters, r, max_iter=10, n_points=1000):
    atoms = protein.positions
    waters = waters.positions
    radius_dict = {'C':1.7,'H':1.2,'O':1.52,'N':1.55,'S':1.8}
    radiuses = np.array([radius_dict[e] for e in protein.elements])
    number = module.func1(atoms, radiuses, waters, r)
    times, results = [], []
    for _ in range(max_iter):
        start = time.perf_counter()
        result = module.func2(atoms, radiuses, r, n_points)
        end = time.perf_counter()
        times.append(end - start)
        results.append(result)
    return np.mean(times), np.mean(results), np.var(results), number

def main():
    xtc = "example/centered.xtc" # Define input files for analysis
    pdb = "example/md.pdb" # Define input files for analysis
    u, protein, waters = load_system(pdb, xtc)
    protein_heavy = protein.select_atoms("not name H*")
    with open("result.txt", "w") as f:
        for ts, _ in enumerate(u.trajectory[:10000]):
            print(ts)
            t, vol, n = monte_volume(protein_heavy, waters, radius=1.4)
            f.write(f"{ts}\t{t}\t{vol}\t{n}\n")

if __name__ == "__main__":
    main()
