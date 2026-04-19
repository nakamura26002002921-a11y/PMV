import numpy as np
import time
from scipy.spatial import Voronoi
import MDAnalysis as mda
import pmv_calc

def load_system(pdb, xtc):
    u = mda.Universe(pdb, xtc)
    u.guess_TopologyAttrs(to_guess=["elements"])
    protein = u.select_atoms("protein")
    waters = u.select_atoms("resname SOL and name O")
    return u, protein, waters

def process_frame(protein, waters, radius):
    points_atoms = protein.positions.copy()
    water_coords = waters.positions.copy()
    min_xyz = np.min(points_atoms - 2 * radius, axis=0)
    max_xyz = np.max(points_atoms + 2 * radius, axis=0)
    corners = np.array([
        [min_xyz[0], min_xyz[1], min_xyz[2]],
        [min_xyz[0], min_xyz[1], max_xyz[2]],
        [min_xyz[0], max_xyz[1], min_xyz[2]],
        [min_xyz[0], max_xyz[1], max_xyz[2]],
        [max_xyz[0], min_xyz[1], min_xyz[2]],
        [max_xyz[0], min_xyz[1], max_xyz[2]],
        [max_xyz[0], max_xyz[1], min_xyz[2]],
        [max_xyz[0], max_xyz[1], max_xyz[2]]
    ])
    points = np.concatenate([points_atoms, corners], axis=0)
    start = time.perf_counter()
    vor = Voronoi(points)
    total_volume = pmv_calc.compute_volume(
        [list(r) for r in vor.regions if r and -1 not in r],
        vor.vertices.astype(np.float64).tolist(),
        points_atoms.astype(np.float64).tolist(),
        float(radius)
    )
    end = time.perf_counter()
    diff = points_atoms[:, None, :] - water_coords[None, :, :]
    distances = np.linalg.norm(diff, axis=2)
    mask = distances < radius
    number = np.sum(np.any(mask, axis=0))
    return end - start, total_volume, number

def main():
    xtc = "example/centered.xtc" # Define input files for analysis
    pdb = "example/md.pdb" # Define input files for analysis
    u, protein, waters = load_system(pdb, xtc)
    protein_heavy = protein.select_atoms("not name H*")
    with open("result.txt", "w") as f:
        for ts, _ in enumerate(u.trajectory[:10000]):
            print(ts)
            t, vol, n = process_frame(protein_heavy, waters, radius=10)
            f.write(f"{ts}\t{t}\t{vol}\t{n}\n")

if __name__ == "__main__":
    main()
