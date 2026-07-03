# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 13:15:51 2025

@author: biswas.161
"""

#!/usr/bin/env python3
import MDAnalysis as mda
from MDAnalysis.coordinates.PDB import PDBWriter

#Input gro file, trajectory file (.xtc or .trr), output prefix for file name
topology = "md_0.gro"     # or system.gro
trajectory = "md_noPBC.xtc"
output_prefix = "snapshot"

#List of time stamps (in ns) where you want snapshots
snapshot_times_ns = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

#Load trajectory
u = mda.Universe(topology, trajectory)
dt_ps = u.trajectory.dt
total_time_ps = u.trajectory.totaltime

print(f"Trajectory info: {len(u.trajectory)} frames, dt = {dt_ps} ps, total = {total_time_ps/1000:.2f} ns")

#Extract snapshot
for t_ns in snapshot_times_ns:
    target_time_ps = t_ns * 1000
    if target_time_ps > total_time_ps:
        print(f"Requested {t_ns:.2f} ns > trajectory length {total_time_ps/1000:.2f} ns")
        continue

    frame_index = int(round(target_time_ps / dt_ps))
    u.trajectory[frame_index]

    outname = f"{output_prefix}_{t_ns}ns.pdb"
    with PDBWriter(outname, multiframe=False) as pdb:
        pdb.write(u.atoms)

    print(f"Wrote {outname} at {t_ns:.2f} ns (frame {frame_index})")
