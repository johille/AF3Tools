from pymol import cmd
import csv
import math
import argparse
import os

#TODO
# which modes for executing this tool. Single structure - just get 
# write help for arguments
# "/home/jhille/ownCloud/Praktikum_Haubrock/Data/Predicted_Structures/FOS_JUN_heterodimer_targetdna_model.cif" 
""" 
pymol_cif_similarity.py - A tool to acquire potential binding sites in cif molecules and comparing multiple structures.
"""

def single_facilitate(cif_path, cutoff):
    obj = "structure"

    cmd.load(cif_path, obj)
    chains = cmd.get_chains(obj)
    print(chains)

def single_result(cif_path, cutoff):
    obj = "m"
    chainA = "A" # the chains should probably not be hardcoded, but taken from a getter()
    chainB = "B"

    out_csv = "/home/jhille/test_angstrom/interchain_residue_pairs_5A_with_mindist.csv" # this should only be done, if the user wants it, otherwise output on cli

    cmd.load(cif_path, obj)

    print(cmd.get_chains(obj))

    selA = f"{obj} and chain {chainA}"
    selB = f"{obj} and chain {chainB}"

    # Preload all atoms from chainB for quick index->atom lookup
    # We'll still only evaluate pairs that are within cutoff (from cmd.index query).
    pairs = {}  # (resiA,resnA,resiB,resnB) -> min_dist

    # Iterate residues in chain A (CA atoms used to enumerate residues)
    modelA = cmd.get_model(f"({selA}) and name CA")
    for a in modelA.atom:
        resiA, resnA = a.resi, a.resn

        # selection string for residue A (all atoms of that residue)
        sel_resA = f"{obj} and chain {chainA} and resi {resiA}"

        # Find nearby atoms from chainB within cutoff of residue A
        near_idx = cmd.index(f"({selB}) within {cutoff} of ({sel_resA})")
        if not near_idx:
            continue

        # Collect candidate residues in chainB that appear among those near atoms
        cand_resB = set()
        for (_o, idx) in near_idx:
            b = cmd.get_model(f"index {idx}").atom[0]
            cand_resB.add((b.resi, b.resn))

        # Get coordinates for all atoms in residue A once
        atomsA = cmd.get_model(sel_resA).atom
        coordsA = [(at.coord[0], at.coord[1], at.coord[2]) for at in atomsA]

        for resiB, resnB in cand_resB:
            sel_resB = f"{obj} and chain {chainB} and resi {resiB}"

            atomsB = cmd.get_model(sel_resB).atom
            coordsB = [(bt.coord[0], bt.coord[1], bt.coord[2]) for bt in atomsB]

            # Compute minimal atom-atom distance between the two residues
            min_d2 = None
            for (x1,y1,z1) in coordsA:
                for (x2,y2,z2) in coordsB:
                    dx = x1-x2; dy = y1-y2; dz = z1-z2
                    d2 = dx*dx + dy*dy + dz*dz
                    if (min_d2 is None) or (d2 < min_d2):
                        min_d2 = d2

            if min_d2 is None:
                continue

            d = math.sqrt(min_d2)

            key = (resiA, resnA, resiB, resnB)
            if key not in pairs or d < pairs[key]:
                pairs[key] = d

    # Write CSV
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["chainA","resiA","resnA","chainB","resiB","resnB","min_dist_A"])
        complete_distance = 0.0
        # sort lexicographically (resi can be string); for numeric-only resi this is fine
        for (resiA,resnA,resiB,resnB), d in sorted(pairs.items(), key=lambda kv: (kv[0][0], kv[0][2])):
            complete_distance += d
            w.writerow([chainA, resiA, resnA, chainB, resiB, resnB, f"{d:.3f}"])
        average_distance = complete_distance / len(pairs)
        print(average_distance)

    print(f"Wrote {len(pairs)} residue pairs with min distance to {out_csv}")
    cmd.quit()


def main():
    argparser = argparse.ArgumentParser(description="Pymol_Cif_Similarity - A tool to acquire potential binding sites in cif molecules and comparing multiple structures.")
    argparser.add_argument("-c", "--cutoff", type=float, default=5.0, required=False)
    mode = argparser.add_subparsers(dest="cmd", required=True)

    single = mode.add_parser("single")
    single.add_argument("cif")

    compare = mode.add_parser("compare")
    compare.add_argument("reference_cif")
    compare.add_argument("input_folder")

    args = argparser.parse_args()

    if args.cmd == "single":
        single_facilitate(args.cif, args.cutoff)
    elif args.cmd == "compare":
        r_cif = args.reference_cif
        input_folder = args.input_folder

if __name__ == "__main__":
    main()