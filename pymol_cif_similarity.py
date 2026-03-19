from pymol import cmd
import csv
import math
import argparse
import os
from itertools import combinations

#TODO
# implement compare mode
# split calculating of things and writing of csv into 2 functions for more modularity
# "/home/jhille/ownCloud/Praktikum_Haubrock/Data/Predicted_Structures/FOS_JUN_heterodimer_targetdna_model.cif" 
""" 
pymol_cif_similarity.py - A tool to acquire potential binding sites in cif molecules and comparing multiple structures.
"""

def write_to_csv_whole(csv_path, pairs, result_pairs):
    result_map = {tuple(pair): dist_dict for pair, dist_dict in result_pairs}

    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["chain1", "resi1", "resn1", "chain2", "resi2", "resn2", "min_dist_A"])

        for chainA, chainB in pairs:
            dist_dict = result_map.get((chainA, chainB), {})
            # sort by residue indices (resi are strings; if numeric, this is OK; otherwise lexicographic)
            for (resiA, resnA, resiB, resnB), d in sorted(dist_dict.items(), key=lambda kv: (kv[0][0], kv[0][2])):
                w.writerow([chainA, resiA, resnA, chainB, resiB, resnB, f"{d:.3f}"])

            w.writerow([])

    print(f"Wrote combined CSV to {csv_path}")

def get_average_dist(pairs, result_pairs):
    result_map = {tuple(pair): dist_dict for pair, dist_dict in result_pairs}

    global_sum = 0.0
    global_n = 0

    for chainA, chainB in pairs:
        dist_dict = result_map.get((chainA, chainB), {})

        pair_n = 0
        pair_sum = 0

        pair_n = len(dist_dict)
        pair_sum = sum(dist_dict.values())


        if pair_n > 0:
            print(f"average binding distance {chainA}-{chainB}: {pair_sum/pair_n:.3f} Å ({pair_n} contacts)")
            global_sum += pair_sum
            global_n += pair_n
        else:
            print(f"no contacts within cutoff for {chainA}-{chainB}")

    if global_n > 0:
        print(f"global average binding distance: {global_sum/global_n:.3f} Å ({global_n} contacts)")
        return global_sum/global_n, global_n
    else:
        print("no contacts within cutoff for any chain pair")
        return 0.0, 0.0


def get_binding_pairs(cutoff, obj, chainA, chainB):
    selA = f"{obj} and chain {chainA}"
    selB = f"{obj} and chain {chainB}"

    # Preload all atoms from chainB for quick index->atom lookup
    # We'll still only evaluate pairs that are within cutoff (from cmd.index query).
    pairs = {}  # (resiA,resnA,resiB,resnB) -> min_dist

    # Iterate residues in chain A (CA atoms used to enumerate residues)
    modelA = cmd.get_model(f"({selA}) and name CA")
    for a in modelA.atom:
        resiA, resnA = a.resi, a.resn #resiA = Index, resnA = amino acid or nucleotide

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
    return pairs

def single_facilitate(cif_path, cutoff, output_path):
    obj = "structure"
    cmd.load(cif_path, obj)
    chains = cmd.get_chains(obj)
    pairs = list(combinations(chains, 2))

    result_pairs = []
    for pair in pairs:
        result_pairs.append([pair, get_binding_pairs(cutoff, obj, pair[0], pair[1])])

    _ = get_average_dist(pairs, result_pairs)
    write_to_csv_whole(output_path, pairs, result_pairs)

    cmd.quit()

def compare_facilitate(cif_path, cutoff, input_folder ,output_path):
    ref_obj = "reference_structure"
    cmd.load(cif_path, ref_obj)
    ref_chains = cmd.get_chains(ref_obj)
    ref_pairs = list(combinations(ref_chains, 2))
    ref_result_pairs = []

    for pair in ref_pairs:
        ref_result_pairs.append([pair, get_binding_pairs(cutoff, ref_obj, pair[0], pair[1])])

    
    ref_avg_dist = get_average_dist(ref_pairs, ref_result_pairs)[0]
    
    input = os.listdir(input_folder)
    avg_distances_structures = []
    for structure in input:
        if ".cif" in structure:
            cmd.load(input_folder + "/" + structure, structure)
            chains = cmd.get_chains(structure)
            pairs = list(combinations(chains, 2))
            result_pairs = []

            for pair in pairs:
                result_pairs.append([pair, get_binding_pairs(cutoff, structure, pair[0], pair[1])])
            avg_dist = get_average_dist(pairs, result_pairs)[0]
            diff = abs(ref_avg_dist - avg_dist)
            avg_distances_structures.append([structure, avg_dist, diff])
    avg_distances_structures.sort(key=lambda x: x[2])
    print(avg_distances_structures)
    pass

def main():
    argparser = argparse.ArgumentParser(description="Pymol_Cif_Similarity - A tool to acquire potential binding sites in cif molecules and comparing multiple structures.")
    argparser.add_argument("-c", "--cutoff", type=float, default=5.0, required=False)
    argparser.add_argument("-o", "--output_path", type=str, required=True)
    mode = argparser.add_subparsers(dest="cmd", required=True)

    single = mode.add_parser("single")
    single.add_argument("cif")

    compare = mode.add_parser("compare")
    compare.add_argument("reference_cif")
    compare.add_argument("--input_folder", type=str, required=True)

    args = argparser.parse_args()

    if args.cmd == "single":
        single_facilitate(args.cif, args.cutoff, args.output_path)
    elif args.cmd == "compare":
        compare_facilitate(args.reference_cif, args.cutoff, args.input_folder ,args.output_path)
        #r_cif = args.reference_cif
        #input_folder = args.input_folder

if __name__ == "__main__":
    main()