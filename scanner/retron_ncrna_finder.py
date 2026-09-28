#!/usr/bin/env python3
"""
retron_ncrna_finder.py
----------------------
Identify candidate retron ncRNA (msr-msd) genes in nucleotide sequences.

Conserved features used for detection:
  1. Inverted repeat (a1/a2 region) ≥ 8 bp with a configurable max loop
  2. Priming guanines at the inner junction of each arm:
       5'-a1-G -- [internal region] -- G-a2-3'
     The G is the last base of a1 and the first base of a2 (where the
     stem breaks); these two Gs are the RT priming sites.
  3. Optional length constraints on the candidate region

Usage:
  python retron_ncrna_finder.py -i genome.fasta [options]

Output:
  FASTA file of candidate ncRNA sequences with annotation in the header.
"""

import argparse
import sys
import re
from pathlib import Path
from itertools import product as iterproduct


# ---------------------------------------------------------------------------
# FASTA I/O
# ---------------------------------------------------------------------------

def parse_fasta(path):
    """Yield (header, sequence) tuples from a FASTA file."""
    header, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks).upper()
                header, chunks = line[1:], []
            elif line:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks).upper()


def reverse_complement(seq):
    comp = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
    return seq.translate(comp)[::-1]


# ---------------------------------------------------------------------------
# Inverted-repeat search
# ---------------------------------------------------------------------------

WOBBLE = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}

def iupac_match(b1, b2):
    """Return True if b2 is the Watson-Crick complement of b1 (IUPAC-aware)."""
    complement = {"A": "T", "T": "A", "G": "C", "C": "G",
                  "R": "Y", "Y": "R", "S": "S", "W": "W",
                  "K": "M", "M": "K", "B": "V", "V": "B",
                  "D": "H", "H": "D", "N": "N"}
    c1 = complement.get(b1, "N")
    # check if b2 could match c1 (both may be ambiguous)
    def expands(b):
        table = {
            "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"},
            "R": {"A","G"}, "Y": {"C","T"}, "S": {"G","C"},
            "W": {"A","T"}, "K": {"G","T"}, "M": {"A","C"},
            "B": {"C","G","T"}, "D": {"A","G","T"},
            "H": {"A","C","T"}, "V": {"A","C","G"}, "N": {"A","C","G","T"},
        }
        return table.get(b, {"N"})
    return bool(expands(c1) & expands(b2))


def find_inverted_repeats(seq, min_arm=8, max_loop=400, max_mismatch=1):
    """
    Find inverted repeats (stem-loops) in seq.

    Yields dicts:
      arm1_start, arm1_end, arm2_start, arm2_end, arm_len, loop_len, mismatches
    where arm1 is the 5' arm and arm2 is the 3' arm.
    """
    n = len(seq)
    results = []

    for i in range(n - 2 * min_arm):
        for arm_len in range(min_arm, (n - i) // 2 + 1):
            # arm2 begins at least min_arm + 1 nt after arm1 ends (loop ≥ 1)
            arm1 = seq[i : i + arm_len]
            arm2_rc_needed = reverse_complement(arm1)

            for j in range(i + arm_len + 1, n - arm_len + 1):
                loop_len = j - (i + arm_len)
                if loop_len > max_loop:
                    break
                arm2 = seq[j : j + arm_len]

                # Count mismatches
                mm = sum(1 for a, b in zip(arm2_rc_needed, arm2) if a != b)
                if mm <= max_mismatch:
                    results.append({
                        "arm1_start": i,
                        "arm1_end":   i + arm_len,
                        "arm2_start": j,
                        "arm2_end":   j + arm_len,
                        "arm_len":    arm_len,
                        "loop_len":   loop_len,
                        "mismatches": mm,
                    })

    # Remove strictly dominated hits (shorter arm, same or worse mm at same position)
    results.sort(key=lambda r: (-r["arm_len"], r["mismatches"]))
    return results


def find_inverted_repeats_fast(seq, min_arm=8, max_loop=400, max_mismatch=1):
    """
    Find inverted repeats (IRs) in seq.

    Structure:
        5'-[arm1]---[loop]---[arm2]-3'
           i..i+L           j..j+L
        where seq[j:j+L] == RC(seq[i:i+L])  (arm2 is the reverse complement of arm1)

    Strategy: iterate over all possible (arm1_inner_end, arm2_inner_start) pairs,
    i.e. the two positions immediately flanking the loop.  Then extend the arm
    outward one base at a time, comparing the new arm1 outer base with the
    complement of the new arm2 outer base.

    At extension step k (0-indexed, growing outward from the loop):
        arm1 base: seq[arm1_inner_end - 1 - k]   (walking left/5' from loop)
        arm2 base: seq[arm2_inner_start + k]      (walking right/3' from loop)
        These must be Watson-Crick complements:
            seq[arm1_inner_end - 1 - k]  ==  complement(seq[arm2_inner_start + k])

    This is correct because in a stem-loop:
        the innermost pair is (arm1[-1], arm2[0])
        the next pair out is  (arm1[-2], arm2[1])
        etc.
    """
    n = len(seq)
    comp_table = str.maketrans("ACGTN", "TGCAN")
    best = {}   # (e, s) -> best hit dict, where e=arm1_inner_end, s=arm2_inner_start

    # e = arm1_inner_end (exclusive): arm1 ends at e, loop starts at e
    # s = arm2_inner_start: loop ends at s-1, arm2 starts at s
    # loop_len = s - e  (must be >= 1)
    # arm_len  = k + 1 for extension step k
    # arm1 = seq[e - arm_len : e]
    # arm2 = seq[s : s + arm_len]

    for e in range(min_arm, n - min_arm - 1):           # arm1 ends here (exclusive)
        for s in range(e + 1, min(e + max_loop + 1, n - min_arm + 1)):  # arm2 starts here
            loop_len = s - e
            if loop_len < 1:
                continue

            mm = 0
            arm_len = 0
            max_arm = min(e, n - s)   # arm1 can grow at most e bases left; arm2 fits in seq

            while arm_len < max_arm:
                # New outermost pair
                a1_base = seq[e - 1 - arm_len]          # next base left in arm1
                a2_base = seq[s + arm_len]               # next base right in arm2
                if a1_base == a2_base.translate(comp_table):
                    arm_len += 1
                elif mm < max_mismatch:
                    mm += 1
                    arm_len += 1
                else:
                    break

                if arm_len >= min_arm:
                    key = (e, s)
                    if key not in best or arm_len > best[key]["arm_len"] or \
                       (arm_len == best[key]["arm_len"] and mm < best[key]["mismatches"]):
                        best[key] = {
                            "arm1_start": e - arm_len,
                            "arm1_end":   e,
                            "arm2_start": s,
                            "arm2_end":   s + arm_len,
                            "arm_len":    arm_len,
                            "loop_len":   loop_len,
                            "mismatches": mm,
                        }

    return list(best.values())


# ---------------------------------------------------------------------------
# Priming-G detection
# ---------------------------------------------------------------------------

def check_junction_gs(seq, ir):
    """
    Check for the conserved priming guanines at the inner junction of the
    inverted-repeat arms.

    Retron ncRNA structure:

        5'-[a1 arm]-G -- [internal msr/msd region] -- G-[a2 arm]-3'
                    ↑                                  ↑
               arm1_end                          arm2_start - 1

    The junction Gs sit immediately OUTSIDE the IR arms:
      - 5′ junction G: seq[ir["arm1_end"]]       (first base after a1 ends)
      - 3′ junction G: seq[ir["arm2_start"] - 1] (last base before a2 starts)

    These Gs are NOT part of the inverted-repeat arms themselves; they mark
    the boundary between the stem and the internal msd region and are the
    sites used for reverse-transcription priming.

    Returns:
        (has_5p_g, has_3p_g, g5_base, g3_base)
          has_5p_g  – True if the 5′ junction base is G
          has_3p_g  – True if the 3′ junction base is G
          g5_base   – actual nucleotide at the 5′ junction position
          g3_base   – actual nucleotide at the 3′ junction position
    """
    pos5 = ir["arm1_end"]           # first base after a1 (0-based in seq)
    pos3 = ir["arm2_start"] - 1     # last base before a2 (0-based in seq)

    g5_base = seq[pos5]
    g3_base = seq[pos3]

    return g5_base == "G", g3_base == "G", g5_base, g3_base


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def score_candidate(ir, has_5p_g, has_3p_g):
    """
    Simple heuristic score (higher = more likely ncRNA).
    Components:
      - arm length (longer → higher)
      - loop length penalty (very long loops are less favourable)
      - mismatch penalty
      - bonus for each junction G present (5′ and/or 3′)
    """
    score = ir["arm_len"] * 10
    score -= ir["loop_len"] * 0.05
    score -= ir["mismatches"] * 5
    if has_5p_g:
        score += 10   # 5′ junction G present
    if has_3p_g:
        score += 10   # 3′ junction G present
    return round(score, 2)


# ---------------------------------------------------------------------------
# Main search
# ---------------------------------------------------------------------------

def search_sequence(seq_id, seq, args):
    """
    Search a single nucleotide sequence for retron ncRNA candidates.
    Returns a list of candidate dicts.
    """
    candidates = []
    seq_len = len(seq)

    # Apply optional search-region constraints
    start_offset = max(0, args.region_start - 1) if args.region_start else 0
    end_offset   = min(seq_len, args.region_end)   if args.region_end   else seq_len
    subseq = seq[start_offset:end_offset]
    sub_len = len(subseq)

    # Search both strands
    for strand, working_seq in (("+", subseq), ("-", reverse_complement(subseq))):
        irs = find_inverted_repeats_fast(
            working_seq,
            min_arm=args.min_arm,
            max_loop=args.max_loop,
            max_mismatch=args.max_mismatch,
        )

        for ir in irs:
            # The candidate ncRNA spans from the start of arm1 to the end of arm2
            cand_start = ir["arm1_start"]
            cand_end   = ir["arm2_end"]
            cand_len   = cand_end - cand_start

            # Length filter
            if cand_len < args.min_length or cand_len > args.max_length:
                continue

            cand_seq = working_seq[cand_start:cand_end]

            # Priming-G check: the conserved Gs are at the inner junctions
            # of the two IR arms within working_seq (not within cand_seq),
            # because ir coordinates are relative to working_seq.
            has_5p_g, has_3p_g, g5_base, g3_base = check_junction_gs(working_seq, ir)

            # --require-g: both junction positions must be G
            if args.require_g and not (has_5p_g and has_3p_g):
                continue

            # Score
            sc = score_candidate(ir, has_5p_g, has_3p_g)

            # Convert coordinates back to original sequence space
            if strand == "+":
                genome_start = start_offset + cand_start + 1  # 1-based
                genome_end   = start_offset + cand_end
            else:
                # On minus strand the working_seq is RC of subseq
                genome_end   = end_offset - cand_start
                genome_start = end_offset - cand_end + 1

            # Junction positions in the candidate sequence (1-based, for reporting)
            # 5' junction G is at arm1_end in working_seq → offset (arm1_end - cand_start) → +1 for 1-based
            junc5_in_cand = ir["arm1_end"] - cand_start + 1      # 1-based pos of 5' G in cand_seq
            junc3_in_cand = ir["arm2_start"] - 1 - cand_start + 1  # 1-based pos of 3' G in cand_seq

            candidates.append({
                "seq_id":        seq_id,
                "strand":        strand,
                "genome_start":  genome_start,
                "genome_end":    genome_end,
                "cand_seq":      cand_seq,
                "cand_len":      cand_len,
                "arm_len":       ir["arm_len"],
                "loop_len":      ir["loop_len"],
                "mismatches":    ir["mismatches"],
                "has_5p_g":      has_5p_g,
                "has_3p_g":      has_3p_g,
                "g5_base":       g5_base,
                "g3_base":       g3_base,
                "junc5_pos":     junc5_in_cand,   # 1-based in cand_seq
                "junc3_pos":     junc3_in_cand,   # 1-based in cand_seq
                "score":         sc,
            })

    # Deduplicate overlapping hits on same (or the opposite) strand, keep highest score.
    # A later hit is discarded if it overlaps at least half of the shorter
    # candidate; this removes many near-duplicates without
    # requiring one hit to be fully contained in the other.
    candidates.sort(key=lambda c: -c["score"])
    kept = []
    for c in candidates:
        overlap = False
        for k in kept:
            #if c["seq_id"] != k["seq_id"] or c["strand"] != k["strand"]:
            if c["seq_id"] != k["seq_id"]: #remove overlapping regions from the same strand as well
                continue

            overlap_start = max(c["genome_start"], k["genome_start"])
            overlap_end = min(c["genome_end"], k["genome_end"])
            overlap_len = max(0, overlap_end - overlap_start + 1)
            shorter_len = min(c["cand_len"], k["cand_len"])

            if shorter_len > 0 and overlap_len >= 0.5 * shorter_len:
                overlap = True
                break

        if not overlap:
            kept.append(c)

    return kept


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def fasta_header(c, idx):
    # Report the two junction bases and whether each is the expected G
    g5 = f"{c['g5_base']}(5'junc,pos={c['junc5_pos']})"
    g3 = f"{c['g3_base']}(3'junc,pos={c['junc3_pos']})"
    both_g = c["has_5p_g"] and c["has_3p_g"]
    g_flag = "both_junc_G=yes" if both_g else (
        "5p_junc_G=yes,3p_junc_G=no" if c["has_5p_g"] else (
        "5p_junc_G=no,3p_junc_G=yes" if c["has_3p_g"] else
        "both_junc_G=no"))
    return (
        f">ncrna{idx} "
        f"seq={c['seq_id']} "
        f"strand={c['strand']} "
        f"start={c['genome_start']} "
        f"end={c['genome_end']} "
        f"len={c['cand_len']} "
        f"arm_len={c['arm_len']} "
        f"loop_len={c['loop_len']} "
        f"mismatches={c['mismatches']} "
        f"5p_junc={g5} "
        f"3p_junc={g3} "
        f"{g_flag} "
        f"score={c['score']}"
    )


def write_fasta(candidates, out_path, wrap=60):
    with open(out_path, "w") as fh:
        for idx, c in enumerate(candidates, 1):
            fh.write(fasta_header(c, idx) + "\n")
            seq = c["cand_seq"]
            for i in range(0, len(seq), wrap):
                fh.write(seq[i : i + wrap] + "\n")
    return len(candidates)


def write_tsv(candidates, out_path):
    """Write a summary TSV alongside the FASTA."""
    cols = [
        "candidate_id", "seq_id", "strand",
        "genome_start", "genome_end", "length",
        "arm_len", "loop_len", "mismatches",
        "5p_junc_base", "5p_junc_pos_in_cand", "5p_junc_is_G",
        "3p_junc_base", "3p_junc_pos_in_cand", "3p_junc_is_G",
        "both_junc_G", "score",
    ]
    with open(out_path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for idx, c in enumerate(candidates, 1):
            row = [
                f"candidate_{idx:04d}",
                c["seq_id"],
                c["strand"],
                c["genome_start"],
                c["genome_end"],
                c["cand_len"],
                c["arm_len"],
                c["loop_len"],
                c["mismatches"],
                c["g5_base"],
                c["junc5_pos"],
                "yes" if c["has_5p_g"] else "no",
                c["g3_base"],
                c["junc3_pos"],
                "yes" if c["has_3p_g"] else "no",
                "yes" if (c["has_5p_g"] and c["has_3p_g"]) else "no",
                c["score"],
            ]
            fh.write("\t".join(str(x) for x in row) + "\n")


def write_gff(candidates, out_path):
    """Write a summary in gff format alongside the FASTA."""
    cols = ["seq_id", "annotator", "region", "start", "end", ".", "strand", "frame", "des"]
    with open(out_path, "w") as fh:
        for idx, c in enumerate(candidates, 1):
            stmp = c["seq_id"].split()
            seq_id = stmp[0]
            des = f"ID={seq_id}_{c["genome_start"]}_{c["genome_end"]}_{c["strand"]};what=ncRNA"
            row = [seq_id, "rscan", "RNA", c["genome_start"], c["genome_end"], '.', c["strand"], '0', des]  
            fh.write("\t".join(str(x) for x in row) + "\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser():
    p = argparse.ArgumentParser(
        prog="retron_ncrna_finder.py",
        description=(
            "Find candidate retron ncRNA (msr-msd) genes using conserved "
            "structural features: an a1/a2 inverted repeat ≥ 8 bp and a "
            "5′ priming guanine."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # --- required ---
    p.add_argument(
        "-i", "--input", required=True, metavar="FASTA",
        help="Input FASTA file (DNA, one or multiple sequences).",
    )

    # --- output ---
    p.add_argument(
        "-o", "--output", default="retron_ncrna_candidates.fasta",
        metavar="FASTA",
        help="Output FASTA file for candidates.",
    )
    p.add_argument(
        "--tsv", default=None, metavar="TSV",
        help="Optional: also write a summary TSV table.",
    )
    p.add_argument(
        "--gff", default=None, metavar="gff",
        help="Optional: also write a summary in gff format.",
    )

    # --- search-region constraints ---
    region = p.add_argument_group("Search-region constraints (1-based, inclusive)")
    region.add_argument(
        "--region-start", type=int, default=None, metavar="INT",
        help="First nucleotide of the region to search (1-based).",
    )
    region.add_argument(
        "--region-end", type=int, default=None, metavar="INT",
        help="Last nucleotide of the region to search (1-based).",
    )

    # --- structural parameters ---
    struct = p.add_argument_group("Structural parameters")
    struct.add_argument(
        "--min-arm", type=int, default=8, metavar="INT",
        help="Minimum inverted-repeat arm length (a1/a2 region).",
    )
    struct.add_argument(
        "--max-loop", type=int, default=400, metavar="INT",
        help="Maximum loop length between the two IR arms.",
    )
    struct.add_argument(
        "--max-mismatch", type=int, default=1, metavar="INT",
        help="Maximum mismatches allowed in the IR arm pairing.",
    )
    struct.add_argument(
        "--min-length", type=int, default=50, metavar="INT",
        help="Minimum total candidate ncRNA length (nt).",
    )
    struct.add_argument(
        "--max-length", type=int, default=500, metavar="INT",
        help="Maximum total candidate ncRNA length (nt).",
    )
    struct.add_argument(
        "--require-g", action="store_true",
        help=(
            "Discard candidates where BOTH inner arm-junction positions "
            "(last base of a1 and first base of a2) are not G. "
            "By default candidates with partial or missing junction Gs are "
            "still reported but scored lower."
        ),
    )

    # --- scoring ---
    p.add_argument(
        "--min-score", type=float, default=0.0, metavar="FLOAT",
        help="Minimum heuristic score to report a candidate.",
    )
    p.add_argument(
        "--top", type=int, default=None, metavar="INT",
        help="Keep only the top N candidates ranked by score.",
    )

    # --- misc ---
    p.add_argument(
        "--both-strands", action="store_true", default=True,
        help="Search both strands (default: True).",
    )
    p.add_argument(
        "--forward-only", action="store_true",
        help="Search the forward strand only.",
    )
    p.add_argument(
        "-v", "--verbose", action="store_true",
        help="Print progress to stderr.",
    )

    return p


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = build_parser()
    args = parser.parse_args()

    if args.forward_only:
        args.both_strands = False

    # Validate
    if args.min_arm < 4:
        parser.error("--min-arm must be ≥ 4.")
    if args.region_start and args.region_end and args.region_start > args.region_end:
        parser.error("--region-start must be ≤ --region-end.")

    if args.verbose:
        print(f"[retron_ncrna_finder] Reading {args.input} …", file=sys.stderr)

    all_candidates = []
    n_seqs = 0

    for seq_id, seq in parse_fasta(args.input):
        n_seqs += 1
        if args.verbose:
            print(f"  Searching {seq_id} ({len(seq):,} nt) …", file=sys.stderr)

        hits = search_sequence(seq_id, seq, args)

        # Score filter
        hits = [h for h in hits if h["score"] >= args.min_score]

        if args.verbose:
            print(f"    → {len(hits)} candidate(s)", file=sys.stderr)

        all_candidates.extend(hits)

    # Global sort + top-N
    all_candidates.sort(key=lambda c: -c["score"])
    if args.top:
        all_candidates = all_candidates[: args.top]

    # Write outputs
    n_written = write_fasta(all_candidates, args.output)

    if args.tsv:
        write_tsv(all_candidates, args.tsv)
        if args.verbose:
            print(f"[retron_ncrna_finder] TSV written to {args.tsv}", file=sys.stderr)

    if args.gff:
        write_gff(all_candidates, args.gff)
        if args.verbose:
            print(f"[retron_ncrna_finder] GFF written to {args.gff}", file=sys.stderr)

    print(
        f"[retron_ncrna_finder] Done. "
        f"Searched {n_seqs} sequence(s). "
        f"Wrote {n_written} candidate(s) to {args.output}.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
