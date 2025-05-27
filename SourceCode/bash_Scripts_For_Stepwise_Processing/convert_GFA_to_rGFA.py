#!/usr/bin/env python3
import argparse
import sys
import re


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert a GFA with potentially fragmented reference walks (W-lines) into rGFA by annotating reference segments"
    )
    parser.add_argument(
        "gfa",
        help="Path to input GFA file"
    )
    parser.add_argument(
        "-r", "--ref-walk",
        required=True,
        help="Name of the W-line(s) to use as reference segments"
    )
    parser.add_argument(
        "-o", "--out",
        default=None,
        help="Path to output rGFA file (defaults to stdout)"
    )
    return parser.parse_args()


def read_segment_lengths(gfa_path):
    """
    First pass: read all S-lines to map segment ID → sequence length.
    """
    print("[1/4] Reading segment lengths...")
    seglens = {}
    with open(gfa_path) as f:
        for line in f:
            if line.startswith("S\t"):
                parts = line.rstrip("\n").split("\t")
                sid = parts[1]
                seq = parts[2]
                seglens[sid] = len(seq)
    print(f"    Found {len(seglens)} segments.")
    return seglens


def read_ref_fragments(gfa_path, ref_walk):
    """
    Second pass: collect all W-lines matching ref_walk name, extracting
    tuples of (start_coord, segment_ids_list, end_coord, hap_index, contig).
    """
    print(f"[2/4] Collecting reference walk fragments for '{ref_walk}'...")
    fragments = []
    with open(gfa_path) as f:
        for line in f:
            if not line.startswith("W\t"):
                continue
            parts = line.rstrip("\n").split("\t")
            name = parts[1]
            if name != ref_walk:
                continue
            if len(parts) < 7:
                sys.exit(f"ERROR: Malformed W-line for '{ref_walk}', expected ≥7 columns, got {len(parts)}")
            hap_index = parts[2]
            contig = parts[3]
            try:
                start = int(parts[4])
                end = int(parts[5])
            except ValueError:
                sys.exit(f"ERROR: W-line for '{ref_walk}' has non-integer start/end: '{parts[4]}','{parts[5]}'")
            mapping_field = parts[6]
            seg_ids = re.findall(r'[><]([^><]+)', mapping_field)
            if not seg_ids:
                sys.exit(f"ERROR: No segment entries found in mapping '{mapping_field}' for '{ref_walk}'")
            fragments.append((start, seg_ids, end, hap_index, contig))
    if not fragments:
        sys.exit(f"ERROR: No W-lines named '{ref_walk}' found in GFA")
    print(f"    Found {len(fragments)} fragments.")
    return fragments


def sanity_check_fragments(fragments, seglens, ref_walk):
    """
    Third pass: for each fragment, assert start + sum(lengths) == end.
    """
    print(f"[3/4] Performing sanity checks on {len(fragments)} fragments...")
    for start, seg_ids, end, hap_index, contig in fragments:
        length_sum = 0
        for sid in seg_ids:
            if sid not in seglens:
                sys.exit(f"ERROR: Segment '{sid}' in ref_walk fragment not found among S-lines")
            length_sum += seglens[sid]
        if start + length_sum != end:
            sys.exit(
                f"ERROR: Fragment sanity failed for '{ref_walk}' (hap{hap_index} on {contig}): "
                f"start({start}) + sum(lengths={length_sum}) != end({end})"
            )
    print("    Sanity checks passed.")


def assemble_ref_path(fragments):
    """
    Fourth pass: sort fragments by start_coord and build full ref_path
    as list of (segment_id, hap_index, contig).
    """
    fragments.sort(key=lambda x: x[0])
    ref_path = []
    for start, seg_ids, end, hap_index, contig in fragments:
        for sid in seg_ids:
            ref_path.append((sid, hap_index, contig))
    print(f"[4/4] Assembled reference path of {len(ref_path)} segments.")
    return ref_path


def compute_offsets(ref_path, seglens):
    """
    After assembly: compute 0-based cumulative offsets for each reference segment,
    mapping segment_id → (offset, hap_index, contig).
    """
    print("Computing offsets for reference segments...")
    offsets = {}
    cum = 0
    for sid, hap_index, contig in ref_path:
        offsets[sid] = (cum, hap_index, contig)
        cum += seglens[sid]
    print(f"    Computed offsets for {len(offsets)} segments.")
    return offsets


def convert_to_rgfa(gfa_path, offsets, ref_walk, out_path=None):
    """
    Final pass: stream through GFA, appending SN, SO, SR tags to reference S-lines,
    and emitting all other lines unchanged.
    """
    target = out_path or 'stdout'
    print(f"Writing rGFA output to {target}...")
    out_f = open(out_path, "w") if out_path else sys.stdout
    with open(gfa_path) as f:
        for line in f:
            if line.startswith("S\t"):
                parts = line.rstrip("\n").split("\t")
                sid = parts[1]
                if sid in offsets:
                    offset, hap_index, contig = offsets[sid]
                    tags = [
                        f"SN:Z:{ref_walk}#{hap_index}#{contig}",
                        f"SO:i:{offset}",
                        "SR:i:0"
                    ]
                    out_f.write("\t".join(parts + tags) + "\n")
                else:
                    out_f.write(line)
            else:
                out_f.write(line)
    if out_path:
        out_f.close()
    print("rGFA generation complete.")


def main():
    args = parse_args()
    seglens = read_segment_lengths(args.gfa)
    fragments = read_ref_fragments(args.gfa, args.ref_walk)
    sanity_check_fragments(fragments, seglens, args.ref_walk)
    ref_path = assemble_ref_path(fragments)
    offsets = compute_offsets(ref_path, seglens)
    convert_to_rgfa(args.gfa, offsets, args.ref_walk, args.out)


if __name__ == "__main__":
    main()