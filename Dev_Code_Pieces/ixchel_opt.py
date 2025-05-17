#!/usr/bin/env python3
import argparse
import sqlite3
import os
import struct
import itertools


def build_db(precomputed_file, db_file, batch_size):
    """Bulk-load the .converted file into SQLite with tuned pragmas."""
    print(f"[build_db] opening {db_file}")
    conn = sqlite3.connect(db_file)
    # speed-tuning pragmas
    conn.execute("PRAGMA synchronous = OFF;")
    conn.execute("PRAGMA journal_mode = MEMORY;")
    conn.execute("PRAGMA temp_store = MEMORY;")
    conn.execute("PRAGMA locking_mode = EXCLUSIVE;")
    conn.execute("BEGIN;")

    insert_sql = """
      INSERT OR REPLACE INTO conversion
        (segment_id, segment_offset, stable_source, start, stop, conversion_code)
      VALUES (?, ?, ?, ?, ?, ?)
    """
    batch = []
    count = 0

    with open(precomputed_file, 'r') as fin:
        for line in fin:
            fields = line.rstrip("\n").split("\t")
            # map 'NA' → None
            start = None if fields[1] == 'NA' else int(fields[1])
            stop = None if fields[2] == 'NA' else int(fields[2])
            code = None if fields[7] == 'NA' else int(fields[7])
            batch.append((fields[8], fields[9], fields[0], start, stop, code))
            count += 1

            if count % batch_size == 0:
                conn.executemany(insert_sql, batch)
                conn.commit()
                print(f"  inserted {count} rows so far")
                batch.clear()

    # final leftovers
    if batch:
        conn.executemany(insert_sql, batch)
        conn.commit()
        print(f"  inserted total {count} rows")
    conn.close()
    print("[build_db] complete")


def convert_methyl_optimized(graph_methyl, db_file, output_file):
    """Stream through .graph.methyl, group by segment, and do in-RAM lookups."""
    print(f"[convert] opening DB {db_file} in WAL mode")
    conn = sqlite3.connect(db_file, check_same_thread=False)
    conn.execute("PRAGMA journal_mode = WAL;")
    cursor = conn.cursor()

    def flush_segment(seg_id, lines, fout):
        # pull every offset→mapping for this segment_id
        cursor.execute(
            "SELECT segment_offset, stable_source, start, stop, conversion_code "
            "FROM conversion WHERE segment_id = ?", (seg_id,)
        )
        mapping = {offset: (src, s, e, c) for offset, src, s, e, c in cursor.fetchall()}
        for L in lines:
            F = L.rstrip("\n").split("\t")
            offset = F[1]
            src, s, e, c = mapping.get(offset, ("NA", "NA", "NA", "NA"))
            out = [
                src,  # stable_source
                s,  # start
                e,  # stop
                F[3],  # context
                F[7],  # methylatedFraction
                F[2],  # strand
                F[6],  # coverage
                c,  # conversionCode
                F[0],  # segmentID
                offset  # segmentOffset
            ]
            fout.write("\t".join(map(str, out)) + "\n")

    with open(graph_methyl, 'r') as fin, open(output_file, 'w') as fout:
        current_seg = None
        buffer = []
        for line in fin:
            seg = line.split("\t", 1)[0]
            if seg != current_seg:
                if buffer:
                    flush_segment(current_seg, buffer, fout)
                current_seg = seg
                buffer = [line]
            else:
                buffer.append(line)
        # last segment
        if buffer:
            flush_segment(current_seg, buffer, fout)

    conn.close()
    print(f"[convert] done, output → {output_file}")


def main():
    parser = argparse.ArgumentParser(prog="ixchel_opt", description="Ixchel w/ optimized DB routines")
    sub = parser.add_subparsers(dest='cmd')

    p1 = sub.add_parser('build_db', help='build optimized SQLite DB')
    p1.add_argument('precomputed_file')
    p1.add_argument('db_file')
    p1.add_argument('--batch', type=int, default=5_000_000,
                    help='rows per insert batch')

    p2 = sub.add_parser('convert_methyl_optimized',
                        help='convert .graph.methyl using optimized lookups')
    p2.add_argument('graph_methyl')
    p2.add_argument('db_file')
    p2.add_argument('output_file')

    args = parser.parse_args()
    if args.cmd == 'build_db':
        # ensure table exists
        conn = sqlite3.connect(args.db_file)
        conn.execute("""
                     CREATE TABLE IF NOT EXISTS conversion
                     (
                         segment_id
                         TEXT,
                         segment_offset
                         TEXT,
                         stable_source
                         TEXT,
                         start
                         INTEGER,
                         stop
                         INTEGER,
                         conversion_code
                         INTEGER,
                         PRIMARY
                         KEY
                     (
                         segment_id,
                         segment_offset
                     )
                         )""")
        conn.close()
        build_db(args.precomputed_file, args.db_file, args.batch)
    elif args.cmd == 'convert_methyl_optimized':
        convert_methyl_optimized(args.graph_methyl, args.db_file, args.output_file)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
