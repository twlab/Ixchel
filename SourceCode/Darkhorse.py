import sys


def parse_gfa(filepath):
    """
    Parse the GFA file to extract segments and metadata.
    Handles cases where metadata may be incomplete or malformed and ensures that all
    data is correctly interpreted, especially the SO tags which denote reference start positions.
    """
    segments = {}
    with open(filepath, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if parts[0] == 'S':  # This identifies a segment line
                segment_id = parts[1]
                sequence = parts[2]
                additional_info = {part.split(':')[0]: part.split(':')[-1] for part in parts[3:] if ':' in part}

                # Extract and validate the reference start position
                ref_start_str = additional_info.get('SO', '')
                ref_start = None
                if ref_start_str:
                    try:
                        # Assuming the format SO:i:<number>
                        ref_start = int(ref_start_str.split(':')[2])
                    except (IndexError, ValueError):
                        print(f"Warning: Invalid reference start position '{ref_start_str}' for segment {segment_id}")

                # Determine the orientation based on SR tag
                orientation = '+'
                sr_tag = additional_info.get('SR', '0')
                if sr_tag == '1':
                    orientation = '-'

                segments[segment_id] = {
                    'sequence': sequence,
                    'ref_start': ref_start,
                    'orientation': orientation
                }
    return segments


def find_cytosines(segments):
    """
    Identify cytosine ('C') positions in each segment considering the orientation (+/-).
    Adjusts positions if the segment's orientation is reversed.
    """
    for segment_id, data in segments.items():
        sequence = data['sequence']
        if data['orientation'] == '-':
            sequence = sequence[::-1]  # Reverse the sequence for '-' orientation
        positions = [i for i, base in enumerate(sequence) if base.upper() == 'C']
        if data['orientation'] == '-':
            positions = [len(sequence) - p - 1 for p in positions]  # Convert positions for '-' orientation
        segments[segment_id]['cytosine_positions'] = positions
    return segments


def map_to_hg38(segments):
    """
    Map cytosine positions to hg38 reference coordinates.
    Adds the reference start position to each cytosine position.
    Only segments with valid reference start positions are mapped.
    """
    for segment_id, data in segments.items():
        if data['ref_start'] is not None:
            data['mapped_positions'] = [data['ref_start'] + pos for pos in data['cytosine_positions']]
        else:
            data['mapped_positions'] = []  # No mapping possible if ref_start is None
            print(f"No valid reference start position for segment {segment_id}; cannot map positions.")
    return segments


def main():
    """
    Main function to process the GFA file and map cytosine positions.
    Accepts a file path as a command line argument.
    """
    if len(sys.argv) != 2:
        print("Usage: python ixchel.py <path_to_gfa_file>")
        sys.exit(1)

    filepath = sys.argv[1]
    print("Processing GFA file:", filepath)
    segments = parse_gfa(filepath)
    segments = find_cytosines(segments)
    segments = map_to_hg38(segments)

    # Output the results
    for seg_id, data in segments.items():
        if data['mapped_positions']:
            print(f"Segment {seg_id}: Mapped Cytosine Positions -> {data['mapped_positions']}")
        else:
            print(f"Segment {seg_id}: Mapped Cytosine Positions -> None (Check warnings)")


if __name__ == "__main__":
    main()
