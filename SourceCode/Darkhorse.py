import sys


def parse_gfa(filepath):
    """
    Parse the GFA file to extract segments and metadata.
    This function extracts segment information including IDs, sequences, and
    reference start positions while considering orientation from GFA formatted text.
    Handles cases where metadata may be incomplete or malformed.
    """
    segments = {}
    with open(filepath, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if parts[0] == 'S':  # This identifies a segment line
                segment_id = parts[1]
                sequence = parts[2]
                additional_info = {part.split(':')[0]: part.split(':')[-1] for part in parts[3:] if ':' in part}
                ref_start = None
                if 'SO' in additional_info:
                    try:
                        ref_start = int(additional_info['SO'][2:])  # Correctly handle the integer conversion
                    except ValueError:
                        print(f"Warning: Invalid reference start position for segment {segment_id}")

                orientation = '+'  # Default orientation
                if 'SR' in additional_info and additional_info['SR'] == '1':
                    orientation = '-'  # Reverse orientation if specified

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
    """
    for segment_id, data in segments.items():
        if data['ref_start'] is not None:
            data['mapped_positions'] = [data['ref_start'] + pos for pos in data['cytosine_positions']]
        else:
            data['mapped_positions'] = []
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
        print(f"Segment {seg_id}: Mapped Cytosine Positions -> {data.get('mapped_positions', [])}")


if __name__ == "__main__":
    main()
