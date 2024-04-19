import random
import string


def generate_random_sequence(length):
    """Generate a random sequence of given length."""
    return ''.join(random.choices('ACTG', k=length))


def create_gfa_file(num_nodes, max_sequence_length, num_edges, filename="TestGraph.gfa"):
    """
    Generate a random GFA file with the specified number of nodes and edges.

    Parameters:
    - num_nodes: int, number of sequences (nodes) in the graph
    - max_sequence_length: int, maximum length of genome sequences
    - num_edges: int, number of edges between the sequences
    - filename: str, the name of the file to save the GFA content
    """
    with open(filename, 'w') as file:
        file.write("H\tVN:Z:1.0\n")  # Header line for GFA version 1.0
        nodes = []

        # Generate nodes with random sequences
        for i in range(1, num_nodes + 1):
            seq_length = random.randint(10, max_sequence_length)  # Sequence length between 10 and max_sequence_length
            sequence = generate_random_sequence(seq_length)
            file.write(f"S\t{i}\t{sequence}\n")
            nodes.append(i)

        # Generate random edges between nodes
        for _ in range(num_edges):
            from_node = random.choice(nodes)
            to_node = random.choice(nodes)
            if from_node != to_node:  # Avoid self-loops
                overlap_length = random.randint(1, 10)  # Overlap length between 1 and 10
                file.write(f"L\t{from_node}\t+\t{to_node}\t+\t{overlap_length}M\n")


if __name__ == "__main__":
    NUM_NODES = 10  # Example: 10 nodes
    MAX_SEQUENCE_LENGTH = 100  # Example: maximum sequence length of 100
    NUM_EDGES = 15  # Example: 15 edges

    create_gfa_file(NUM_NODES, MAX_SEQUENCE_LENGTH, NUM_EDGES, "RandomTestGraph.gfa")
