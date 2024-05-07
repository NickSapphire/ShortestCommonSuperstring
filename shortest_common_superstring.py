import random
import itertools

def overlap(S: str, P: str) -> int:
    """
    Given two strings S and P, returns the length of the largest overlap of a suffix
    of S and a prefix of P.
    """
    overlap_size = 0
    for i in range(1, min(len(S), len(P)) + 1):  # Tries all suffix-prefix overlaps of P and S
        if S[-i:] == P[:i]:  # Slice S from the end and P from the start
            overlap_size = i
    return overlap_size

def generate_genome(length: int) -> str:
    """
    Generates a string of a given length with evenly distributed random nucleotides.
    """
    if length < 0:
        raise ValueError("length must be a positive integer")
    seq = "".join(random.choices(population="ACGT", k=length))
    return seq


def get_substrings(seq: str, k: int) -> list:
    """
    Returns the list of overlapping k-mers from the string seq.
    """
    kmers = []
    for i in range(len(seq) - k + 1):
        kmers.append(seq[i:i+k])
    return kmers


def create_overlap_matrix(seq: str, k: int, threshold: int) -> tuple:
    """
    Given a sequence and the size of the k-mers initializes the overlap matrix.
    Explained further in the README file.
    :param seq: Input sequence
    :param k: Size of the overlapping substrings
    :param threshold: Minimum overlap length required to get added to the overlap matrix
    :return: A tuple of the overlap matrix and the substrings dictionary
    """
    size = len(seq) - k + 1
    IDs = [(i,) for i in range(size)]  # Store k-mer IDs as tuples for merging
    matrix = {ID:{} for ID in IDs}     # Each inner dict stores overlaps with ID as the suffix

    kmers = get_substrings(seq, k)
    substrings = {(i,):kmers[i] for i in range(size)}  # Keys: IDs from matrix,
                                                       # Values: Corresponding substrings
    for s in range(len(kmers)):      # Suffix IDs
        for p in range(len(kmers)):  # Prefix IDs
            if s != p:  # Exclude alignments of a k-mer with itself
                overlap_length = overlap(kmers[s], kmers[p])
                if overlap_length >= threshold:  # Exclude overlaps smaller than threshold
                    matrix[(s,)][(p,)] = overlap_length
    return (matrix, substrings)


def merge(matrix: dict, substrings: dict, suf_id: tuple, pre_id: tuple):
    """
    Performs in-place modification on matrix and substrings, merging the overlapping
    strings with IDs suf_id and pre_id.
    :param matrix: Overlap matrix (generated from create_overlap_matrix)
    :param substrings: Dictionary mapping tuple IDs to substrings of the sequence
                       (generated in create_overlap_matrix)
    :param suf_id: ID of the substring whose suffix is taken for the overlap
    :param pre_id: ID of the substring whose prefix is taken for the overlap
    """
    # Retrieve overlap and merge the strings
    overlap_size = matrix[suf_id][pre_id]
    merged_string = substrings[suf_id] + substrings[pre_id][overlap_size:]

    # Remove prefix and suffix strings from substrings and matrix (in place)
    del substrings[suf_id]
    del substrings[pre_id]
    del matrix[suf_id]
    del matrix[pre_id]
    for dic in matrix.values():
        dic.pop(suf_id, None)  # Will be ignored if the key is not found
        dic.pop(pre_id, None)

    # Add the new substring ID to substrings and matrix, and compute the overlap lengths
    merged_ID = suf_id + pre_id
    substrings[merged_ID] = merged_string

    new_row = dict()
    for id, dic in matrix.copy().items():  # doesn't let you modify a dict in a loop
        matrix[id][merged_ID] = overlap(substrings[id], merged_string)  # merged_string as prefix
        new_row[id] = overlap(merged_string, substrings[id]) # merged_string as suffix
    matrix[merged_ID] = new_row


def find_SCS_greedy(matrix: dict, substrings: dict) -> str:
    """
    Finds the shortest common superstring from an array of overlapping k-mers
    by recursively merging the two substrings with the highest overlap length
    until it gets to one substring. In case more than one pair has the largest
    overlap, it picks a random pair among them.
    """
    if len(substrings) == 1:
        return [*substrings.values()][0]

    # Find the maximum overlap and all the positions at which it occurs
    max_overlap = max(i for dic in matrix.values() for i in dic.values())
    max_positions = [(s, p) for s, dic in matrix.items() for p, val in dic.items()
                     if val == max_overlap]

    # Select a random (suffix, prefix) pair from max_position and merge them
    s, p = random.choice(max_positions)
    merge(matrix, substrings, s, p)
    return find_SCS_greedy(matrix, substrings)


def find_SCS_bf(seq: str, k: int) -> str:
    """
    Iterates over all possible permutations of a given string's k-mers and finds
    a common superstring with the smallest possible length.
    """
    substrings = get_substrings(seq, k)
    min_length = float("inf")  # initialized with infinity since it's a minimization problem

    for perm in itertools.permutations(substrings):
        merged_string = ""
        for substring in perm:
            overlap_size = overlap(merged_string, substring)
            merged_string += substring[overlap_size:]
        if len(merged_string) < min_length:
            scs = merged_string
            min_length = len(merged_string)
    return scs
