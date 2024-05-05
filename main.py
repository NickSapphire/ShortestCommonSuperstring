from shortest_common_superstring import *

if __name__ == "__main__":
    seq = input("Enter the sequence (leave blank for a randomly generated one): ")
    if seq == "":
        size = int(input("Enter the length of the random sequence: "))
        seq = generate_genome(size)
        print("Sequence generated:", seq)

    k = int(input("Enter the length of the k-mers: "))
    if k > len(seq):
        raise ValueError("k cannot be higher than the length of the sequence")

    threshold = int(input("Enter the minimum overlap length (< k): "))
    if k <= threshold:
        raise ValueError("Overlap threshold must be smaller than k")

    print("\nGreedy approach:")
    greedy = find_SCS_greedy(*create_overlap_matrix(seq, k, threshold))
    print(f"{greedy}\nLength: {len(greedy)}")

    do_bf = input("\nTest brute force approach? (WARNING: Very long computation "
                  "time for long sequences) [y/n]: ").lower()
    if do_bf == "y":
        bf = find_SCS_bf(seq, k)
        print("Brute force approach:")
        print(f"{bf}\nLength: {len(bf)}")
