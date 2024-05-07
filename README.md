# ShortestCommonSuperstring

This is an implementation of the Shortest Common Superstring problem with greedy
and brute force approaches.

## Library dependencies

- random: https://docs.python.org/3/library/random.html

- itertools: https://docs.python.org/3/library/itertools.html

This program uses the random module to generate random sequences with an even
distribution and to decide which substring pair to merge in the greedy function,
and itertools to iterate over every permutation of the k-mer array in the brute
force solution.

## How to use

The shortest_common_superstring script contains all the functions which are imported into
the main script for testing with user input or randomly generated sequences.

Once it's run, the program will prompt the user for the following input, respectively:
- ```seq```: The original sequence the program will partition into overlapping k-mers and try to
reconstruct. If it is left blank, it will ask for another argument ```size``` and create
a randomly generated sequence of that length.
- ```k```: Length of the overlapping substrings ```seq``` is broken into.
- ```threshold```: Minimum suffix-prefix overlap length required for a pair of substrings to
get added to the overlap matrix.

Parameter restrictions: ```len(seq)``` >= ```k``` > ```threshold```

Afterwards, it will output the greedy solution to the shortest common superstring problem for the
generated k-mers. The function picks randomly from the substring pairs with maximum length
in each recursive call, so it may return different results each time it's run.

Finally, since the brute force approach is very inefficient, even for relatively short sequences,
I have decided to leave it up to the user whether to test it or not to avoid having to stop
the program every time.


## The overlap matrix in detail

Each k-mer is assigned a unique ID and it is stored within the substrings dictionary.
The IDs are stored as tuples of one value to keep track of merging, e.g. if substrings with
IDs (1,3) and (2,4) respectively are merged, the new substring will have the ID (1,3,2,4),
obtained by simply summing the tuples.

The overlap matrix is a dictionary which has tuples of k-mer IDs as keys and
inner dictionaries containing the size of the overlaps of the outer key (as suffix)
and the inner keys (as prefixes) as values. None of the inner dictionaries contain
overlaps of any given substring with itself.


## Contact

Nicolas Schiappa

nicolas.schiappa@studio.unibo.it