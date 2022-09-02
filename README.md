# The Nussinov Algorithm for RNA secondary structure prediction.

## Idea (biological)
- Stacked base pairs of helical regions are considered to stabilize an RNA molecule.
→ maximize the number of base pairs.
- Input: RNA sequence S
- Output: a non-crossing RNA structure P of S that maximizes
- |P| (i.e. the number of base pairs in P).

## Nussinov considering:
- i j pair
- i being unpaired
- j being unpair
- bifurcation

S(i, j) = max (S(i + 1, j - 1) + 1 [if i, j base pair], S(i + 1, j), S(i, j + 1), max S(i, k) + S(K + 1, j))
