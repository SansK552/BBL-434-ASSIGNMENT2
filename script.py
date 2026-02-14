import sys
import math

def read_fasta(filename):
    with open(filename) as f:
        lines = f.readlines()
    seq = ""
    for line in lines:
        if not line.startswith(">"):
            seq += line.strip()
    return seq

def affine_global_alignment(seq1, seq2, match, mismatch, gap_open, gap_extend):

    n = len(seq1)
    m = len(seq2)

    M  = [[-math.inf]*(m+1) for _ in range(n+1)]
    Ix = [[-math.inf]*(m+1) for _ in range(n+1)]
    Iy = [[-math.inf]*(m+1) for _ in range(n+1)]

    M[0][0] = 0

    for i in range(1, n+1):
        Ix[i][0] = gap_open + (i-1)*gap_extend
        M[i][0] = Ix[i][0]

    for j in range(1, m+1):
        Iy[0][j] = gap_open + (j-1)*gap_extend
        M[0][j] = Iy[0][j]

    for i in range(1, n+1):
        for j in range(1, m+1):

            score = match if seq1[i-1] == seq2[j-1] else mismatch

            M[i][j] = max(
                M[i-1][j-1],
                Ix[i-1][j-1],
                Iy[i-1][j-1]
            ) + score

            Ix[i][j] = max(
                M[i-1][j] + gap_open,
                Ix[i-1][j] + gap_extend
            )

            Iy[i][j] = max(
                M[i][j-1] + gap_open,
                Iy[i][j-1] + gap_extend
            )

    i, j = n, m
    matrices = {"M": M[n][m], "Ix": Ix[n][m], "Iy": Iy[n][m]}
    current = max(matrices, key=matrices.get)

    align1, align2 = "", ""

    while i > 0 or j > 0:

        if current == "M":
            score = match if seq1[i-1] == seq2[j-1] else mismatch
            if M[i][j] == M[i-1][j-1] + score:
                current = "M"
            elif M[i][j] == Ix[i-1][j-1] + score:
                current = "Ix"
            else:
                current = "Iy"

            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1

        elif current == "Ix":
            if Ix[i][j] == M[i-1][j] + gap_open:
                current = "M"
            else:
                current = "Ix"

            align1 = seq1[i-1] + align1
            align2 = "-" + align2
            i -= 1

        else:
            if Iy[i][j] == M[i][j-1] + gap_open:
                current = "M"
            else:
                current = "Iy"

            align1 = "-" + align1
            align2 = seq2[j-1] + align2
            j -= 1

    final_score = max(M[n][m], Ix[n][m], Iy[n][m])

    return align1, align2, final_score


if __name__ == "__main__":

    seq1 = read_fasta(sys.argv[1])
    seq2 = read_fasta(sys.argv[2])

    match = int(sys.argv[3])
    mismatch = int(sys.argv[4])
    gap_open = int(sys.argv[5])
    gap_extend = int(sys.argv[6])

    align1, align2, score = affine_global_alignment(
        seq1, seq2, match, mismatch, gap_open, gap_extend
    )

    with open("alignment_output.txt", "w") as f:
        f.write("Score: " + str(score) + "\n\n")
        f.write(align1 + "\n")
        f.write(align2 + "\n")

    print("Alignment complete. Output saved to alignment_output.txt")
