import sys

def read_sequences(file1, file2):
    with open(file1, 'r') as f:
        seq1 = f.readline().strip()
    with open(file2, 'r') as f:
        seq2 = f.readline().strip()
    return seq1, seq2

def get_score_matrix(match_score=1, mismatch_score=-1):
    score_matrix = {}
    for base1 in 'ACGT':
        for base2 in 'ACGT':
            if base1 == base2:
                score_matrix[(base1, base2)] = match_score
            else:
                score_matrix[(base1, base2)] = mismatch_score
    return score_matrix

def init_matrix(seq1_len, seq2_len, gap_score=-1):
    matrix = [[0 for _ in range(seq2_len+1)] for _ in range(seq1_len+1)]
    for i in range(1, seq1_len+1):
        matrix[i][0] = matrix[i-1][0] + gap_score
    for j in range(1, seq2_len+1):
        matrix[0][j] = matrix[0][j-1] + gap_score
    return matrix

def fill_matrix(seq1, seq2, score_matrix, gap_score=-1):
    seq1_len, seq2_len = len(seq1), len(seq2)
    matrix = init_matrix(seq1_len, seq2_len, gap_score)
    for i in range(1, seq1_len+1):
        for j in range(1, seq2_len+1):
            match_score = matrix[i-1][j-1] + score_matrix[(seq1[i-1], seq2[j-1])]
            delete_score = matrix[i-1][j] + gap_score
            insert_score = matrix[i][j-1] + gap_score
            matrix[i][j] = max(match_score, delete_score, insert_score)
    return matrix

def traceback(matrix, seq1, seq2, score_matrix, gap_score=-1):
    seq1_len, seq2_len = len(seq1), len(seq2)
    aligned_seq1, aligned_seq2 = '', ''
    i, j = seq1_len, seq2_len
    while i > 0 or j > 0:
        if i > 0 and j > 0 and matrix[i][j] == matrix[i-1][j-1] + score_matrix[(seq1[i-1], seq2[j-1])]:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        elif i > 0 and matrix[i][j] == matrix[i-1][j] + gap_score:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = '-' + aligned_seq2
            i -= 1
        else:
            aligned_seq1 = '-' + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1
    return aligned_seq1, aligned_seq2

def main(file1, file2):
    seq1, seq2 = read_sequences(file1, file2)
    score_matrix = get_score_matrix()
    matrix = fill_matrix(seq1, seq2, score_matrix)
    aligned_seq1, aligned_seq2 = traceback(matrix, seq1, seq2, score_matrix)
    print(f"Aligned Sequence 1: {aligned_seq1}")
    print(f"Aligned Sequence 2: {aligned_seq2}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python pairwise_alignment.py <sequence1_file> <sequence2_file>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])