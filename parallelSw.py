import numpy as np
from concurrent.futures import ThreadPoolExecutor

def smith_waterman_parallel(seq1, seq2, match_score=2, mismatch_score=-1, gap_penalty=-1, num_threads=8):
    m, n = len(seq1), len(seq2)
    score_matrix = np.zeros((m + 1, n + 1), dtype=int)
    max_score = 0
    max_pos = (0, 0)

    def fill_score_matrix(i):
        nonlocal max_score, max_pos
        for j in range(1, n + 1):
            score = score_matrix[i - 1, j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score)
            score = max(0, score, score_matrix[i - 1, j] + gap_penalty, score_matrix[i, j - 1] + gap_penalty)
            score_matrix[i, j] = score
            if score > max_score:
                max_score = score
                max_pos = (i, j)
        return max_score, max_pos

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        results = list(executor.map(fill_score_matrix, range(1, m + 1)))

    max_score, max_pos = max(results, key=lambda x: x[0])

    return score_matrix, max_score, max_pos


def traceback_parallel(score_matrix, seq1, seq2, max_pos, match_score=2, mismatch_score=-1, gap_penalty=-1):
    i, j = max_pos
    align1, align2 = "", ""

    while i > 0 and j > 0 and score_matrix[i, j] != 0:
        current_score = score_matrix[i, j]
        diagonal_score = score_matrix[i - 1, j - 1]
        up_score = score_matrix[i - 1, j]
        left_score = score_matrix[i, j - 1]

        if current_score == diagonal_score + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score):
            align1 = seq1[i - 1] + align1
            align2 = seq2[j - 1] + align2
            i -= 1
            j -= 1
        elif current_score == up_score + gap_penalty:
            align1 = seq1[i - 1] + align1
            align2 = "-" + align2
            i -= 1
        else:
            align1 = "-" + align1
            align2 = seq2[j - 1] + align2
            j -= 1

    return align1, align2
