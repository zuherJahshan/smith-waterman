import numpy as np
import numba
from numba import njit

import numpy as np
import numba
from numba import njit

@njit
def smith_waterman_numba(seq1, seq2, match_score=2, mismatch_score=-1, gap_penalty=-1):
    m, n = len(seq1), len(seq2)
    score_matrix = np.zeros((m + 1, n + 1), dtype=np.int32)

    max_score = 0
    max_pos = (0, 0)

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score_matrix[i - 1, j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score)
            delete = score_matrix[i - 1, j] + gap_penalty
            insert = score_matrix[i, j - 1] + gap_penalty
            score_matrix[i, j] = max(match, delete, insert, 0)

            if score_matrix[i, j] > max_score:
                max_score = score_matrix[i, j]
                max_pos = (i, j)

    return score_matrix, max_score, max_pos

@njit
def traceback_numba(score_matrix, seq1, seq2, max_pos, match_score=2, mismatch_score=-1, gap_penalty=-1):
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