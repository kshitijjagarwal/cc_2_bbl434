import numpy as np
import sys

def global_alignment(v, w, match_score, mismatch_penalty, gap_penalty):
    """
    Computes the Global Alignment (Needleman-Wunsch).
    """
    n = len(v)
    m = len(w)
    
    # Initialize the scoring matrix with gap penalties in first row/col
    s = np.zeros((n + 1, m + 1), dtype=int)
    
    for i in range(1, n + 1):
        s[i, 0] = i * gap_penalty
    for j in range(1, m + 1):
        s[0, j] = j * gap_penalty
        
    # Fill the matrix (Forward Pass)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if v[i - 1] == w[j - 1]:
                diagonal = s[i - 1, j - 1] + match_score
            else:
                diagonal = s[i - 1, j - 1] + mismatch_penalty
                
            vertical = s[i - 1, j] + gap_penalty
            horizontal = s[i, j - 1] + gap_penalty
            
            s[i, j] = max(diagonal, vertical, horizontal)

    # Backtrack (Reverse Pass)
    align_v = ""
    align_w = ""
    i, j = n, m
    
    while i > 0 or j > 0:
        if i > 0 and j > 0 and (v[i-1] == w[j-1]):
             if s[i, j] == s[i-1, j-1] + match_score:
                align_v = v[i-1] + align_v
                align_w = w[j-1] + align_w
                i -= 1; j -= 1
                continue
        
        if i > 0 and j > 0 and (v[i-1] != w[j-1]):
            if s[i, j] == s[i-1, j-1] + mismatch_penalty:
                align_v = v[i-1] + align_v
                align_w = w[j-1] + align_w
                i -= 1; j -= 1
                continue
                
        if i > 0 and s[i, j] == s[i-1, j] + gap_penalty:
            align_v = v[i-1] + align_v
            align_w = "-" + align_w
            i -= 1
        else:
            align_v = "-" + align_v
            align_w = w[j-1] + align_w
            j -= 1
            
    return s[n, m], align_v, align_w

if __name__ == "__main__":
    try:
        seq1 = input("Enter Sequence 1: ").strip()
        seq2 = input("Enter Sequence 2: ").strip()
        
        match = int(input("Enter Match Score (e.g. 1): "))
        mismatch = int(input("Enter Mismatch Penalty (e.g. -1): "))
        gap = int(input("Enter Gap Penalty (e.g. -2): "))
        
        score, res_v, res_w = global_alignment(seq1, seq2, match, mismatch, gap)
        
        print(f"\nAlignment Score: {score}")
        print(f"Seq1: {res_v}")
        print(f"Seq2: {res_w}")
        
    except ValueError:
        print("Error: Please enter valid integers for scores.")
