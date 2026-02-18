import numpy as np
import sys

def local_alignment_affine(v, w, match, mismatch, open_penalty, extend_penalty):
    """
    Computes Local Alignment using Affine Gap Penalties (Gotoh's Algorithm variant).
    Uses 3 Matrices:
    - M: Best score ending in a Match/Mismatch
    - X: Best score ending in a Gap in w (Vertical)
    - Y: Best score ending in a Gap in v (Horizontal)
    """
    n, m = len(v), len(w)
    
    # Initialize matrices with -infinity logic, but for Local Alignment we floor at 0
    M = np.zeros((n + 1, m + 1))
    X = np.zeros((n + 1, m + 1))
    Y = np.zeros((n + 1, m + 1))
    
    max_score = -1
    max_pos = (0, 0)
    max_matrix = 'M'
    
    # --- Forward Pass (Fill Matrices) ---
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # 1. Update X (Vertical Gap): Extend existing X or Open new from M
            x_ext = X[i-1, j] + extend_penalty
            x_open = M[i-1, j] + open_penalty
            X[i, j] = max(0, x_ext, x_open) # Floor at 0 for Local
            
            # 2. Update Y (Horizontal Gap): Extend existing Y or Open new from M
            y_ext = Y[i, j-1] + extend_penalty
            y_open = M[i, j-1] + open_penalty
            Y[i, j] = max(0, y_ext, y_open) # Floor at 0 for Local
            
            # 3. Update M (Match/Mismatch)
            score = match if v[i-1] == w[j-1] else mismatch
            m_match = M[i-1, j-1] + score
            m_from_x = X[i-1, j-1] + score
            m_from_y = Y[i-1, j-1] + score
            
            M[i, j] = max(0, m_match, m_from_x, m_from_y)
            
            # Track the global maximum score
            if M[i, j] >= max_score:
                max_score = M[i, j]
                max_pos = (i, j)
                max_matrix = 'M'

    # --- Backward Pass (Traceback) ---
    align_v, align_w = "", ""
    i, j = max_pos
    curr_matrix = max_matrix
    
    # Continue until we hit a zero-score cell (Start of Local Match)
    while i > 0 and j > 0:
        if curr_matrix == 'M':
            if M[i, j] == 0: break # Stop condition for Local Alignment
            
            score = match if v[i-1] == w[j-1] else mismatch
            
            # Did we come from M, X, or Y?
            if M[i, j] == M[i-1, j-1] + score:
                align_v = v[i-1] + align_v
                align_w = w[j-1] + align_w
                i -= 1; j -= 1
                curr_matrix = 'M'
            elif M[i, j] == X[i-1, j-1] + score:
                align_v = v[i-1] + align_v
                align_w = w[j-1] + align_w
                i -= 1; j -= 1
                curr_matrix = 'X'
            elif M[i, j] == Y[i-1, j-1] + score:
                align_v = v[i-1] + align_v
                align_w = w[j-1] + align_w
                i -= 1; j -= 1
                curr_matrix = 'Y'
            else:
                break 

        elif curr_matrix == 'X':
            # Vertical Gap: Did we extend or open?
            if X[i, j] == 0: break
            if X[i, j] == X[i-1, j] + extend_penalty:
                align_v = v[i-1] + align_v
                align_w = "-" + align_w
                i -= 1
                curr_matrix = 'X'
            else: # Opened from M
                align_v = v[i-1] + align_v
                align_w = "-" + align_w
                i -= 1
                curr_matrix = 'M'
                
        elif curr_matrix == 'Y':
            # Horizontal Gap: Did we extend or open?
            if Y[i, j] == 0: break
            if Y[i, j] == Y[i, j-1] + extend_penalty:
                align_v = "-" + align_v
                align_w = w[j-1] + align_w
                j -= 1
                curr_matrix = 'Y'
            else: # Opened from M
                align_v = "-" + align_v
                align_w = w[j-1] + align_w
                j -= 1
                curr_matrix = 'M'
                
    return max_score, align_v, align_w

if __name__ == "__main__":
    
    try:
        seq1 = input("Enter Sequence 1: ").strip()
        seq2 = input("Enter Sequence 2: ").strip()
        
        # User explicitly asked for opening vs extension
        match = float(input("Enter Match Score (e.g. 1): "))
        mismatch = float(input("Enter Mismatch Penalty (e.g. -1): "))
        gap_open = float(input("Enter Gap Opening Penalty (e.g. -2): "))
        gap_extend = float(input("Enter Gap Extension Penalty (e.g. -1): "))
        
        score, a_v, a_w = local_alignment_affine(seq1, seq2, match, mismatch, gap_open, gap_extend)
        
        print("\nBest Alignment:")
        print(f"Score: {score}")
        print(f"Seq1: {a_v}")
        print(f"Seq2: {a_w}")
        
    except ValueError:
        print("Error: Please enter valid numbers for scores.")
