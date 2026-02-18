import numpy as np
import sys

def longest_common_subsequence(v, w):
    """
    Computes the Longest Common Subsequence (LCS) between strings v and w.
    """
    n = len(v)
    m = len(w)
    
    # Initialize the scoring matrix
    s = np.zeros((n + 1, m + 1), dtype=int)
    
    # Fill the matrix (Forward Pass)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if v[i - 1] == w[j - 1]:
                s[i, j] = s[i - 1, j - 1] + 1
            else:
                s[i, j] = max(s[i - 1, j], s[i, j - 1])

    # Backtrack (Reverse Pass)
    lcs_str = ""
    i, j = n, m
    
    while i > 0 and j > 0:
        if v[i - 1] == w[j - 1]:
            lcs_str = v[i - 1] + lcs_str
            i -= 1
            j -= 1
        elif s[i - 1, j] > s[i, j - 1]:
            i -= 1
        else:
            j -= 1
            
    return s[n, m], lcs_str

if __name__ == "__main__":
    # Clean User Interface
    try:
        seq1 = input("Enter Sequence 1: ").strip()
        seq2 = input("Enter Sequence 2: ").strip()
        
        score, alignment = longest_common_subsequence(seq1, seq2)
        
        print(f"\nLCS Length: {score}")
        print(f"LCS String: {alignment}")
        
    except Exception as e:
        print(f"An error occurred: {e}")
