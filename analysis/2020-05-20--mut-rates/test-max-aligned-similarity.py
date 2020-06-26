
def max_aligned_similarity(a, b):
    if len(a) != len(b):
        print("Must be equal length!")
        exit(-1)
    # perfectly aligned
    num_bits = len(a)
    max_similarity = -1
    # right_shift_alignments = {}
    # left_shift_alignments = {}
    for shift in range(0, num_bits):
        right_alignment_similarity = 0
        left_alignment_similarity = 0

        for i in range(0, len(a)-shift):
            right_alignment_similarity += int(a[i+shift] == b[i])
        # right_shift_alignments[shift] = right_alignment_similarity

        if shift == 0:
            max_similarity = right_alignment_similarity
            continue

        for i in range(0, len(b)-shift):
            left_alignment_similarity += int(a[i] == b[i+shift])

        max_similarity = max(max_similarity, right_alignment_similarity, left_alignment_similarity)
        # left_shift_alignments[shift] = left_alignment_similarity
    # print("Right:", right_shift_alignments)
    # print("Left:", left_shift_alignments)
    print(f"Max alignment similarity score for {a} and {b}:", max_similarity)
    return max_similarity

max_aligned_similarity("0011", "0000")
max_aligned_similarity("1111", "0000")
max_aligned_similarity("1111", "1111")
max_aligned_similarity("0101", "1010")
max_aligned_similarity("1100", "1010")
max_aligned_similarity("1111", "0010")