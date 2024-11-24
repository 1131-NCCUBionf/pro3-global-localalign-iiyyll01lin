import numpy as np
from Bio import SeqIO


def parse_score_matrix(score_path):
    with open(score_path, "r") as f:
        lines = [
            line.strip() for line in f if not line.startswith("#") and line.strip()
        ]
    headers = lines[0].split()
    matrix = {char: {} for char in headers}
    for line in lines[1:]:
        parts = line.split()
        row_char = parts[0]
        for i, score in enumerate(parts[1:]):
            matrix[row_char][headers[i]] = int(score)
    return matrix


def initialize_matrix(rows, cols, global_alignment, gap):
    score_matrix = np.zeros((rows, cols), dtype=int)
    traceback_matrix = np.zeros((rows, cols), dtype=int)
    if global_alignment:
        for i in range(rows):
            score_matrix[i][0] = gap * i
            traceback_matrix[i][0] = 1  # Up
        for j in range(cols):
            score_matrix[0][j] = gap * j
            traceback_matrix[0][j] = 2  # Left
    return score_matrix, traceback_matrix


def alignment(input_path, score_path, output_path, aln, gap):
    records = list(SeqIO.parse(input_path, "fasta"))
    seq1_id, seq1 = records[0].id, str(records[0].seq)
    seq2_id, seq2 = records[1].id, str(records[1].seq)

    # print(seq1)
    # print(seq2)

    score_matrix = parse_score_matrix(score_path)
    rows, cols = len(seq1) + 1, len(seq2) + 1
    global_alignment = aln == "global"
    score, traceback = initialize_matrix(rows, cols, global_alignment, gap)

    max_score = float('-inf')
    # print(max_score)
    max_positions = []

    for i in range(1, rows):
        for j in range(1, cols):
            match_score = score_matrix[seq1[i - 1]].get(seq2[j - 1], gap)
            match = score[i - 1, j - 1] + match_score
            delete = score[i - 1, j] + gap
            insert = score[i, j - 1] + gap

            if not global_alignment:
                # print(global_alignment)
                # Reset score to 0 if all options are negative
                current_score = max(0, match, delete, insert)
                print(current_score)
                score[i, j] = current_score

                # Track positions with max score
                if current_score > max_score:
                    max_score = current_score
                    max_positions = [(i, j)]
                elif current_score == max_score and current_score > 0:
                    # Only add if it's not overlapping with existing positions
                    is_distinct = True
                    for pos_i, pos_j in max_positions:
                        if abs(i - pos_i) < 6 and abs(j - pos_j) < 6:  # Adjust threshold as needed
                            is_distinct = False
                            break
                    if is_distinct:
                        max_positions.append((i, j))
            else:
                score[i, j] = max(match, delete, insert)

            # Update traceback matrix
            if score[i, j] == 0 and not global_alignment:
                traceback[i, j] = 3
            elif score[i, j] == match:
                traceback[i, j] = 0
            elif score[i, j] == delete:
                traceback[i, j] = 1
            else:
                traceback[i, j] = 2

    if global_alignment:
        # Handle global alignment as before
        i, j = len(seq1), len(seq2)
        aligned_seq1, aligned_seq2 = traceback_global(seq1, seq2, i, j, traceback)
        with open(output_path, "w") as f:
            f.write(f">{seq1_id}\n{aligned_seq1}\n")
            f.write(f">{seq2_id}\n{aligned_seq2}\n")
    else:
        # Find all local alignments
        alignments = set()
        
        for end_i, end_j in max_positions:
            print(end_i, end_j)
            curr_i, curr_j = end_i, end_j
            temp_seq1, temp_seq2 = "", ""

            # Trace back until hitting a zero or the start
            while curr_i > 0 and curr_j > 0 and score[curr_i, curr_j] > 0:
                if traceback[curr_i, curr_j] == 0:  # Diagonal
                    temp_seq1 = seq1[curr_i - 1] + temp_seq1
                    temp_seq2 = seq2[curr_j - 1] + temp_seq2
                    curr_i -= 1
                    curr_j -= 1
                elif traceback[curr_i, curr_j] == 1:  # Up
                    curr_i -= 1
                elif traceback[curr_i, curr_j] == 2:  # Left
                    curr_j -= 1
                else:  # Stop
                    break

            if temp_seq1 and temp_seq2:
                start_i = curr_i
                start_j = curr_j
                aligned_seq1 = seq1[start_i:end_i]
                aligned_seq2 = seq2[start_j:end_j]

                # print(seq1, seq2)
                # print(aligned_seq1, aligned_seq2)

                alignments.add((
                    aligned_seq1,
                    aligned_seq2,
                    # seq1, 
                    # seq2,
                    len(aligned_seq1),
                    seq1.index(aligned_seq1),
                    seq2.index(aligned_seq2)
                ))

        alignments = list(alignments)
        # print(alignments)
        max_len = max(aln[2] for aln in alignments)
        max_alignments = [aln for aln in alignments if aln[2] == max_len]
        max_alignments.sort(key=lambda x: (x[3], x[4]))

        print("\nAll maximum length alignments:")
        for aligned_seq1, aligned_seq2, length, pos1, pos2 in max_alignments:
            print(f"\nAlignment (length {length}):")
            print(f"Sequence 1 (position {pos1}): {aligned_seq1}")
            print(f"Sequence 2 (position {pos2}): {aligned_seq2}")

        with open(output_path, "w") as f:
            for aligned_seq1, aligned_seq2, _, _, _ in max_alignments:
                f.write(f">{seq1_id}\n{aligned_seq1}\n")
                f.write(f">{seq2_id}\n{aligned_seq2}\n")


def traceback_global(seq1, seq2, i, j, traceback):
    aligned_seq1, aligned_seq2 = "", ""
    while i > 0 or j > 0:
        if traceback[i, j] == 0:  # Diagonal
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            i -= 1
            j -= 1
        elif traceback[i, j] == 1:  # Up
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        else:  # Left
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            j -= 1
    return aligned_seq1, aligned_seq2

