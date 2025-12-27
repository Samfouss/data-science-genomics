def identify_codon(seq: str, n_seq: int):
    res = {}
    for i in range(0, len(seq) - n_seq + 1):
        fixed_seq = seq[i : i + n_seq]
        res[fixed_seq] = None
        n = 0
        for j in range(0, len(seq), n_seq):
            current_seq = seq[j : j + n_seq]
            if fixed_seq == current_seq:
                n = n + 1
        res[fixed_seq] = n

    return res


print(identify_codon("AAAAAAAAAAATG", 3))

codons = identify_codon("AAAAAAAAAAATG", 3)
max_keys = [k for k, v in codons.items() if v == 3]
print(max_keys)
