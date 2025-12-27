import sys

filename = sys.argv[1]


def read_fasta_file_and_process(
    filename: str,
    n_seq: int,
):
    start_seq = "ATG"
    ends_seq = ("TAA", "TAG", "TGA")

    seqs = {}
    seqs_len = {}
    codons_seq = {}
    codons_stat = {}
    repeat_seq = {}

    try:
        with open(filename, "r") as file:
            nb_records = 0
            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    nb_records += 1
                    words = line.split()
                    name = words[0][1:]
                    # On garde ici le nom de la séquence
                    seqs[name] = ""
                    seqs_len[name] = 0
                else:
                    # On prend la séquence
                    seqs[name] = seqs[name] + line
                    seqs_len[name] = seqs_len[name] + len(line)

        max_value = max(seqs_len.values())
        min_value = min(seqs_len.values())

        max_keys = [k for k, v in seqs_len.items() if v == max_value]
        min_keys = [k for k, v in seqs_len.items() if v == min_value]

        ######## Identification des parties codantes
        for k, v in seqs.items():
            codons_seq[k] = ""
            codons_stat[k] = ""
            # codons = set()
            codons = {}
            for pos in [0, 1, 2]:
                seq = v[pos:]
                for i in range(0, len(seq), 3):
                    if seq[i : i + 3] == start_seq:
                        for j in range(i + 3, len(seq), 3):
                            if seq[j : j + 3] in ends_seq:
                                codons[i + 1] = ""
                                codons[i + 1] = [seq[i : j + 3], len(seq[i : j + 3])]
                                break

            if codons:
                max_value = max(v[1] for v in codons.values())
                max_keys = [k for k, v in codons.items() if v[1] == max_value]
                codons_stat[k] = [max_keys, max_value]
            else:
                codons_stat[k] = [[], 0]

        all_codons = {}
        for k1, v1 in codons_seq.items():
            for k2, v2 in v1:
                if k2 in all_codons.keys():
                    all_codons[k2] = all_codons[k2] + v2
                else:
                    all_codons[k2] = v2
        if all_codons:
            max_codons_len = max([v for v in all_codons.values()])
            max_codons_keys = [v for v in all_codons.values() if v == max_codons_len]
        else:
            max_codons_len = None
            max_codons_keys = []

        ##################### Repeat seq
        if n_seq > 1:
            for k, v in seqs.items():
                repeat_seq[k] = ""
                res = {}
                seq = v

                for i in range(0, len(seq) - n_seq + 1):
                    fixed_seq = seq[i : i + n_seq]
                    res[fixed_seq] = ""
                    n = 0
                    for j in range(0, len(seq), n_seq):
                        current_seq = seq[j : j + n_seq]
                        if fixed_seq == current_seq:
                            n = n + 1
                    res[fixed_seq] = n

                repeat_seq[k] = res

            all_repeat_seq = {}
            if repeat_seq:
                for k1, v1 in repeat_seq.items():
                    for k2, v2 in v1.items():
                        if k2 in all_repeat_seq.keys():
                            all_repeat_seq[k2] = all_repeat_seq[k2] + v2
                        else:
                            all_repeat_seq[k2] = v2
            if all_repeat_seq:
                max_repeat_seq = max([v for v in all_repeat_seq.values()])
                max_repeat_seq_keys = [
                    v for v in all_repeat_seq.values() if v == max_repeat_seq
                ]
            else:
                max_repeat_seq = None
                max_repeat_seq_keys = []

        res = {
            "nb_records": nb_records,
            "seq": seqs,
            "seq_len": seqs_len,
            "codons_seq": codons_seq,
            "all_codons": all_codons,
            "repeat_seq": repeat_seq,
            "all_repeat_seq": all_repeat_seq,
            "stat": {
                "min_value": min_value,
                "max_value": max_value,
                "max_keys": max_keys,
                "min_keys": min_keys,
                "codons_stat": codons_stat,
                "max_codons_len": max_codons_len,
                "max_codons_keys": max_codons_keys,
                "max_repeat_seq": max_repeat_seq,
                "max_repeat_seq_keys": max_repeat_seq_keys,
            },
        }
        return res

    except IOError:
        print("Le fichier n'existe pas")
        return None
