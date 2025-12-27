def longestCommonPrefix(s1, s2):
    i = 0
    while i < len(s1) and i < len(s2) and s1[i] == s2[i]:
        i += 1
    return s1[:i]


def match(s1, s2):
    if not len(s1) == len(s2):
        return False
    for i in range(0, len(s1)):
        if not s1[i] == s2[i]:
            return False
    return True


def reverseComplement(s):
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    t = ""
    for base in s:
        t = complement[base] + t
    return t


############################ Parsing genom
def readGenome(filename):
    genome = ""
    with open(filename, "r") as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == ">":
                genome += line.rstrip()
    return genome


# !wget http://d28rh4a8wq0iu5.cloudfront.net/ads1/data/lambda_virus.fa

############################ Sequences reads


def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip()  # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities


def phred33ToQ(qual):
    return ord(qual) - 33


def createHist(qualityStrings):
    # Create a histogram of quality scores
    hist = [0] * 50
    for read in qualityStrings:
        for phred in read:
            q = phred33ToQ(phred)
            hist[q] += 1
    return hist


######################## Analyzing reads by positions


def findGCByPos(reads):
    """Find the GC ratio at each position in the read"""
    # Keep track of the number of G/C bases and the total number of bases at each position
    gc = [0] * 100
    totals = [0] * 100
    for read in reads:
        for i in range(len(read)):
            if read[i] == "C" or read[i] == "G":
                gc[i] += 1
            totals[i] += 1
    # Divide G/C counts by total counts to get the average at each position
    for i in range(len(gc)):
        if totals[i] > 0:
            gc[i] /= float(totals[i])
    return gc


##### Naive matching


def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):
        match = True
        for j in range(len(p)):
            if t[i + j] != p[j]:
                match = False
                break
        if match:
            occurrences.append(i)
    return occurrences
