import time

text_fasta = "C:\\Users\\felix\\Downloads\\AMBIPrak\\Praktikum_1_Data\\text.fasta"
virus_fasta = "C:\\Users\\felix\\Downloads\\AMBIPrak\\Praktikum_1_Data\\Virus.fasta"
gen_fasta = "C:\\Users\\felix\\Downloads\\AMBIPrak\\Praktikum_1_Data\\gen.fasta"


def result_print(pattern_matches, successful_shift, name, pattern):
    print('\n', name, '\n Pattern:', pattern, '\n matches:', pattern_matches,
          '\n Shifts:', successful_shift,
          '\n', )


def rabin(text, pattern):
    name = "RabinKarpAlgorithmus"
    pattern_matches = 0
    d = 256
    p = 0  # hash value for pattern
    t = 0  # hash value for txt
    h = 1
    q = int(input("Primzahl:"))  # modulo
    successful_shift = []

    for i in range(len(pattern) - 1):
        h = (h * d) % q

    for i in range(len(pattern)):
        p = (d * p + ord(pattern[i])) % q  # ascii code of letter
        t = (d * t + ord(text[i])) % q

    for s in range(len(text) - len(pattern) + 1):  # move the window one step to the right

        if p == t:
            j = 0
            count = 0
            while True:
                if text[s + j] == pattern[j] and j <= len(pattern):
                    j += 1
                    count += 1
                else:
                    break
                if count == len(pattern):
                    pattern_matches += 1
                    successful_shift.append(s)
                    break
            # Calculate hash value for next window of text: Remove
            # leading digit, add trailing digit
        if s < len(text) - len(pattern):
            t = (d * (t - ord(text[s]) * h) + ord(text[s + len(pattern)])) % q
            if t < 0:
                t = t + q

    result_print(pattern_matches, successful_shift, name, pattern)


def naive(text, pattern):
    start = time.time()
    name = "NaivePatternMatcher"
    successful_shift = []
    pattern_matches = 0
    for s in range(len(text) - len(pattern)):
        count = 0
        j = 0
        while True:
            if text[s + j] == pattern[j] and j <= len(pattern):
                j += 1
                count += 1
            else:
                break
            if count == len(pattern):
                pattern_matches += 1
                successful_shift.append(s)
                break
    end = time.time()
    exec_time = end - start
    print(exec_time)
    result_print(pattern_matches, successful_shift, name, pattern)


def compute_prefix(pattern):
    # Longest Proper Prefix that is suffix array
    M = len(pattern)
    pi = [0] * M

    k = 0
    for i in range(1, M):

        while k > 0 and pattern[k + 1] != pattern[i]:
            k = pi[k]

        if pattern[k + 1] == pattern[i]:
            k += 1
            pi[i] = k

    return pi


def knuth(text, pattern):
    name = "knuth morris"
    successful_shift = []
    pattern_matches = 0
    n = len(text)
    m = len(pattern)
    pi = compute_prefix(pattern)
    q = 0
    for i in range(1, n):
        while q > 0 and pattern[q + 1] != text[i]:
            q = pi[q]
        if pattern[q + 1] == text[i]:
            q += 1
        if q == m:
            successful_shift.append(i)
            pattern_matches += 1
            q = pi[q]
        print(i)
    print(pattern_matches)
    result_print(pattern_matches, successful_shift, name, pattern)


def fasta_reader(fasta_name):
    deflines = []
    sequences = ''
    f = open(fasta_name)
    morelines = True
    while morelines:
        seq = ""
        moreseq = True
        while moreseq:
            nxtline = f.readline()
            if not nxtline:
                morelines = False
                break
            elif nxtline[0] != ">":
                seq += nxtline.strip()
            else:
                deflines.append(nxtline.strip())
                moreseq = False
    sequences += seq

    return sequences


def main():
    path_arr = [text_fasta, virus_fasta, gen_fasta]

    pattern = input('Pattern:')
    f_user = int(input("text.fasta = 1\nvirus.fasta = 2\ngen.fasta = 3\n"))
    algo_user = input("naive = 1 \nrabin karp = 2\nknuthmorris = 3\n")

    match f_user, algo_user:
        case 1, '1':
            fh = fasta_reader(path_arr[0])
            naive(fh, pattern)
        case 1, '2':
            fh = fasta_reader(path_arr[0])
            rabin(fh, pattern)
        case 1, '3':
            fh = fasta_reader(path_arr[0])
            knuth(fh, pattern)
        case 2, '1':
            fh = fasta_reader(path_arr[1])
            naive(fh, pattern)
        case 2, '2':
            fh = fasta_reader(path_arr[1])
            rabin(fh, pattern)
        case 2, '3':
            fh = fasta_reader(path_arr[1])
            knuth(fh, pattern)
        case 3, '1':
            fh = fasta_reader(path_arr[2])
            naive(fh, pattern)
        case 3, '2':
            fh = fasta_reader(path_arr[2])
            rabin(fh, pattern)
        case 3, '3':
            fh = fasta_reader(path_arr[2])
            knuth(fh, pattern)


if __name__ == '__main__':
    main()
