import time


def result_print(pattern_matches, successful_shift, name, pattern):
    print('Die Funktion', name, 'hat', pattern_matches, 'Matches an den Stellen:', successful_shift,
          'mit dem Muster:', pattern)


def rabin(text, pattern):
    start = time.time()
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
    end = time.time()
    exec_time = end - start
    print(exec_time)
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
    text_fasta = "C:\\Users\\felix\\Downloads\\AMBIPrak\\Praktikum_1_Data\\text.fasta"
    virus_fasta = "C:\\Users\\felix\\Downloads\\AMBIPrak\\Praktikum_1_Data\\Virus.fasta"
    gen_fasta = "C:\\Users\\felix\\Downloads\\AMBIPrak\\Praktikum_1_Data\\gen.fasta"
    pattern = input('Pattern:')
    fh = fasta_reader(text_fasta)
    rabin(fh, pattern)
    naive(fh, pattern)


if __name__ == '__main__':
    main()
