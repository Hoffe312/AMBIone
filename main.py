from timeit import default_timer as time

text_fasta = "C:\\Users\\felix\\Downloads\\AMBIPrak\\Praktikum_1_Data\\text.fasta"
virus_fasta = "C:\\Users\\felix\\Downloads\\AMBIPrak\\Praktikum_1_Data\\Virus.fasta"
gen_fasta = "C:\\Users\\felix\\Downloads\\AMBIPrak\\Praktikum_1_Data\\gen.fasta"
fna_fasta = "C:\\Users\\felix\\Downloads\\AMBIPrak\\Praktikum_1_Data\\BA000002.fna"


def result_print(pattern_matches, successful_shift, name, pattern, exec_time, steps):
    print('\n', name, '\n Pattern:', pattern, '\n matches:', pattern_matches,
          '\n Shifts:', successful_shift,
          '\n time needed:', exec_time, 's\n steps needed:', steps)


def fasta_reader(fasta_name):
    deflines = []
    sequences = ''
    f = open(fasta_name)
    bool_flag = True
    while bool_flag:
        seq = ""
        moreseq = True
        while moreseq:
            nxtline = f.readline()
            if not nxtline:
                bool_flag = False
                break
            elif nxtline[0] != ">":
                seq += nxtline.strip()
            else:
                deflines.append(nxtline.strip())
                moreseq = False
    sequences += seq

    return sequences


def match_options(algo_user, pattern, choice):  # if fasta is used
    fasta_input = text_fasta
    if choice == 'n':
        match algo_user:
            case '1':
                fh = fasta_reader(fasta_input)
                naive(fh, pattern)
            case '2':
                fh = fasta_reader(fasta_input)
                rabin(fh, pattern)
            case '3':
                fh = fasta_reader(fasta_input)
                knuth(fh, pattern)
            case '4':
                fh = fasta_reader(fasta_input)
                boyer(fh, pattern)
    else:
        text = input('Your Text:')
        match algo_user:
            case '1':
                naive(text, pattern)
            case '2':
                rabin(text, pattern)
            case '3':
                knuth(text, pattern)
            case '4':
                boyer(text, pattern)


def naive(text, pattern):
    start = time()
    name = "NaivePatternMatcher"
    successful_shift = []
    pattern_matches = 0
    steps = 0
    n = len(text)
    m = len(pattern)

    for s in range(n - m):
        count = 0
        j = 0
        while True:
            if text[s + j] == pattern[j] and j <= m:
                j += 1
                count += 1
                steps += 1
            else:
                steps += 1
                break
            if count == m:
                pattern_matches += 1
                successful_shift.append(s)
                break

    exec_time = time() - start

    result_print(pattern_matches, successful_shift, name, pattern, exec_time, steps)


def rabin(text, pattern):
    q = int(input("Primzahl:"))  # modulo
    start = time()
    name = "RabinKarpAlgorithm"
    pattern_matches = 0
    d = 256
    p = 0  # hash value for pattern
    t = 0  # hash value for txt
    h = 1
    n = len(text)
    m = len(pattern)
    steps = 0
    successful_shift = []

    for i in range(m - 1):
        h = (h * d) % q

    for i in range(m):
        p = (d * p + ord(pattern[i])) % q  # ascii code of letter
        t = (d * t + ord(text[i])) % q

    for s in range(n - m + 1):  # move the window one step to the right
        if p == t:
            j = 0
            count = 0
            while True:
                if text[s + j] == pattern[j] and j <= m:
                    j += 1
                    count += 1
                    steps += 1
                else:
                    steps += 1
                    break
                if count == m:
                    pattern_matches += 1
                    successful_shift.append(s)
                    break
            # Calculate hash value for next window of text: Remove
            # leading digit, add trailing digit
        if s < n - m:
            t = (d * (t - ord(text[s]) * h) + ord(text[s + m])) % q
            if t < 0:
                t = t + q

    exec_time = time() - start

    result_print(pattern_matches, successful_shift, name, pattern, exec_time, steps)


def compute_prefix(pattern):
    # Longest Proper Prefix that is suffix array
    m = len(pattern)
    pi = [0] * m

    k = 0
    for i in range(1, m):

        while k > 0 and pattern[k + 1] != pattern[i]:
            k = pi[k - 1]

        if pattern[k + 1] == pattern[i]:
            k += 1
            pi[i] = k

    return pi


def knuth(text, pattern):
    start = time()
    name = "KnuthMorrisAlgorithm"
    successful_shift = []  # array of successful shifts
    pattern_matches = 0  # pattern match counter
    n = len(text)  # length of text
    m = len(pattern)  # length of pattern
    pi = compute_prefix(pattern)
    q = 0  # count of similarities
    steps = 0
    for i in range(1, n):
        while q > 0 and pattern[q] != text[i]:
            q = pi[q - 1]  # if it is not equal => jump
            steps += 1
        if pattern[q] == text[i]:
            q += 1  # if equal => compare next sign
            steps += 1
        if q == m:  # successful => next match
            successful_shift.append(i + 1)
            pattern_matches += 1
            q = pi[q - 1]
            steps += 1
    exec_time = time() - start

    result_print(pattern_matches, successful_shift, name, pattern, exec_time, steps)


def last_occurence(pattern, text):
    phi = {}
    for b in text:
        phi[b] = -1
    for a in pattern:
        last = pattern.rfind(a) + 1
        phi[a] = last
    return phi


def good_suffix(pattern, m):
    gamma = [0] * m
    pi = compute_prefix(pattern)
    pi_reverse = compute_prefix(pi[::-1])
    for j in range(m):
        gamma[j] = m - pi[m - 1]
    for l in range(m):
        j = m - pi_reverse[l] - 1
        if gamma[j] > l - pi_reverse[l]:
            gamma[j] = l + 1 - pi_reverse[l]
    return gamma


def boyer(text, pattern):
    start = time()
    name = 'BoyerMooreAlgorithm'
    successful_shift = []
    pattern_matches = 0

    n = len(text)
    m = len(pattern)
    phi = last_occurence(pattern, text)
    gamma = good_suffix(pattern, m)
    s = 0
    steps = 0
    while s <= n - m:
        j = m-1
        while j >= 0 and pattern[j] == text[s+j]:
            j = j - 1
            steps += 1
        if j == -1:
            successful_shift.append(s)
            pattern_matches += 1
            s = s + gamma[0]
            steps += 1
        else:
            s = s + max(gamma[j], j - phi[text[s+j]])
            steps += 1

    exec_time = time() - start
    result_print(pattern_matches, successful_shift, name, pattern, exec_time, steps)


def main():
    pattern = input('Pattern:')
    algo_user = input('naive = 1 \nrabin karp = 2\nknuth morris = 3\nboyer moore = 4\n')
    text_choice = input('Own text = y  or fasta data = n :\n')

    match_options(algo_user, pattern, text_choice)


if __name__ == '__main__':
    while True:
        main()
        print('\nNew patternmatch:\n')
