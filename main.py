from timeit import default_timer as time

__author__ = "7245710, Hoffmann - 7238099, Markmann "
__credits__ = ""
__email__ = "hoffe312@gmail.com - mmaxmarkmann@gmail.com"


# Fasta paths for testing only
#  text_fasta = "C:\\Users\\felix\\Downloads\\AMBIPrak\\Praktikum_1_Data\\text.fasta"
#  virus_fasta = "C:\\Users\\felix\\Downloads\\AMBIPrak\\Praktikum_1_Data\\Virus.fasta"
#  gen_fasta = "C:\\Users\\felix\\Downloads\\AMBIPrak\\Praktikum_1_Data\\gen.fasta"
#  fna_fasta = "C:\\Users\\felix\\Downloads\\AMBIPrak\\Praktikum_1_Data\\BA000002.fna"


# prints all results in formatted matter
def result_print(pattern_matches, successful_shift, name, pattern, exec_time, steps):
    print(f'\n {name}\n Pattern: {pattern}\n matches: {pattern_matches}\n Shifts: {successful_shift}\n time needed: '
          f'{exec_time}s\n steps needed: {steps}\n')


def fasta_reader(fasta_name):
    sequences = ''
    f = open(fasta_name, "r")  # read only
    seq = ""
    bool_flag = True
    while bool_flag:
        moreseq = True
        while moreseq:
            nxtline = f.readline()  # if next line
            if not nxtline:  # break clause
                bool_flag = False
                break
            elif nxtline[0] != ">":  # while there is another line which is not the start of the fasta
                seq += nxtline.strip()  # => append line to empty string
            else:
                moreseq = False  # iterates one line further
    sequences += seq

    return sequences


# matches user choices and hands them over to the pattern matcher of choice
def match_options(algo_user, pattern, choice, fasta_text):
    if choice == 'n':  # if fasta file from path
        match algo_user:
            case 1:
                fh = fasta_reader(fasta_text)
                naive(fh, pattern)
            case 2:
                fh = fasta_reader(fasta_text)
                rabin(fh, pattern)
            case 3:
                fh = fasta_reader(fasta_text)
                knuth(fh, pattern)
            case 4:
                fh = fasta_reader(fasta_text)
                boyer(fh, pattern)
    else:  # own fasta file
        text = input('Your Text: ')
        match algo_user:
            case 1:
                naive(text, pattern)
            case 2:
                rabin(text, pattern)
            case 3:
                knuth(text, pattern)
            case 4:
                boyer(text, pattern)


def naive(text, pattern):
    start = time()  # runtime variable
    name = "NaivePatternMatcher"
    successful_shift = []
    pattern_matches = 0
    steps = 0
    n = len(text)
    m = len(pattern)

    for s in range(n - m):
        count = 0  # counts the matches of pattern[j] and text [s+j]
        j = 0
        while True:
            if text[s + j] == pattern[j] and j <= m:
                j += 1
                count += 1
                steps += 1
            else:  # iterate the text
                steps += 1
                break
            if count == m:  # if match
                pattern_matches += 1
                successful_shift.append(s)
                break

    exec_time = time() - start

    result_print(pattern_matches, successful_shift, name, pattern, exec_time, steps)  # print function


def rabin(text, pattern):
    q = int(input("Primzahl: "))  # modulo
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
            while True:  # while loop from naive pattern matcher
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
            steps += 1
    exec_time = time() - start

    result_print(pattern_matches, successful_shift, name, pattern, exec_time, steps)


def compute_prefix(pattern):
    """prefix algorithm as a part of the Knuth Morris Pratt algorithm"""
    # Longest Proper Prefix that is suffix array
    m = len(pattern)
    pi = [0] * m

    k = 0
    for i in range(2, m + 1):

        while k > 0 and pattern[k] != pattern[i - 1]:
            k = pi[k - 1]

        if pattern[k] == pattern[i - 1]:
            k += 1
            pi[i - 1] = k

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
    for i in range(n):
        while q > 0 and pattern[q] != text[i]:
            q = pi[q - 1]  # if it is not equal => jump
            steps += 1
        if pattern[q] == text[i]:
            q += 1  # if equal => compare next sign
            steps += 1
        if q == m:  # successful => next match
            successful_shift.append((i + 1) - m)
            pattern_matches += 1
            q = pi[q - 1]

    exec_time = time() - start

    result_print(pattern_matches, successful_shift, name, pattern, exec_time, steps)


def last_occurrence(pattern, text):
    """last occurence function as a part of the Boyer Moore algorithm"""
    phi = {}
    for b in text:
        phi[b] = -1
    for a in pattern:
        last = pattern.rfind(a) + 1
        phi[a] = last
    return phi


def good_suffix(pattern, m):
    """good suffix function as a part of the Boyer Moore algorithm"""
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
    phi = last_occurrence(pattern, text)
    gamma = good_suffix(pattern, m)
    s = 0
    steps = 0
    while s <= n - m:
        j = m - 1
        while j >= 0 and pattern[j] == text[s + j]:
            j = j - 1
            steps += 1
        if j == -1:
            successful_shift.append(s)
            pattern_matches += 1
            s = s + gamma[0]
        else:
            s = s + max(gamma[j], j - phi[text[s + j]])
            steps += 1

    exec_time = time() - start
    result_print(pattern_matches, successful_shift, name, pattern, exec_time, steps)


def main():
    pattern_input = input('Pattern as a .fasta = 1 or own pattern = 2: ')
    if pattern_input == '1':
        pattern = input('Pattern.fasta: ')
        pattern = fasta_reader(pattern)
    else:
        pattern = input('Pattern:')
    algo_user = int(input('naive = 1 \nrabin karp = 2\nknuth morris = 3\nboyer moore = 4 :\n'))
    text_choice = input('Own text = y  or fasta data = n :\n')
    fasta_input = input("Fasta path: ")
    match_options(algo_user, pattern, text_choice, fasta_input)


if __name__ == '__main__':
    print("This program searches patterns in fasta files or own input text.\n"
          "You have got the choice out of four different algorithms:\n"
          "1.Naive Pattern Matcher, 2.Rabin Karp Algorithm, 3.Knuth Morris Pratt Algorithm, 4.Boyer Moore Algorithm\n"
          "The program guides you through each step, but it needs a few user inputs like: "
          "pattern (.fasta or typed string), text (.fasta or typed string) and your algorithm choice.")
    while True:
        main()
        print('*' * 20, 'FINISHED', '*' * 20, '\n')
