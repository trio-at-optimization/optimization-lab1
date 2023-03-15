import random, sys

def gen(n, k):
    randomizer = random.Random(int(random.uniform(0, sys.maxsize)))
    mx = randomizer.randint(k, sys.maxsize)
    mn = mx // k
    matrix = [0] * n
    matrix[0] = mn
    matrix[-1] = mx
    for i in range(1, n - 1):
        matrix[i] = randomizer.randint(mn + 1, mx - 1)
    return matrix
