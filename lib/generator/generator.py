import random, sys

class Generator:
    randomizer = random.Random(int(random.uniform(0, sys.maxsize)))
    
    def __init__(self, n, k):
        self.n = n
        self.k = k
    
    def gen(self):
        mx = Generator.randomizer.randint(self.k, sys.maxsize)
        mn = mx // self.k
        matrix = [0] * self.n
        matrix[0] = mn
        matrix[-1] = mx
        for i in range(1, self.n - 1):
            matrix[i] = Generator.randomizer.randint(mn + 1, mx - 1)
        return matrix

def gen(n, k):
    return Generator(n, k).gen()
