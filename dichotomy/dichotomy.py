from typing import Callable

class Dichotomy:
    def __init__(self, func: Callable[[float], float], a: float, b: float, eps: float):
        self.func = func
        self.a = a
        self.b = b
        self.eps = eps
    
    def gen(self) -> float:
        while (self.b - self.a) > self.eps:
            x = (self.a + self.b) / 2
            f1 = self.func(x - self.eps)
            f2 = self.func(x + self.eps)
            if f1 < f2:
                self.b = x
            else:
                self.a = x
        return self.b

def gen(func: Callable[[float], float], a: float, b: float, eps: float) -> float:
    return Dichotomy(func, a, b, eps)
