from sage.all import GF, vector, matrix
import random
import itertools

from binteger import Bin


class TSHash:
    def __init__(self, n, seed=None):
        if seed is None:
            seed = random.randrange(2**128)
        state0 = random.getstate()
        random.seed(seed)
        self.seed = seed

        self.n = int(n)
        self.p0_raw = self._gen_poly()
        self.p1_raw = self._gen_poly()
        # set LSB = coef of x^(n-1) to avoid initial skips
        self.s0_raw = Bin.random(n) | 1

        self.h_raw = self.p0_raw ^ self.p1_raw

        X = GF(2).polynomial_ring().gen()
        modulus = X**self.n + sum(X**(i) for i in self.p0_raw.support)
        self.F = GF(2**self.n, name='X', modulus=modulus)

        self.g = self.F.gen()
        self.h = self.toF(self.h_raw)

        self.s0 = self.toF(self.s0_raw)
        self.alpha = self._get_alpha()

        random.setstate(state0)

        self._last_pos = None

    def _gen_poly(self):
        while True:
            p = Bin.random(self.n) | (1 << (self.n-1))  # MSB = coef of x^0
            x = GF(2)['x'].gen()
            pol = sum(x**(i) for i in p.support) + x**self.n
            if pol.is_irreducible():
                return p

    def _get_alpha(self):
        n = self.n
        mat = matrix(GF(2), [
            [(self.g**(i+j)).trace() for i in range(n)]
            for j in range(n)
        ])
        alpha = mat.solve_left(vector([0] * (n-1) + [1]))
        alpha = sum(take * self.g**i for i, take in enumerate(alpha))
        assert [(alpha * self.g**i).trace() for i in range(n)] == [0] * (n-1) + [1]
        return alpha

    def toF(self, v: Bin) -> "F2n":
        v = Bin(v, self.n)[::-1].int
        return self.F.from_integer(v)

    def fromF(self, v: "F2n") -> Bin:
        v = v.to_integer()
        return Bin(v, self.n)[::-1]

    def calc_raw(self, msg: Bin) -> Bin:
        s = self.s0_raw
        poly = [self.p0_raw, self.p1_raw]
        for mi in msg:
            while s & 1 == 0:
                s >>= 1
            s = (s >> 1) ^ poly[mi]
        return s

    def calc_ext(self, emsg: Bin) -> "F2n":
        s = self.s0
        for mi in Bin(emsg):
            s = self.g * s + mi * self.h * (self.alpha * s).trace()
        return s

    def calc_ext_fast(self, emsg: Bin) -> "F2n":
        s = self.s0_raw
        for mi in Bin(emsg):
            if s & 1:
                s = (s >> 1) ^ (self.p1_raw if mi else self.p0_raw)
            else:
                s = (s >> 1)
        return self.toF(s)

    def calc(self, msg: Bin) -> "F2n":
        s = self.s0
        for mi in Bin(msg):
            while not (self.alpha * s).trace():
                s = self.g * s
            s = self.g * s + mi * self.h
        return s

    def extend(self, msg: Bin) -> Bin:
        s = self.s0
        emsg = []
        self._last_pos = []
        for mi in Bin(msg):
            while not (self.alpha * s).trace():
                s = self.g * s
                emsg.append(0)
            s = self.g * s + mi * self.h
            self._last_pos.append(len(emsg))
            emsg.append(mi)
        return Bin(emsg)

    def compress(self, emsg: Bin) -> Bin:
        s = self.s0
        msg = []
        self._last_pos = []
        for imi, mi in enumerate(Bin(emsg)):
            if (self.alpha * s).trace():
                s = self.g * s + mi * self.h
                msg.append(mi)
                self._last_pos.append(imi)
            else:
                s = self.g * s
        return Bin(msg)

    def tr(self, v: "F2n") -> "F2":
        """same as (self.alpha * v).trace()"""
        return v.polynomial()[self.n-1]
