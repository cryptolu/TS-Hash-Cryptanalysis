from binteger import Bin


class Linearization:
    def __init__(self, T):
        self.T = T

        self.N = None
        self.zero = None
        self.basis = None

        self.conds = None

        self.prefix_N = None
        self.prefix_start = None
        self.prefix_step = None
        self.prefix_zero = None
        self.prefix_basis = None

    def linearize_prefix(self, N):
        T = self.T
        self.prefix_N = int(N)

        for s in range(1, 2**32):
            if all(
                (T.alpha * T.h * T.g**(i*s-1)).trace() == 0
                for i in range(1, self.prefix_N)
            ):
                break
        else:
            raise RuntimeError()
        print("found s", s)

        for l in range(1, 2**32):
            if all(
                (T.alpha * T.s0 * T.g**(l+i*s-1)).trace() == 1
                for i in range(0, self.prefix_N)
            ):
                break
        else:
            raise RuntimeError()
        print("found l", l)
        self.init_prefix(start=l, step=s)

    def init_prefix(self, start, step):
        self.prefix_step = step
        self.prefix_start = start

        self.prefix_zero = self.eval_prefix_raw(Bin(0, self.prefix_N))
        self.prefix_basis = [
            self.eval_prefix_raw(Bin.unit(i, self.prefix_N)) ^ self.prefix_zero
            for i in range(self.prefix_N)
        ]

        for _ in range(50):
            x = Bin.random(self.prefix_N)
            h = self.eval_prefix_raw(x)
            h2 = self.prefix_zero
            for i in x.support:
                h2 ^= self.prefix_basis[i]
            assert h == h2

    def encode_prefix(self, msg: Bin) -> Bin:
        x0 = Bin(msg, self.prefix_N)
        xx = []
        lens = [self.prefix_start-1] + [self.prefix_step-1] * (self.prefix_N-1)
        for ln, bit in zip(lens, x0[::-1]):
            xx.extend([0] * ln)
            xx.append(bit)
        return Bin(xx)

    def eval_prefix_raw(self, msg: Bin):
        msg = self.encode_prefix(msg)
        return self.T.fromF(self.T.calc_ext_fast(msg))

    def simulate(self, x: Bin):
        pass
        # h = self.prefix_zero
        # assert x.n == len(basis)
        # for i in x.support:
        #     h ^= basis[i]
        # return h

    def step(self, state: Bin, mbit: int) -> Bin:
        assert mbit in (0, 1)
        if state & 1:
            return (state >> 1) ^ (self.T.p1_raw if mbit else self.T.p0_raw)
        else:
            return state >> 1

    def linearize_quadratic(self, N, tests_per_bit=100):
        self.N = N
        k = self.prefix_N

        zero = self.prefix_zero
        basis = self.prefix_basis
        conds = {}

        for l in range(self.prefix_N, self.N):
            tri = [int(v) & 1 for v in basis] + [int(zero) & 1]
            print("bit", l + 1, "/", self.N)
            if tri[-1]:
                new = self.step(zero, 1) ^ self.step(zero, 0)
            else:
                off = tri.index(1)
                new = self.step(basis[off], 1) ^ self.step(basis[off], 0)

            zero = self.step(zero, 0)
            basis = [self.step(v, 0) for v in basis] + [new]

            conds[l] = tri
            cnt = 0
            for _ in range(tests_per_bit):
                x = Bin.random(l+1)
                pref, post = x[:self.prefix_N], x[self.prefix_N:]
                msg = self.encode_prefix(pref).list + post.list
                h = self.T.fromF(self.T.calc_ext_fast(Bin(msg)))

                h2 = zero
                for i in x.support:
                    h2 ^= basis[i]
                cnt += h == h2
                if h != h2:
                    found = 0
                    for li, cond in conds.items():
                        assert len(cond) == li + 1
                        tr = sum(coef*x[i] for i, coef in enumerate(cond[:li])) + cond[-1]
                        tr &= 1
                        found |= (x[li] == 1 and tr == 0)
                        if found:
                            break
                    assert found
            if tests_per_bit:
                print(cnt, "/", 100, "   nxt", cnt * 3 // 4)
        self.basis = basis
        self.zero = zero
        self.conds = conds
