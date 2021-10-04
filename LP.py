import numpy.linalg
from decimal import Decimal, getcontext


def solve(b, a, s, ct):
    getcontext().prec = 10
    ls = [0 for i in range(len(b))]
    for i in range(len(b)):
        ls[i] = float(Decimal(b[i]) - Decimal(s[i]))


def main():
    b = [5.0, 12.9, 2.0]
    a = [[3.2, 8.7, 5.9],
         [2.4, 3.1, 1.1],
         [9.7, 6.1, 0.3]]
    s = [3.7, 3.8, 1.0]
    ct = [8.2, 9.7, 1.1]
    print(solve(b, a, s, ct))


if __name__ == '__main__':
    main()
