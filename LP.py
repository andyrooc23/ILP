import numpy.linalg
from decimal import Decimal, getcontext
import json


def solve(b, a, s, ct):
    try:
        getcontext().prec = 10
        n = len(b)
        ls = [0] * n
        for i in range(n):
            ls[i] = float(Decimal(b[i]) - Decimal(s[i]))
        xvals = numpy.linalg.solve(a, ls)

        '''
                Another way to solve this system manually
                ainv = numpy.linalg.inv(a)
                multiply left side by inverse of a to find values of x and return
                xvals = numpy.dot(ainv, ls)
                '''
        if isinstance(xvals[0], list):
            biggest = xvals[0]
            for j in range(1, len(xvals)):
                if numpy.dot(ct, xvals[j]) > numpy.dot(ct, biggest):
                    biggest = xvals[j]
            print("Verifying solution")
            print("True" if (numpy.dot(a, xvals) is b) else "False")
            xvals = biggest

        print("Verifying solution")
        ax = numpy.dot(a, xvals).tolist()
        ax = [round(num) for num in ax]
        print(ax, "  is the dot product of solution with original matrix")
        print(ls, "  is the calculated b - s")
        print("Valid solution" if ax == ls else "Invalid Solution")
        return xvals

    except ValueError:
        print("Opps! The matrix isn't solvable")


def main():
    file = json.load(open('matrix3.json'))
    matrices = file['Matrix']
    a = matrices['a']
    b = matrices['b']
    s = matrices['s']
    ct = matrices['ct']
    answers = solve(b, a, s, ct)
    out = "Maximum value answers is "
    for num in answers:
        out += str(round(num, 5)) + "   "
    print(out)
    print(str(round(numpy.dot(ct, answers), 5)) + "   is the value of cT and X dot product")


if __name__ == '__main__':
    main()
