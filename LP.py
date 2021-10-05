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

        if not isinstance(xvals[0], list):
            print("Verifying solution")
            print("True" if (numpy.dot(a, xvals) is b) else "False")
            return xvals

        else:
            biggest = xvals[0]
            for j in range(1, len(xvals)):
                if numpy.dot(ct, xvals[j]) > numpy.dot(ct, biggest):
                    biggest = xvals[j]
            print("Verifying solution")
            print("True" if (numpy.dot(a, xvals) is b)  else "False")
            return biggest


    except ValueError:
        print("Opps! The matrix isn't singular")

    #return numpy.dot(ls, ct)


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
        out += str(num) + "   "
    print(out)
    print(str(numpy.dot(ct, answers)) + "   is the value of cT.X")


if __name__ == '__main__':
    main()
