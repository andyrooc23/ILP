import numpy.linalg
import json


def solve(b, a, s, ct):
    try:
        n = len(b)
        ls = [0] * n
        for i in range(n):
            ls[i] = b[i] - s[i]
        print(ls)
        xvals = numpy.linalg.solve(a, ls)


        # ainv = numpy.linalg.inv(a)
        # multiply left side by inverse of a to find values of x and return
        # xvals = numpy.dot(ainv, ls)
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
            print("True" if (numpy.dot(a, xvals) is b) else "False")
            return biggest


    except ValueError:
        print("Opps! The matrix isn't singular")

    #return numpy.dot(ls, ct)


def main():
    # b = [5.0, 12.9, 2.0]
    #
    # a = [[3.2, 8.7, 5.9],
    #      [2.4, 3.1, 1.1],
    #      [9.7, 6.1, 0.3]]
    #
    # s = [3.7, 3.8, 1.0]
    #
    # ct = [8.2, 9.7, 1.1]
    file = json.load(open('matrix2.json'))
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
