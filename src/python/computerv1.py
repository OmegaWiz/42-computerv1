import sys

# UTILS

def fmt(x):
    return f"{x+0.0:g}"

def ft_sqrt(x):
    if x != x or x == float('inf'):
        raise ValueError("sqrt of NaN or inf")
    if x < 0:
        raise ValueError("sqrt of negative")
    if x == 0:
        return 0.0
    g = x if x >= 1 else 1.0
    while True:
        nxt = (g + x / g) / 2
        if nxt >= g:
            return g
        g = nxt

MAX_POWER = 1000
EPS = 1e-12

def is_finite(x):
    return x - x == 0

def is_near_zero(x, scale):
    return abs(x) <= EPS * scale

def is_valid_expo(exp):
    # plain ASCII digits only: rejects "-1", "+2", "1_0", non-ASCII digits
    if not (exp.isascii() and exp.isdigit()):
        return False
    return int(exp) <= MAX_POWER

def is_valid_coeff(coeff):
    # optional sign, ASCII digits, at most one '.', at least one digit:
    # rejects "nan", "inf", "1e3", "1_000", non-ASCII digits
    body = coeff[1:] if coeff[:1] in ('-', '+') else coeff
    if not body.isascii() or body.count('.') > 1:
        return False
    digits = body.replace('.', '')
    return digits != '' and digits.isdigit()

class Polynomial:
    def reduce(self):
        while self.coeffs and self.coeffs[-1] == 0:
            self.coeffs.pop()
        if not self.coeffs:
            self.coeffs = [0]

    def __init__(self, src = None):
        self.coeffs = [0]
        if src is None:
            pass
        elif isinstance(src, list):
            self.coeffs = list(src)
        elif isinstance(src, str):
            args = src.strip().split()
            i = 0
            while i < len(args):
                if (args[i] == '-'):
                    sign = -1
                else:
                    sign = 1
                if (not args[i] in ('-', '+')):
                    if (i != 0):
                        raise ValueError("Invalid polynomial string")
                else:
                    i += 1
                if (i >= len(args)):
                    raise ValueError("Invalid polynomial string")

                if (args[i].startswith(('X', 'x'))):
                    coeff = 1.0
                else:
                    if (not is_valid_coeff(args[i])):
                        raise ValueError("Invalid polynomial string")
                    coeff = float(args[i])
                    if (not is_finite(coeff)):
                        raise ValueError("Coefficient out of range")
                    if (i + 1 >= len(args) or args[i + 1] in ('-', '+')):
                        self.add_term(sign * coeff, 0)
                        i += 1
                        continue
                    i += 1
                    if (i >= len(args)):
                        raise ValueError("Invalid polynomial string")
                    if (args[i] != '*'):
                        raise ValueError("Invalid polynomial string")
                    i += 1
                    if (i >= len(args)):
                        raise ValueError("Invalid polynomial string")
                    if (not args[i].startswith(('X', 'x'))):
                        raise ValueError("Invalid polynomial string")
                if (args[i] == 'X' or args[i] == 'x'):
                    power = 1
                elif (args[i].startswith(('X^', 'x^'))):
                    if (not is_valid_expo(args[i][2:])):
                        raise ValueError("Invalid exponent")
                    power = int(args[i][2:])
                else:
                    raise ValueError("Invalid polynomial string")
                self.add_term(sign * coeff, power)
                i += 1

        else:
            raise TypeError("Invalid type for Polynomial initialization")
        self.reduce()

    def degree(self):
        return len(self.coeffs) - 1

    def __str__(self):
        terms = []
        for i in range(len(self.coeffs)):
            coeff = self.coeffs[i]
            if (i > 0):
                if (coeff >= 0):
                    term = f"+ {fmt(coeff)} * X^{i}"
                else:
                    term = f"- {fmt(abs(coeff))} * X^{i}"
            else:
                term = f"{fmt(coeff)} * X^{i}"
            terms.append(term)
        return ' '.join(terms) if terms else '0'

    def __sub__(self, other):
        max_len = max(len(self.coeffs), len(other.coeffs))
        self.coeffs.extend([0] * (max_len - len(self.coeffs)))
        other.coeffs.extend([0] * (max_len - len(other.coeffs)))
        res = [0] * max_len
        for i in range(max_len):
            l, r = self.coeffs[i], other.coeffs[i]
            res[i] = l - r
            if is_near_zero(res[i], max(abs(l), abs(r))):
                res[i] = 0
        return Polynomial(res)

    def add_term(self, coeff, power):
        self.coeffs.extend([0] * (power + 1 - len(self.coeffs)))
        old = self.coeffs[power]
        new = old + coeff
        if is_near_zero(new, max(abs(old), abs(coeff))):
            new = 0
        self.coeffs[power] = new
        self.reduce()

# main
def main():
    if len(sys.argv) > 2:
        print(f'usage: python3 {sys.argv[0]} "<equation>"')
        sys.exit(1)
    elif len(sys.argv) < 2:
        try:
            arg  = input("Enter a polynomial equation in the form of sums of 'n * X^k' (e.g., '5 * X^0 + 4 * X^1 - 9.3 * X^2 = 1 * X^0'):\n")
        except (EOFError, KeyboardInterrupt):
            sys.exit(1)
    else:
        arg = sys.argv[1]
    if (arg.count('=') != 1):
        print("Error: Invalid argument, must contain exactly one '='")
        sys.exit(1)
    l, r = arg.split('=')
    if (l.strip() == '' or r.strip() == ''):
        print("Error: Invalid argument, must contain both sides of the equation")
        sys.exit(1)
    try:
        lp = Polynomial(l)
        rp = Polynomial(r)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)
    p = lp - rp
    if (not all(is_finite(c) for c in p.coeffs)):
        print("Error: coefficient out of floating-point range")
        sys.exit(1)

    print(f"Reduced form: {p} = 0")
    if p.degree() > 0:
        print(f"Polynomial degree: {p.degree()}")
    if p.degree() > 2:
        print("The polynomial degree is strictly greater than 2, I can't solve.")
        sys.exit(1)

    match p.degree():
        case 0:
            if p.coeffs[0] == 0:
                print("Any real number is a solution.")
            else:
                print("No solution.")
        case 1:
            a = p.coeffs[1]
            b = p.coeffs[0]
            sol = -b / a
            if (not is_finite(sol)):
                print("Error: solution out of floating-point range")
                sys.exit(1)
            print(f"The solution is:\n{fmt(sol)}")
        case 2:
            a = p.coeffs[2]
            b = p.coeffs[1]
            c = p.coeffs[0]
            delta = b * b - 4 * a * c
            if (not is_finite(delta)):
                print("Error: discriminant out of floating-point range")
                sys.exit(1)
            if (is_near_zero(delta, b * b + abs(4 * a * c))):
                delta = 0
            if (delta > 0):
                sol1 = (-b - ft_sqrt(delta)) / (2 * a)
                sol2 = (-b + ft_sqrt(delta)) / (2 * a)
                print(f"Discriminant is strictly positive, the two solutions are:\n{fmt(sol1)}\n{fmt(sol2)}")
            elif (delta == 0):
                sol = -b / (2 * a)
                print(f"Discriminant is zero, the solution is:\n{fmt(sol)}")
            else:
                solr = -b / (2 * a)
                soli = abs(ft_sqrt(abs(delta)) / (2 * a))
                sol1 = f"{fmt(solr)} - {fmt(soli)}i"
                sol2 = f"{fmt(solr)} + {fmt(soli)}i"
                print(f"Discriminant is strictly negative, the two complex solutions are:\n{sol1}\n{sol2}")

if __name__ == "__main__":
    main()
