"""Exact rational verification of the fixed kernel and loaded coefficient intervals."""
from fractions import Fraction as F
import json, math
from interval_engine import require

def verify_kernel(root):
    """Certify the compact kernel expansion and every coefficient enclosure."""
    directory = root / "source/guarded/robust_745/bulk_r1_exact"
    cycles = json.loads((directory / "cycles.json").read_text())
    compact = json.loads((directory / "kernel_compact.json").read_text())
    step = F(2, 5)
    theta = F(999999999, 1000000000)
    groups = [0, 1, 5, 8, 16]
    require(cycles["source_groups"] == compact["source_groups"] == groups,
            "Kernel source groups mismatch")
    require(F(cycles["source_step"]) == F(cycles["output_step"]) == step,
            "Kernel mesh mismatch")
    require(cycles["source_cells"] == 16 and cycles["output_cells"] == 25
            and cycles["hub_row"] == 0 and cycles["hub_col"] == 24,
            "Kernel dimensions or hub mismatch")
    for description in (cycles, compact):
        require(F(description["theta_numerator"], description["theta_denominator"]) == theta,
                "Kernel perturbation factor mismatch")
    require(compact["denominator"] == 10**21 and len(compact["table"]) == 23,
            "Compact kernel dimensions mismatch")
    expanded = {}
    for column in compact["table"]:
        j, numerators = column["column"], column["group_numerators"]
        require(type(j) is int and j in range(1, 24) and len(numerators) == 3,
                "Invalid compact kernel column")
        require(all(type(n) is int for n in numerators), "Noninteger compact coefficient")
        for group in range(1, 4):
            for i in range(groups[group], groups[group + 1]):
                require((i, j) not in expanded, "Duplicate compact kernel coefficient")
                expanded[i, j] = F(numerators[group - 1], compact["denominator"])
    require(len(expanded) == 345, "Incomplete compact kernel expansion")
    observed = {}
    for cycle in cycles["cycles"]:
        key = cycle["i"], cycle["j"]
        require(key not in observed, "Duplicate kernel cycle")
        observed[key] = F(cycle["numerator"], cycle["denominator"])
    require(observed == expanded, "The 69 compact integers do not match all 345 cycles")
    delta = [[F(0) for _ in range(25)] for _ in range(16)]
    for (i, j), coefficient in expanded.items():
        amplitude = theta * coefficient
        for ii, jj, sign in [(i, j, 1), (0, 24, 1), (i, 24, -1), (0, j, -1)]:
            delta[ii][jj] += sign * amplitude
    require(all(sum(row) == 0 for row in delta), "Kernel row marginal mismatch")
    require(all(sum(delta[i][j] for i in range(16)) == 0 for j in range(25)),
            "Kernel column marginal mismatch")

    def exp_interval(x, degree=160):
        require(0 <= x <= 10, "Taylor argument outside certified range")
        total = term = F(1)
        for k in range(1, degree + 1):
            term *= x / k
            total += term
        remainder = (term * x / (degree + 1)) / (1 - x / (degree + 2))
        return total, total + remainder

    def add(a, b):
        return a[0] + b[0], a[1] + b[1]

    def negate(a):
        return -a[1], -a[0]

    def subtract(a, b):
        return add(a, negate(b))

    def multiply(a, b):
        products = [x * y for x in a for y in b]
        return min(products), max(products)

    def inverse(a):
        require(a[0] > 0, "Kernel rational denominator not positive")
        return F(1) / a[1], F(1) / a[0]

    def point(x):
        return x, x

    positive_exp = [exp_interval(j * step) for j in range(26)]
    negative_exp = [inverse(item) for item in positive_exp]
    omega = [subtract(negative_exp[j], negative_exp[j + 1]) for j in range(25)]
    require(all(item[0] > 0 for item in omega), "Output-cell mass not positive")
    expected, positive_lower_bounds, zero_rectangles = [], [], []
    for i in range(16):
        for j in range(25):
            coefficient = delta[i][j]
            mass = (add(multiply(point(step), omega[j]), point(coefficient))
                    if i < j else point(coefficient))
            require(mass[0] >= 0, f"Kernel positivity fails at ({i},{j})")
            if mass[0] > 0:
                require(mass[0] > F(698, 10**14), "Kernel positive-mass margin fails")
                positive_lower_bounds.append(mass[0])
            else:
                zero_rectangles.append((i, j))
            expected.append(multiply(point(coefficient / step), inverse(omega[j])))
    for i in range(16):
        total = F(0)
        for j in range(26):
            expected.append(multiply(point(total / step), positive_exp[j]))
            if j < 25:
                total += delta[i][j]
    pairs = []
    for line in (directory / "kernel_intervals.txt").read_text().splitlines():
        tokens = line.split()
        require(len(tokens) == 2, "Malformed kernel interval")
        values = [float(token) for token in tokens]
        require(all(math.isfinite(value) for value in values), "Nonfinite kernel interval")
        pairs.append(tuple(map(F, values)))
    require(len(expected) == len(pairs) == 816, "Unexpected kernel interval count")
    for index, ((lower, upper), (a, b)) in enumerate(zip(pairs, expected)):
        require(lower <= a <= b <= upper, "Kernel interval enclosure fails: " + str(index))
    triangle = subtract(point(F(1)), multiply(point(1 + step), negative_exp[1]))
    require(triangle[0] > F(6155, 100000), "Positive diagonal-cell lower bound fails")
    print("PASS_EXACT_KERNEL", "69_INTEGERS", "345_CYCLES", "816_ENCLOSURES", flush=True)
    return {
        "status": "PASS_EXACT_RATIONAL_KERNEL",
        "compact_integers_checked": 69,
        "cycles_checked": 345,
        "coefficient_intervals_certified": 816,
        "row_column_marginals_exact": True,
        "pointwise_positivity_certified": True,
        "strict_positive_rectangle_lower_bounds": len(positive_lower_bounds),
        "zero_rectangle_lower_bounds": len(zero_rectangles),
        "smallest_positive_lower_exact": str(min(positive_lower_bounds)),
        "positive_diagonal_00_lower_gt": "1231/20000",
        "method": "Exact rational degree-160 positive Taylor sums with an explicit geometric tail; no numerical transcendental library.",
    }

