"""Exact dyadic partitions of a closed unit cube.

Coordinates are converted to exact rationals and embedded in an integer grid.
Each accepted terminal node equals its one input box. Every internal node
splits its parent into complementary children and assigns each input box to
exactly one child. Induction therefore proves exact coverage and disjoint
interiors. The preliminary volume check is not used alone to infer coverage.
Neighboring closed boxes may share their faces.
"""

from fractions import Fraction as F
import math


def partition(data, dimension):
    assert __debug__, "Run without Python -O"
    assert type(dimension) is int and dimension > 0
    denominators = [1] * dimension
    coordinates = []
    for lower, upper in data:
        lower, upper = tuple(map(F, lower)), tuple(map(F, upper))
        assert len(lower) == len(upper) == dimension
        for axis in range(dimension):
            assert 0 <= lower[axis] < upper[axis] <= 1
            for coordinate in (lower[axis], upper[axis]):
                denominator = coordinate.denominator
                assert denominator & (denominator - 1) == 0
                denominators[axis] = max(denominators[axis], denominator)
        coordinates.append((lower, upper))
    assert coordinates, "Empty cover"
    boxes = [
        (tuple(int(x * denominators[j]) for j, x in enumerate(lower)),
         tuple(int(x * denominators[j]) for j, x in enumerate(upper)))
        for lower, upper in coordinates
    ]
    del coordinates
    assert len(set(boxes)) == len(boxes), "Duplicate boxes"
    assert sum(math.prod(b - a for a, b in zip(lower, upper))
               for lower, upper in boxes) == math.prod(denominators)
    pending = [(list(range(len(boxes))), (0,) * dimension,
                tuple(denominators))]
    nodes = 0
    while pending:
        indices, lower, upper = pending.pop()
        nodes += 1
        assert indices
        if len(indices) == 1:
            assert boxes[indices[0]] == (lower, upper)
            continue
        split = False
        for axis in range(dimension):
            if (lower[axis] + upper[axis]) % 2:
                continue
            midpoint = (lower[axis] + upper[axis]) // 2
            if midpoint in (lower[axis], upper[axis]):
                continue
            left, right = [], []
            for index in indices:
                box_lower, box_upper = boxes[index]
                if box_upper[axis] <= midpoint:
                    left.append(index)
                elif box_lower[axis] >= midpoint:
                    right.append(index)
                else:
                    break
            else:
                if left and right:
                    left_upper, right_lower = list(upper), list(lower)
                    left_upper[axis] = right_lower[axis] = midpoint
                    pending.extend([(left, lower, tuple(left_upper)),
                                    (right, tuple(right_lower), upper)])
                    split = True
                    break
        assert split, "No complete exact dyadic bisection"
    assert nodes == 2 * len(boxes) - 1
    return {
        "dimension": dimension,
        "leaves": len(boxes),
        "volume_exact": "1",
        "coverage_method": "exact_recursive_bisection",
        "binary_tree_nodes": nodes,
        "coordinate_denominators": denominators,
        "closed_box_coverage": True,
    }
