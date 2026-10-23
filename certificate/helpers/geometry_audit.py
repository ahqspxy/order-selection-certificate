"""Exact normalized-shape/scale mappings for the sixteen retained charts."""

from fractions import Fraction as F

from .partition_audit import partition


def validate_mapping(row):
    pivot, exponent = row["pivot"], row["exponent"]
    assert type(pivot) is int and pivot in range(4)
    assert type(exponent) is int and exponent in range(-14, -10)
    lower, upper = row["param_lo"], row["param_hi"]
    shape_lower, shape_upper = row["alo"], row["ahi"]
    assert len(lower) == len(upper) == len(shape_lower) == len(shape_upper) == 4
    assert shape_lower[pivot] == shape_upper[pivot] == 1
    assert all(0 <= F(a) < F(b) <= 1 for a, b in zip(lower, upper))
    for k, axis in enumerate(j for j in range(4) if j != pivot):
        assert (0 <= F(shape_lower[axis]) <= F(lower[k])
                < F(upper[k]) <= F(shape_upper[axis]) <= 1)
    scale_lower, scale_upper = map(F, row["scale"])
    unit = F(2) ** exponent
    mapped_lower = unit * (1 + F(lower[3]))
    mapped_upper = unit * (1 + F(upper[3]))
    assert 0 < scale_lower <= mapped_lower < mapped_upper <= scale_upper <= 1
    resolution = row["resolution"]
    assert len(resolution) == 2
    assert all(type(n) is int and n > 0 for n in resolution)
    assert resolution[0] % 5 == 0
    return scale_lower == mapped_lower and scale_upper == mapped_upper


def geometry(rows, pivot, exponent):
    assert rows
    exact_scales = 0
    for row in rows:
        assert (row["pivot"], row["exponent"]) == (pivot, exponent)
        exact_scales += int(validate_mapping(row))
    result = partition([(row["param_lo"], row["param_hi"])
                        for row in rows], 4)
    return {
        "pivot": pivot,
        "exponent": exponent,
        **result,
        "all_leaf_mappings_valid": True,
        "exact_scale_mappings": exact_scales,
        "scale_interval_exact": [str(F(2) ** exponent),
                                 str(F(2) ** (exponent + 1))],
    }
