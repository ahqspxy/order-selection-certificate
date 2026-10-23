"""Guarded interval calls shared by the complete replay. No stored outputs are read."""
from __future__ import annotations

import argparse
from collections import Counter
import ctypes as C
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
import json
import lzma
import math
import os
from pathlib import Path
import platform
import struct
import subprocess
import sys
import time
import traceback


GAMMA = F(74500001, 100000000)
UNIFORM = F(999999999886, 10**12)
SOURCES = {
    "normalized_dp": "normalized_dp.cpp",
    "refined_dp": "normalized_dp_refined.cpp",
    "fasttight": "fasttight4.cpp",
    "jointtight": "jointtight4.cpp",
    "jointtime": "jointtime4.cpp",
    "normalized_finite": "normalized_finite.cpp",
}
POTENTIAL = {"fasttight", "jointtight", "jointtime"}
PTR = C.POINTER(C.c_double)
FLAGS = ["-std=c++17", "-O3", "-fno-fast-math", "-ffp-contract=off", "-fPIC", "-shared"]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_hash(value):
    raw = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(raw.encode()).hexdigest()


def binary(value):
    return struct.pack(">d", value).hex()


def save(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def analytic_zero_bound():
    term = total = F(1)
    for k in range(1, 15):
        term *= F(32, 5) / k
        total += term
    require(total > 600, "Zero-mass exponential lower bound failed")
    value = GAMMA + GAMMA * GAMMA / (600 * (1 - GAMMA))
    require(value < F(749, 1000), "Zero-mass bound failed")
    return value


def build(cxx, source, target, include, grouped, out):
    command = [cxx, *FLAGS]
    if grouped:
        command.append("-DGROUP_COUNT=4")
    command += [str(source), "-I", str(include), "-o", str(target)]
    print("BUILD", target.name, flush=True)
    started = time.monotonic()
    result = subprocess.run(command, capture_output=True, text=True)
    (out / (target.stem + ".build.stdout.log")).write_text(result.stdout)
    (out / (target.stem + ".build.stderr.log")).write_text(result.stderr)
    record = {
        "command": command,
        "exit_code": result.returncode,
        "elapsed_seconds": time.monotonic() - started,
        "source_sha256": sha(source),
        "binary_sha256": sha(target) if result.returncode == 0 else None,
    }
    require(result.returncode == 0, "Build failed for " + target.name + ": " + result.stderr[-5000:])
    return record


def wrapper(method):
    gamma_expression = "source_tail::GAMMA" if method in POTENTIAL else "GAMMA"
    return f'''#include "{SOURCES[method]}"
extern "C" int replay_kernel(double *out) {{
  int q=0;
  for(int i=0;i<16;i++)for(int j=0;j<25;j++){{out[q++]=wi[i][j].lo;out[q++]=wi[i][j].hi;}}
  for(int i=0;i<16;i++)for(int j=0;j<26;j++){{out[q++]=ci[i][j].lo;out[q++]=ci[i][j].hi;}}
  return q;
}}
extern "C" int replay_gamma(double *out) {{
  I g=gamma_interval(.74500001);out[0]=g.lo;out[1]=g.hi;
  out[2]={gamma_expression}.lo;out[3]={gamma_expression}.hi;
  return strict_mode?1:0;
}}
'''


class Engine:
    def __init__(self, method, binary_path, kernel, environment):
        self.method = method
        self.environment = environment
        require(environment.certification_environment() == 0,
                "Unsupported floating-point environment before initialization")
        self.lib = lib = C.CDLL(str(binary_path))
        if method in POTENTIAL:
            lib.set_gain_beta.argtypes = [C.c_double]
            lib.set_gain_beta.restype = None
            lib.set_gain_beta(0.0)
            for name, path in (("init_kernel", kernel / "constants.txt"),
                               ("init_strict", kernel / "kernel_intervals.txt")):
                function = getattr(lib, name)
                function.argtypes = [C.c_char_p]
                function.restype = C.c_int
                require(function(os.fsencode(path)) == 0, name + " rejected the kernel")
            self.function = lib.evaluate
            self.function.argtypes = [PTR, PTR, C.c_double, C.c_int, C.c_int, C.c_int, PTR]
            self.diagnostic_count = 2
        else:
            lib.initialize.argtypes = [C.c_char_p]
            lib.initialize.restype = C.c_int
            require(lib.initialize(os.fsencode(kernel / "kernel_intervals.txt")) == 0,
                    "initialize rejected the kernel")
            name = {"normalized_dp": "evaluate_dp",
                    "refined_dp": "evaluate_refined", "normalized_finite": "evaluate_finite"}[method]
            self.function = getattr(lib, name)
            if method == "normalized_finite":
                self.function.argtypes = [PTR, PTR, C.c_int, C.c_double, C.c_double, PTR]
            else:
                self.function.argtypes = [PTR, PTR, C.c_int, C.c_int, C.c_double, C.c_double, PTR]
            self.diagnostic_count = 3
        self.function.restype = C.c_double
        lib.replay_kernel.argtypes = [PTR]
        lib.replay_kernel.restype = C.c_int
        lib.replay_gamma.argtypes = [PTR]
        lib.replay_gamma.restype = C.c_int
        actual = (C.c_double * 1632)()
        require(lib.replay_kernel(actual) == 1632, "Wrong kernel readback length")
        expected = list(map(float, (kernel / "kernel_intervals.txt").read_text().split()))
        require(len(expected) == 1632 and
                all(binary(a) == binary(b) for a, b in zip(actual, expected)),
                "Kernel readback mismatch")
        gamma = (C.c_double * 4)()
        require(lib.replay_gamma(gamma) == 1, "Strict interval mode is inactive")
        require(F(gamma[0]) <= GAMMA <= F(gamma[1]) and
                F(gamma[2]) <= GAMMA <= F(gamma[3]),
                "Observed Gamma intervals do not contain the exact target")
        require(environment.certification_environment() == 0,
                "Unsupported floating-point environment after initialization")
        self.identity = {
            "method": method,
            "binary_sha256": sha(binary_path),
            "kernel_interval_sha256": sha(kernel / "kernel_intervals.txt"),
            "strict_mode_observed": True,
            "kernel_readback_endpoints": 1632,
            "kernel_readback_bitwise_exact": True,
            "gamma_intervals": [[str(F(gamma[0])), str(F(gamma[1]))],
                                [str(F(gamma[2])), str(F(gamma[3]))]],
        }
        lib.audit_reset.argtypes = []
        lib.audit_reset.restype = None
        lib.audit_get_counts.argtypes = [C.POINTER(C.c_ulonglong)]
        lib.audit_get_counts.restype = None
        lib.audit_get_ranges.argtypes = [PTR]
        lib.audit_get_ranges.restype = None
        lib.audit_reset()

    def call(self, lo, hi, resolution, scale_lower=1.0, scale_upper=1.0):
        require(len(lo) == len(hi) == 4 and
                all(math.isfinite(a) and math.isfinite(b) and 0 <= a <= b <= 1
                    for a, b in zip(lo, hi)), "Invalid parameter box")
        require(math.isfinite(scale_lower) and math.isfinite(scale_upper) and
                0 <= scale_lower <= scale_upper <= 1 and scale_upper > 0,
                "Invalid scale interval")
        if self.method == "normalized_finite":
            require(type(resolution) is int and resolution > 0 and resolution % 5 == 0,
                    "Invalid scalar resolution")
        else:
            require(isinstance(resolution, list) and len(resolution) == 2 and
                    all(type(n) is int and n > 0 for n in resolution) and
                    resolution[0] % 5 == 0,
                    "Invalid paired resolution: positive integers and time resolution divisible by 5 required")
        require(self.environment.certification_environment() == 0,
                "Environment check failed before numerical call")
        a, b = (C.c_double * 4)(*lo), (C.c_double * 4)(*hi)
        diagnostic = (C.c_double * self.diagnostic_count)()
        started = time.monotonic()
        if self.method in POTENTIAL:
            value = self.function(a, b, float(GAMMA), resolution[0], resolution[1], 10, diagnostic)
        elif self.method == "normalized_finite":
            value = self.function(a, b, resolution, scale_lower, scale_upper, diagnostic)
        else:
            value = self.function(a, b, resolution[0], resolution[1], scale_lower, scale_upper, diagnostic)
        elapsed = time.monotonic() - started
        require(self.environment.certification_environment() == 0,
                "Environment check failed after numerical call")
        require(math.isfinite(value) and all(math.isfinite(x) for x in diagnostic),
                "Nonfinite result or diagnostic")
        require(F(value) <= UNIFORM, "Fresh value exceeds the exact uniform bound: " + str(F(value)))
        return {
            "method": self.method, "upper": value, "upper_exact": str(F(value)),
            "upper_hex": value.hex(), "diagnostic": list(diagnostic),
            "diagnostic_hex": [x.hex() for x in diagnostic], "seconds": elapsed,
            "input_lo": lo, "input_hi": hi, "resolution": resolution,
            "scale_lower_exact": str(F(scale_lower)), "scale_upper_exact": str(F(scale_upper)),
            "environment_before_after_passed": True,
        }

    def guard_statistics(self):
        counts, ranges = (C.c_ulonglong * 5)(), (C.c_double * 4)()
        self.lib.audit_get_counts(counts)
        self.lib.audit_get_ranges(ranges)
        return {"counts": list(counts),
                "ranges": [x if math.isfinite(x) else None for x in ranges],
                "count_order": ["intervals", "divisions", "exponentials", "square_roots", "max_sqrt_steps"],
                "range_order": ["minimum_divisor", "maximum_exp_argument", "minimum_sqrt_argument", "maximum_sqrt_argument"]}

