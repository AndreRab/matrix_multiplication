"""
Implementation of recursive matrix multiplication algorithms: classical (Binét) and
Strassen.  The module exposes a `MatrixMultiplier` class that can multiply two
square matrices of size ``n × n`` where ``n`` is a power of two.  The
``MatrixMultiplier`` chooses between Strassen's algorithm and the classical
block multiplication (here referred to as Binét) depending on a user-supplied
threshold size.  When the current matrix size is less than or equal to this
threshold, Strassen's algorithm is applied; otherwise, the classical block
multiplication is used.  Both algorithms are implemented recursively.

In addition to computing the product matrix, the multiplier counts the number
of floating‑point operations (additions/subtractions and multiplications)
performed during the computation.  This allows empirical comparison of the
efficiency of different thresholds.  A convenience function ``run_experiment``
is provided to generate random matrices, run the multiplication for a range
of matrix sizes and thresholds, and collect timing and operation count data.

The code is written using plain Python lists to explicitly account for each
arithmetic operation.  It is not optimized for speed but is suitable for
educational purposes, small problem sizes and to illustrate differences
between algorithms.
"""

from __future__ import annotations

import random
import time
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple


Matrix = List[List[float]]


def _is_power_of_two(n: int) -> bool:
    """Return True when ``n`` is a positive power of two."""
    return n > 0 and (n & (n - 1)) == 0


def _validate_matrix(matrix: Matrix, name: str) -> None:
    """Validate that ``matrix`` is a non-empty square matrix."""
    if not matrix or not matrix[0]:
        raise ValueError(f"{name} must be a non-empty square matrix")

    n = len(matrix)
    if any(len(row) != n for row in matrix):
        raise ValueError(f"{name} must be square (got {n} rows with inconsistent column counts)")

    if not _is_power_of_two(n):
        raise ValueError(f"{name} must have size equal to a power of two (got {n}x{n})")


@dataclass
class MatrixMultiplier:
    """Performs recursive matrix multiplication with operation counting.

    Parameters
    ----------
    threshold : int
        The maximum matrix dimension at which Strassen's algorithm will be
        applied.  When the current subproblem dimension ``n`` is greater than
        ``threshold``, the classical block multiplication (Binét) is used
        instead.  Note that the threshold is an absolute dimension (e.g., 8
        corresponds to matrices up to 8×8).  For matrix sizes ``n <= 1`` the
        direct multiplication is performed without further recursion.
    """
    threshold: int

    def __post_init__(self) -> None:
        self.operations: int = 0

    def reset_operations(self) -> None:
        """Reset the operation counter to zero."""
        self.operations = 0

    # ------------------------------------------------------------------
    # Helpers for matrix operations with operation counting
    def _add(self, A: Matrix, B: Matrix) -> Matrix:
        """Element‑wise addition of two matrices with operation counting."""
        n = len(A)
        C = [[0.0] * n for _ in range(n)]
        for i in range(n):
            rowA = A[i]
            rowB = B[i]
            rowC = C[i]
            for j in range(n):
                # one floating point addition
                rowC[j] = rowA[j] + rowB[j]
                self.operations += 1
        return C

    def _sub(self, A: Matrix, B: Matrix) -> Matrix:
        """Element‑wise subtraction of two matrices with operation counting."""
        n = len(A)
        C = [[0.0] * n for _ in range(n)]
        for i in range(n):
            rowA = A[i]
            rowB = B[i]
            rowC = C[i]
            for j in range(n):
                # one floating point subtraction is counted as an addition
                rowC[j] = rowA[j] - rowB[j]
                self.operations += 1
        return C

    def _combine_quadrants(self, C11: Matrix, C12: Matrix, C21: Matrix, C22: Matrix) -> Matrix:
        """Combine four submatrices into a single matrix."""
        n2 = len(C11)
        n = n2 * 2
        C = [[0.0] * n for _ in range(n)]
        for i in range(n2):
            # copy top rows
            C[i][:n2] = C11[i]
            C[i][n2:] = C12[i]
        for i in range(n2):
            # copy bottom rows
            C[i + n2][:n2] = C21[i]
            C[i + n2][n2:] = C22[i]
        return C

    def _split_quadrants(self, A: Matrix) -> Tuple[Matrix, Matrix, Matrix, Matrix]:
        """Split a square matrix into four equal quadrants."""
        n2 = len(A) // 2
        return (
            [row[:n2] for row in A[:n2]],
            [row[n2:] for row in A[:n2]],
            [row[:n2] for row in A[n2:]],
            [row[n2:] for row in A[n2:]],
        )

    # ------------------------------------------------------------------
    # Core multiplication routine
    def multiply(self, A: Matrix, B: Matrix) -> Matrix:
        """Validate inputs and multiply matrices ``A`` and ``B``."""
        _validate_matrix(A, "A")
        _validate_matrix(B, "B")
        if len(A) != len(B):
            raise ValueError(
                f"A and B must have the same dimensions, got {len(A)}x{len(A)} and {len(B)}x{len(B)}"
            )
        return self._multiply_recursive(A, B)

    def _multiply_recursive(self, A: Matrix, B: Matrix) -> Matrix:
        """Multiply matrices ``A`` and ``B`` using the selected algorithm.

        The algorithm used depends on the current dimension and the ``threshold``
        attribute.  This method updates the ``operations`` counter.
        """
        n = len(A)
        # Base case: 1×1 matrix multiplication
        if n == 1:
            # one multiplication
            self.operations += 1
            return [[A[0][0] * B[0][0]]]

        # Decide algorithm based on threshold
        if n <= self.threshold:
            return self._strassen(A, B)
        else:
            return self._binet(A, B)

    # ------------------------------------------------------------------
    # Classical recursive multiplication (Binét)
    def _binet(self, A: Matrix, B: Matrix) -> Matrix:
        """Classical block multiplication (referred to as Binét)."""
        A11, A12, A21, A22 = self._split_quadrants(A)
        B11, B12, B21, B22 = self._split_quadrants(B)

        # Compute subproducts recursively (8 multiplications)
        M1 = self._multiply_recursive(A11, B11)
        M2 = self._multiply_recursive(A12, B21)
        M3 = self._multiply_recursive(A11, B12)
        M4 = self._multiply_recursive(A12, B22)
        M5 = self._multiply_recursive(A21, B11)
        M6 = self._multiply_recursive(A22, B21)
        M7 = self._multiply_recursive(A21, B12)
        M8 = self._multiply_recursive(A22, B22)

        # Combine results: C11 = M1 + M2; C12 = M3 + M4; C21 = M5 + M6; C22 = M7 + M8
        C11 = self._add(M1, M2)
        C12 = self._add(M3, M4)
        C21 = self._add(M5, M6)
        C22 = self._add(M7, M8)
        return self._combine_quadrants(C11, C12, C21, C22)

    # ------------------------------------------------------------------
    # Strassen's recursive multiplication
    def _strassen(self, A: Matrix, B: Matrix) -> Matrix:
        """Strassen's matrix multiplication algorithm."""
        n = len(A)
        if n == 1:
            return [[A[0][0] * B[0][0]]]

        A11, A12, A21, A22 = self._split_quadrants(A)
        B11, B12, B21, B22 = self._split_quadrants(B)

        # Compute the seven products recursively
        # P1 = (A11 + A22)*(B11 + B22)
        S1 = self._add(A11, A22)
        S2 = self._add(B11, B22)
        P1 = self._multiply_recursive(S1, S2)
        # P2 = (A21 + A22)*B11
        S3 = self._add(A21, A22)
        P2 = self._multiply_recursive(S3, B11)
        # P3 = A11*(B12 - B22)
        S4 = self._sub(B12, B22)
        P3 = self._multiply_recursive(A11, S4)
        # P4 = A22*(B21 - B11)
        S5 = self._sub(B21, B11)
        P4 = self._multiply_recursive(A22, S5)
        # P5 = (A11 + A12)*B22
        S6 = self._add(A11, A12)
        P5 = self._multiply_recursive(S6, B22)
        # P6 = (A21 - A11)*(B11 + B12)
        S7 = self._sub(A21, A11)
        S8 = self._add(B11, B12)
        P6 = self._multiply_recursive(S7, S8)
        # P7 = (A12 - A22)*(B21 + B22)
        S9 = self._sub(A12, A22)
        S10 = self._add(B21, B22)
        P7 = self._multiply_recursive(S9, S10)

        # Compute the four quadrants of C
        # C11 = P1 + P4 - P5 + P7
        C11 = self._add(self._sub(self._add(P1, P4), P5), P7)
        # C12 = P3 + P5
        C12 = self._add(P3, P5)
        # C21 = P2 + P4
        C21 = self._add(P2, P4)
        # C22 = P1 - P2 + P3 + P6
        C22 = self._add(self._add(self._sub(P1, P2), P3), P6)

        return self._combine_quadrants(C11, C12, C21, C22)


def _generate_random_matrix(n: int) -> Matrix:
    """Generate an n×n matrix with random floating‑point entries."""
    return [[random.random() for _ in range(n)] for _ in range(n)]


def classical_reference(A: Matrix, B: Matrix) -> Matrix:
    """Reference triple-loop multiplication used for correctness checks."""
    _validate_matrix(A, "A")
    _validate_matrix(B, "B")
    if len(A) != len(B):
        raise ValueError("A and B must have the same dimensions")

    n = len(A)
    C = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for k in range(n):
            for j in range(n):
                C[i][j] += A[i][k] * B[k][j]
    return C


def matrices_close(A: Matrix, B: Matrix, tol: float = 1e-9) -> bool:
    """Return True when matrices are equal up to ``tol``."""
    if len(A) != len(B):
        return False
    return all(
        abs(A[i][j] - B[i][j]) <= tol
        for i in range(len(A))
        for j in range(len(A))
    )


def run_experiment(
    k_values: List[int],
    l_values: List[int],
    *,
    seed: Optional[int] = None,
) -> Dict[int, Dict[int, Tuple[float, int]]]:
    """Run timing and operation count experiments.

    Parameters
    ----------
    k_values : list of ints
        Exponents defining matrix sizes ``2^k`` to test (e.g. k=2 tests
        4×4 matrices).
    l_values : list of ints
        Threshold exponents ``l``.  For a given ``l``, the threshold
        dimension used in the multiplier is ``2**l``.  For each k and l
        combination the experiment generates two random matrices of size
        ``2**k``, multiplies them using a new ``MatrixMultiplier`` with
        threshold ``2**l``, and records execution time and number of
        floating‑point operations.

    Returns
    -------
    results : dict
        A nested dictionary ``results[l][k] = (elapsed_time, operations)``
        with the measured execution time in seconds and the counted number
        of floating‑point operations for each combination of threshold and
        matrix size.
    """
    if seed is not None:
        random.seed(seed)

    results: Dict[int, Dict[int, Tuple[float, int]]] = {}
    for l in sorted(l_values):
        threshold = 2 ** l
        results[l] = {}
        for k in sorted(k_values):
            n = 2 ** k
            # Generate random matrices
            A = _generate_random_matrix(n)
            B = _generate_random_matrix(n)
            # Create a multiplier with the current threshold
            multiplier = MatrixMultiplier(threshold)
            start = time.perf_counter()
            multiplier.multiply(A, B)
            end = time.perf_counter()
            elapsed = end - start
            results[l][k] = (elapsed, multiplier.operations)

    return results


def self_test() -> None:
    """Run a small deterministic correctness check."""
    test_cases = [
        (
            [[1.0]],
            [[7.0]],
        ),
        (
            [[1.0, 2.0], [3.0, 4.0]],
            [[5.0, 6.0], [7.0, 8.0]],
        ),
        (
            [
                [1.0, 2.0, 3.0, 4.0],
                [5.0, 6.0, 7.0, 8.0],
                [9.0, 10.0, 11.0, 12.0],
                [13.0, 14.0, 15.0, 16.0],
            ],
            [
                [16.0, 15.0, 14.0, 13.0],
                [12.0, 11.0, 10.0, 9.0],
                [8.0, 7.0, 6.0, 5.0],
                [4.0, 3.0, 2.0, 1.0],
            ],
        ),
    ]

    for threshold in (1, 2, 4):
        multiplier = MatrixMultiplier(threshold)
        for A, B in test_cases:
            product = multiplier.multiply(A, B)
            reference = classical_reference(A, B)
            if not matrices_close(product, reference):
                raise AssertionError(
                    f"Incorrect result for threshold={threshold}: {product!r} != {reference!r}"
                )
            multiplier.reset_operations()


if __name__ == "__main__":
    self_test()

    # Example usage for quick manual testing. Runs a small experiment
    # after a correctness check and prints the collected measurements.
    ks = [2, 3]
    ls = [1, 2]
    res = run_experiment(ks, ls, seed=42)
    for l, rdict in res.items():
        print(f"Threshold 2^{l} (n={2**l})")
        for k, (t, ops) in rdict.items():
            print(f"  n=2^{k} -> time {t:.6f}s, ops {ops}")
