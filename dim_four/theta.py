from typing import List, Self


def hadamard(coordinates):
    """
    Compute the Hadamard transformation using a recursive strategy
    """
    n = len(coordinates)

    # Base case, compute the sum and difference
    if n == 2:
        a, b = coordinates
        return (a + b, a - b)

    # Otherwise compute hadamard over each half
    m = n // 2
    a = hadamard(coordinates[:m])
    b = hadamard(coordinates[m:])

    return tuple(a[i] + b[i] for i in range(m)) + tuple(a[i] - b[i] for i in range(m))


class ThetaPoint:
    """
    A theta point in the level-2 theta model with 16 projective coordinates.
    """

    def __init__(self, coords: List):
        assert len(coords) == 16, "ThetaPoint needs exactly 16 coordinates"
        self.coords = list(coords)

    def __getitem__(self, i: int):
        return self.coords[i]

    def __setitem__(self, i: int, value):
        self.coords[i] = value

    def __repr__(self):
        return f"ThetaPoint({self.coords})"

    def square(self) -> Self:
        return ThetaPoint([c**2 for c in self.coords])

    def hadamard(self) -> Self:
        return ThetaPoint(hadamard(self.coords))

    def pointwise_mul(self, other: Self) -> Self:
        return ThetaPoint([a * b for a, b in zip(self.coords, other.coords)])

    def compress(self):
        """
        Keep only the 10 independent coordinates: [x0,x1,x2,x3,x5,x6,x7,x9,x11,x15].
        """
        idx = [0, 1, 2, 3, 5, 6, 7, 9, 11, 15]
        return ThetaPointCompressed([self.coords[i] for i in idx])


class ThetaPointCompressed:
    """
    A compressed theta point using only the 10 independent coordinates
    of a dimension-4 theta null point.
    """

    def __init__(self, coords: List):
        assert len(coords) == 10, "ThetaPointCompressed needs exactly 10 coordinates"
        self.coords = list(coords)

    def __getitem__(self, i: int):
        return self.coords[i]

    def __setitem__(self, i: int, value):
        self.coords[i] = value

    def __repr__(self):
        return f"ThetaPointCompressed({self.coords})"

    def decompress(self) -> ThetaPoint:
        """
        Expand back to the full 16-coordinate ThetaPoint.
        """
        z0, z1, z2, z3, z4, z5, z6, z7, z8, z9 = self.coords
        return ThetaPoint(
            [z0, z1, z2, z3, z2, z4, z5, z6, z1, z7, z4, z8, z3, z8, z6, z9]
        )

    def square(self) -> Self:
        return ThetaPointCompressed([c**2 for c in self.coords])

    def proj_batch_pseudo_inversion(self) -> Self:
        """
        Projective pseudo-inverse: coordinate i becomes (product of all
        coordinates except i), which is proportional to 1/x_i without
        ever performing a field inversion.
        """
        x = self.coords
        n = len(x)

        prefix = [x[0]]
        for i in range(1, n):
            prefix.append(prefix[-1] * x[i])

        suffix = [1]
        for i in range(1, n):
            suffix.append(suffix[-1] * x[n - i])

        result = [suffix[-1]]
        for i in range(1, n):
            result.append(suffix[n - 1 - i] * prefix[i - 1])

        return ThetaPointCompressed(result)

    def hadamard(self) -> Self:
        """
        Hadamard transform exploiting the symmetry of the compressed
        representation.
        """
        z0, z1, z2, z3, z4, z5, z6, z7, z8, z9 = self.coords

        t0, t1 = z0 + z9, z0 - z9
        t2, t3 = z5 + z7, z5 - z7
        t4, t5 = z2 + z8, z2 - z8
        t6, t7 = z1 + z6, z1 - z6
        t8, t9 = z3 + z4, z3 - z4

        s0, s1 = t0 + t2, t0 - t2
        s2, s3 = t1 + t3, t1 - t3
        s4, s5 = t6 + t4, t6 - t4

        r0 = t5 + t5
        r1 = t7 + t7
        r2 = t9 + t9
        r3 = t8 + t8
        r4 = s4 + s4
        r5 = s5 + s5

        u0, u1 = s0 + r3, s0 - r3

        return ThetaPointCompressed(
            [
                u0 + r4,
                s2 + r0,
                s3 + r1,
                s1 + r2,
                s1 - r2,
                u1 + r5,
                s2 - r0,
                u1 - r5,
                s3 - r1,
                u0 - r4,
            ]
        )


class ThetaStructure:
    def __init__(
        self,
        null_point: ThetaPoint,
        inv_null_point: ThetaPoint,
        inv_null_point_dual_sq: ThetaPoint,
    ):
        self._null_point = null_point
        self._inv_null_point = inv_null_point
        self._inv_null_point_dual_sq = inv_null_point_dual_sq

    @classmethod
    def from_null_point(cls, null_point: ThetaPoint) -> Self:
        null_point_comp = null_point.compress()

        inv_null_point_comp = null_point_comp.proj_batch_pseudo_inversion()
        inv_null_point_dual_sq_comp = (
            null_point_comp.square().hadamard().proj_batch_pseudo_inversion()
        )

        return cls(
            null_point,
            inv_null_point_comp.decompress(),
            inv_null_point_dual_sq_comp.decompress(),
        )

    @classmethod
    def from_inv_dual(cls, inv_null_point_dual: ThetaPoint) -> Self:
        inv_null_point_dual_comp = inv_null_point_dual.compress()

        null_point_dual_comp = inv_null_point_dual_comp.proj_batch_pseudo_inversion()
        null_point_comp = null_point_dual_comp.hadamard()

        inv_null_point_dual_sq_comp = (
            null_point_comp.square().hadamard().proj_batch_pseudo_inversion()
        )
        inv_null_point_comp = null_point_comp.proj_batch_pseudo_inversion()

        return cls(
            null_point_comp.decompress(),
            inv_null_point_comp.decompress(),
            inv_null_point_dual_sq_comp.decompress(),
        )

    def null_point(self) -> ThetaPoint:
        return self._null_point

    def double(self, point: ThetaPoint) -> ThetaPoint:
        P = point.square()
        P = P.hadamard()
        P = P.square()
        P = P.pointwise_mul(self._inv_null_point_dual_sq)
        P = P.hadamard()
        P = P.pointwise_mul(self._inv_null_point)
        return P

    def double_iter(self, point: ThetaPoint, n: int) -> ThetaPoint:
        for _ in range(n):
            point = self.double(point)
        return point

    def special_moduli(self):
        """
        1I + 4S + 9M + 29a
        """
        x0, _, _, x3, _, x5, x6, _, _, x9, _, _, _, _, _, x15 = self._null_point.coords

        # u0, u1, u2, u3 = hadamard(x0, x6, x9, x15)
        (t0, t1) = (x0 + x6, x0 - x6)
        (t2, t3) = (x9 + x15, x9 - x15)
        (u0, u2) = (t0 + t2, t0 - t2)
        (u1, u3) = (t1 + t3, t1 - t3)

        s0 = u0**2
        s1 = u1**2
        s2 = u2**2
        s3 = u3**2

        t = 16 * (x3 * x5)

        # u0, u1, u2, u3 = hadamard(z0, z1, z2, z3)
        (t0, t1) = (s0 + s1, s0 - s1)
        (t2, t3) = (s2 + s3, s2 - s3)
        u2 = t0 - t2
        (u1, u3) = (t1 + t3, t1 - t3)

        (s4, s5) = (u1 + t, u1 - t)
        (s6, s7) = (u2 + t, u2 - t)
        (s8, s9) = (u3 + t, u3 - t)

        # 1I + 9M
        theta = [s0, s1, s2, s3, s4, s5, s6, s7, s8, s9]
        theta_n = [theta[i] / theta[0] for i in range(1, len(theta))]

        # 3a
        theta_n[0] = 2 * theta_n[0]
        theta_n[1] = 2 * theta_n[1]
        theta_n[2] = 2 * theta_n[2]

        return theta_n

    def absolute_invariants(self):
        """
        1I + 4S + 9M + 29a + 18 S + 18 M + 24 a
        """
        th = self.special_moduli()

        s1, s2, s3, s4, s5, s6, s7, s8, s9 = th
        t1, t2, t3, t4, t5, t6, t7, t8, t9 = (
            s1**2,
            s2**2,
            s3**2,
            s4**2,
            s5**2,
            s6**2,
            s7**2,
            s8**2,
            s9**2,
        )
        r1, r2, r3, r4, r5, r6, r7, r8, r9 = (
            s1 * t1,
            s2 * t2,
            s3 * t3,
            s4 * t4,
            s5 * t5,
            s6 * t6,
            s7 * t7,
            s8 * t8,
            s9 * t9,
        )
        u1, u2, u3, u4, u5, u6, u7, u8, u9 = (
            t1**2,
            t2**2,
            t3**2,
            t4**2,
            t5**2,
            t6**2,
            t7**2,
            t8**2,
            t9**2,
        )
        v1, v2, v3, v4, v5, v6, v7, v8, v9 = (
            s1 * u1,
            s2 * u2,
            s3 * u3,
            s4 * u4,
            s5 * u5,
            s6 * u6,
            s7 * u7,
            s8 * u8,
            s9 * u9,
        )

        J2 = t1 + t2 + t3 + t4 + t5 + t6 + t7 + t8 + t9
        J3 = r1 + r2 + r3 + r4 + r5 + r6 + r7 + r8 + r9
        J4 = u1 + u2 + u3 + u4 + u5 + u6 + u7 + u8 + u9
        J5 = v1 + v2 + v3 + v4 + v5 + v6 + v7 + v8 + v9

        absolute_invariants = (J2, J3, J4, J5)

        return absolute_invariants
