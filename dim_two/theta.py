class ThetaPointDimTwo:
    """
    A theta point in the level-2 theta model with 4 projective coordinates.
    """

    def __init__(self, coords):
        assert len(coords) == 4, "ThetaPointDimTwo needs exactly 4 coordinates"
        self.coords = list(coords)

    def __getitem__(self, i):
        return self.coords[i]

    def __setitem__(self, i, value):
        self.coords[i] = value

    def __iter__(self):
        return iter(self.coords)

    def __repr__(self):
        return f"ThetaPointDimTwo({self.coords})"

    def square(self):
        return ThetaPointDimTwo([c * c for c in self.coords])

    def hadamard(self):
        x, y, z, t = self.coords
        t1, t2 = x + y, x - y
        t3, t4 = z + t, z - t
        return ThetaPointDimTwo([t1 + t3, t2 + t4, t1 - t3, t2 - t4])

    def scholten_change_of_basis(self):
        x, y, z, t = self.coords
        t1, t2 = x + y, x - y
        t3, t4 = z + t, z - t
        return ThetaPointDimTwo([t1 + t4, t2 + t3, t1 - t4, t2 - t3])

    def pointwise_mul(self, other):
        return ThetaPointDimTwo([a * b for a, b in zip(self.coords, other.coords)])

    def pointwise_sub(self, other):
        return ThetaPointDimTwo([a - b for a, b in zip(self.coords, other.coords)])


class ThetaStructureDimTwo:
    def __init__(self, null_point: ThetaPointDimTwo):
        self._null_point = null_point
        self.arithmetic_precomputation(null_point)

    def arithmetic_precomputation(self, null_point):
        """
        Projective ratios [X0, Y0, Z0, T0] used by `diff_addition`
        """
        AA, BB, CC, DD = null_point.square().hadamard().coords

        t1 = AA * BB
        t2 = CC * DD

        T0 = t1 * CC
        Z0 = t1 * DD
        Y0 = t2 * AA
        X0 = t2 * BB

        self._precomputation = [X0, Y0, Z0, T0]

    def diff_addition(self, P, Q, PQ):
        """
        Given P, Q and PQ = P - Q, return P + Q.
        """
        X0, Y0, Z0, T0 = self._precomputation

        p1, p2, p3, p4 = P.square().hadamard().coords
        q1, q2, q3, q4 = Q.square().hadamard().coords

        xp = X0 * p1 * q1
        yp = Y0 * (p2 * q2)
        zp = Z0 * (p3 * q3)
        tp = T0 * (p4 * q4)

        X, Y, Z, T = ThetaPointDimTwo([xp, yp, zp, tp]).hadamard().coords

        # Four divisions by the coordinates of PQ, cleared with 10 multiplications
        PQx, PQy, PQz, PQt = PQ.coords
        PQxy = PQx * PQy
        PQzt = PQz * PQt

        return ThetaPointDimTwo(
            [
                X * PQzt * PQy,
                Y * PQzt * PQx,
                Z * PQxy * PQt,
                T * PQxy * PQz,
            ]
        )
