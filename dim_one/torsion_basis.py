def full_even_point_from_nqr(E):
    """
    Helper function which finds an element suitable
    as an x-coordinate with full even torsion when
    the Montgomery coefficient A is a square in F
    """
    A = E.a2()
    Fp2 = E.base_ring()
    Fp = Fp2.base_ring()
    h = 0
    while True:
        h = h + 1
        if not Fp(1 + h**2).is_square():
            xP = -A / Fp2([1, h])
            if E.is_x_coord(xP):
                return E.lift_x(xP)


def full_even_point_from_A(E):
    """
    Helper function which finds an element suitable
    as an x-coordinate with full even torsion when
    the Montgomery coefficient A is not a square
    """
    A = E.a2()
    xP = 0
    while True:
        xP = xP + A
        if E.is_x_coord(xP):
            return E.lift_x(xP)


def torsion_basis(E, e):
    """
    Based off Algorithm 2.1 of the SQISign spec.

    Compute a deterministic basis for E[2^e]  <P, Q>, assumes A != 0
    """
    A = E.a2()
    assert A != 0

    # Compute a point P of full even order not above (0 : 0)
    if A.is_square():
        P = full_even_point_from_nqr(E)
    else:
        P = full_even_point_from_A(E)

    # Compute a second point of full even torsion not above (0 : 0)
    Q = E.lift_x(-(P.x() + A))

    # Clear the cofactor to recover points of order 2^f
    p = E.base_ring().characteristic()
    cofactor = (p + 1) // (2**e)
    P, Q = cofactor * P, cofactor * Q

    # We set Q = P - Q to ensure Q is above (0 : 0)
    Q = P - Q

    # Ensure we have a good basis
    P2 = 2 ** (e - 1) * P
    Q2 = 2 ** (e - 1) * Q
    assert P2 != Q2
    assert not Q2.x()
    assert not (2 * P2)

    return P, Q
