def add_components(E, P, Q):
    """
    For distinct points P, Q on the Montgomery curve y^2 = x^3 + A*x^2 + x,
    return (u, v, w) with

        x(P + Q) = (u - v) / w
        x(P - Q) = (u + v) / w

    NOTE:

    The C/Rust version works with Jacobian coordinates (x = X/Z^2, y = Y/Z^3),
    but being lazy and using SageMath points means we actually have normalised
    points with Z = 1, so the weighting between Jacobian or Projective doesn't
    matter for now.

    If Z != 1 though, we need to be more careful. TODO for future Jack I suppose.
    """
    assert P[2] == 1 and Q[2] == 1, "add_components needs affine (Z = 1) points"
    A = E.a2()

    xP, yP = P[0], P[1]
    xQ, yQ = Q[0], Q[1]

    dx = xP - xQ
    w = dx * dx
    v = 2 * (yP * yQ)
    u = yP * yP + yQ * yQ - (A + xP + xQ) * w

    return u, v, w
