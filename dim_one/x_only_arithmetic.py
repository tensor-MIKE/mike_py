# TODO: we can hand write the xLADDER if we want to show how it "should" be
#       done instead of using P + [x]Q with SageMath arithmetic. The same
#       could be done for the repeated doubling and Jaocbian arithmetic for
#       shared secret generation


def xdbl(X, Z, A24, C24):
    """
    Montgomery doubling with projective curve constant
    """
    t0 = (X + Z) ** 2
    t1 = (X - Z) ** 2
    t2 = t0 - t1
    t1 = t1 * C24
    X_new = t0 * t1
    Z_new = (t2 * A24 + t1) * t2
    return X_new, Z_new


def xdbl_iter(X, Z, A24, C24, n):
    """
    Compute [2^n](X : Z) using a projective curve constant
    """
    for _ in range(n):
        X, Z = xdbl(X, Z, A24, C24)
    return X, Z
