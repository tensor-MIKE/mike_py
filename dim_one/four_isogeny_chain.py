from .x_only_arithmetic import xdbl_iter


def four_isogeny_codomain(X, Z):
    """
    Compute the codomain of the 4-isogeny E -> E/<P> such that [2]P != (0 : 1)
    together with three constants (c0, c1, c2) used for efficient evaluation
    """
    c1 = X - Z
    c2 = X + Z

    t = 2 * Z**2
    C24 = t**2
    c0 = 2 * t

    A24 = (2 * X**2) ** 2

    return A24, C24, c0, c1, c2


def four_isogeny_eval(c0, c1, c2, X, Z):
    """
    Evaluate a point Q (XQ : ZQ) under the action of the 4-isogeny E -> E/<P>
    for [2]P != (0 : 1) using precomputed constants (c0, c1, c2)
    """
    s = X + Z
    d = X - Z

    X1 = s * c1
    Z1 = d * c2
    u = s * d * c0

    v = (X1 + Z1) ** 2
    w = (X1 - Z1) ** 2

    X_new = (u + v) * v
    Z_new = w * (w - u)
    return X_new, Z_new


def four_isogeny_chain(domain, kernel, e):
    """
    Compute a 2^e-isogeny from `domain` as a chain of 4-isogenies, using a
    balanced strategy, and return the image of the four-torsion above the
    kernel on the codomain.
    """
    assert e % 2 == 0, "e must be even"

    A = domain.a2()
    A24 = (A + 2) / 4
    C24 = domain.base_ring().one()

    X, Z = kernel[0], kernel[2]

    # Stack of (point, order) pairs, order tracked in a parallel list.
    # The top of the stack is always the most recently doubled point.
    strategy = [(X, Z)]
    orders = [e + 2]

    for _ in range(e // 2):
        # Double the top of the stack until it has order exactly 4.
        while orders[-1] != 2:
            prev_X, prev_Z = strategy[-1]
            prev_order = orders[-1]
            m = 2 * (prev_order // 4) + (prev_order & 1)

            strategy.append(xdbl_iter(prev_X, prev_Z, A24, C24, m))
            orders.append(prev_order - m)

        # Pop the order-4 point and use it as the kernel of a 4-isogeny.
        kernel_X, kernel_Z = strategy.pop()
        orders.pop()
        A24, C24, c0, c1, c2 = four_isogeny_codomain(kernel_X, kernel_Z)

        # Push every remaining point through the isogeny.
        strategy = [four_isogeny_eval(c0, c1, c2, X, Z) for X, Z in strategy]
        orders = [order - 2 for order in orders]

    ((X, Z),) = strategy
    return X, Z
