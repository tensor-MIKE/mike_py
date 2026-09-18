from dim_two.isogeny import two_isogeny


def diagonal_isogeny(diagonal_kernel, gluing_kernel, image_plus, image_minus):
    n = len(image_plus)
    assert len(image_minus) == n
    assert len(gluing_kernel) == 3

    # First of each pair into U, second into V.
    all_points = list(gluing_kernel) + list(image_plus) + list(image_minus)
    U_images = [u for (u, _) in all_points]
    V_images = [v for (_, v) in all_points]

    U1, U2, V1, V2 = diagonal_kernel
    codomain_u, U_images = two_isogeny(U1, U2, U_images, [True, False])
    codomain_v, V_images = two_isogeny(V1, V2, V_images, [False, True])

    # Extend the gluing basis with a final differential addition
    u_21 = codomain_u.diff_addition(U_images[2], U_images[1], U_images[0])
    v_21 = codomain_v.diff_addition(V_images[2], V_images[1], V_images[0])

    gluing_kernel = [
        [U_images[0], V_images[0]],
        [U_images[1], V_images[1]],
        [U_images[2], V_images[2]],
        [u_21, v_21],
    ]
    image_plus = [[U_images[3 + i], V_images[3 + i]] for i in range(n)]
    image_minus = [[U_images[3 + n + i], V_images[3 + n + i]] for i in range(n)]

    return gluing_kernel, image_plus, image_minus
