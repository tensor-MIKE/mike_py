from .diagonal import diagonal_isogeny
from .dim_four import dim_four_gluing_isogeny
from .scholten import scholten_isogeny


def gluing_isogeny(E, bases, chain_kernel):
    scholten_basis, diag_basis, gluing_basis = bases

    # E x E^sigma -> A2 x A2
    diagonal_kernel, gluing_kernel, image_plus, image_minus = scholten_isogeny(
        E, scholten_basis, diag_basis, gluing_basis, chain_kernel
    )

    # A2 x A2 -> A2 x A2
    gluing_kernel, image_plus, image_minus = diagonal_isogeny(
        diagonal_kernel, gluing_kernel, image_plus, image_minus
    )

    # A2 x A2 -> A4
    return dim_four_gluing_isogeny(gluing_kernel, image_plus, image_minus)
