

from methods import MidPoint
from kernels.abb import ABBSW, ABBAW

flavours = {
    'sw_mort_model': { 'sw_mort': 'model' },
    'sw_mort_const': { 'sw_mort': 'const' },
}


def abb_init_kernels(L, U, N, flavour):
    """Return SW and AW kernels according to parameters and flavour."""

    sw = ABBSW()
    aw = ABBAW()

    for k in [ sw, aw ]:
        k.L, k.U = L, U
        k.setup(MidPoint(), N)

        flags = flavours[flavour]
        for flag in flags:
            setattr(k, flag, flags[flag])

    return sw, aw