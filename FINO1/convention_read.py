import numpy as np

def decompose(x: np.float32):
    """Decomposes a float32 into negative, exponent, and significand

    See https://stackoverflow.com/questions/46093123 for context.

    """
    negative = x < 0
    n = np.abs(x).view(np.int32)  # discard sign (MSB now 0),
                                  # view bit string as int32
    exponent = (n >> 23) - 127  # drop significand, correct exponent offset
                                # 23 and 127 are specific to float32
    significand = n & np.int32(2**23 - 1)  # second factor provides mask
                                           # to extract significand
    return (negative, exponent, significand)


def accuracy(x: np.float32) -> int:
    """Extracts the accuracy of a float32 as an integer power of two"""
    negative, exponent, significand = decompose(x)
    b = significand != 0
    nz_significand = significand[b]
    reverse_index = 23 * np.ones(x.shape)
    reverse_index[b] = np.log2(nz_significand & -nz_significand)
        # for v & -v trick, see https://stackoverflow.com/questions/18806481
    return reverse_index - 23 + exponent
