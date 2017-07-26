import pytest
from sparcc_fast import sparcc_functions as sf
import pandas as pd
import numpy as np
from numpy.testing import assert_allclose


@pytest.fixture
def frame():
    arr = np.array([[1, 3, 2],
                   [6, 5, 4]], dtype=float)
    r_labels = ['r1', 'r2', ]
    c_labels = ['a', 'b', 'c']
    return pd.DataFrame(arr, columns=c_labels, index=r_labels)


def test_to_fractions(frame):
    # dirichlet random sample
    result = sf.to_fractions(frame)
    assert_allclose(result.sum(axis=1), 1)
    assert frame.index.equals(result.index)
    assert frame.columns.equals(result.columns)


def test_permute_w_replacement(frame):
    perm = sf.permute_w_replacement(frame)
    assert perm.shape == frame.shape
    for c, permvals in perm.iteritems():
        assert set(permvals).issubset(set(frame[c]))
