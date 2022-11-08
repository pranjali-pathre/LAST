#!/usr/bin/env python

"""Pick Testing
"""

import numpy as np
from pick import rmsd

def test_rmsd():
    """Test RMSD
    """
    p_1 = np.array([[1.0,0.35,3.1],[1.06,0.08,2.7]])
    p_2 = np.array([[0.06,0.98,1.0],[8.06,0.7,4.5]])

    assert (rmsd(p_1, p_2, val=2) - 5.399763883726769) < 1e-9

test_rmsd()
