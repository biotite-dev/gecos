# This source code is part of the Gecos package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

import pytest
import inspect
import gecos
import gecos.cli as gecli


def test_default_args_consistency():
    """
    Check whether the default arguments of :meth:`ColorOptimzer.optimize()`
    match the default parameters of the CLI.
    """
    signature = inspect.signature(
        gecos.ColorOptimizer.optimize
    )
    ref_args = {
        key: val.default for key, val in signature.parameters.items()
        if val.default is not inspect.Parameter.empty
    }
    
    # Mock that reports whether the given arguments are the defaults
    def mock_optimize(self, n_steps,
                 beta_start, beta_end, stepsize_start, stepsize_end):
        nonlocal ref_args
        assert n_steps        == ref_args["n_steps"]
        assert beta_start     == ref_args["beta_start"]
        assert beta_end       == ref_args["beta_end"]
        assert stepsize_start == ref_args["stepsize_start"]
        assert stepsize_end   == ref_args["stepsize_end"]
    
    try:
        # Temporarily replace color optimizer with mock
        temp = gecos.ColorOptimizer.optimize
        gecos.ColorOptimizer.optimize = mock_optimize
        gecli.main(args=[
            "--seed", "0",
            "--nruns", "1"
        ])
    finally:
        gecos.ColorOptimizer.optimize = temp
