# Tests

All test cases and benchmarks live in this directory, grouped by module.
Random-polynomial generation APIs remain in their implementation modules.

With Make:

```sh
make test-minor-dp
make test-components
./build/component_tests --list
./build/component_tests rational-solver
```

With CMake (`BUILD_TESTING=ON`, the default):

```sh
cmake --build <build-dir> --target det_minor_dp_test component_tests
ctest --test-dir <build-dir> -R 'drsolve_(minor_dp|component_)' --output-on-failure
```

`det_minor_dp.c` contains automated determinant and scheduling regressions.
The `*_test.c` files extracted from implementation modules retain their
original diagnostic output; many are examples or benchmarks rather than
assertion-based tests. The component runner exposes them individually.
The sparse interpolation case and roots benchmark are opt-in because they
can be expensive. Existing CLI suites such as `--test 1` and `--test 5`
also perform large computations and are not part of the smoke tests above.

`dixon_test.c` and `polynomial_system_solver_test.c` remain linked into the
library to preserve existing CLI test entry points and exported functions.
The former also provides random-system generation used by the CLI's `-r`
mode. The existing `./drsolve --test <n>` interface is unchanged.
