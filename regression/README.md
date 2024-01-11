Regression tests
----------------

These are whole-simulation tests, testing the system's ability to perform a given simulation.

To create/run a new regression test:

1. Write and save your regression test file in this folder (with a `rt_` prefix).
2. In the main directory re-compile your regression test with `make -j # <build_dir>/regression/<test_file_name>`, where # is some responsibly chosen number of cores, `<build_dir>` is either `build` or `cuda-build` (use the latter if on a GPU-enabled machine).
3. Run the test's executable in `build/unit` (or in `cuda-build/unit` if running on GPU).

To see what runtime options are available, including how to run on the GPU, use `<executable> -h`.
