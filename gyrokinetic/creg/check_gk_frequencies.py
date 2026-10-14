"""Check the sheath frequency diagnostic, including restart and optional MPI.

Build rt_gk_sheath_1x2v_p1, then run this script from the repository root:
  python3 gyrokinetic/creg/check_gk_frequencies.py [--mpirun /path/to/mpirun]
Requires Python msgpack. Use LD_LIBRARY_PATH if build libraries are elsewhere.
"""

import argparse
import math
import os
from pathlib import Path
import struct
import subprocess
import tempfile

import msgpack


def read_dynvec(path):
    metadata, records = None, []
    with path.open("rb") as stream:
        while magic := stream.read(5):
            assert magic == b"gkyl0"
            version, kind, meta_size = struct.unpack("=QQQ", stream.read(24))
            assert (version, kind) == (1, 2)
            if meta_size:
                metadata = msgpack.unpackb(stream.read(meta_size))
            dtype, stride, count = struct.unpack("=QQQ", stream.read(24))
            assert (dtype, stride) == (2, 48)
            times = struct.unpack(f"={count}d", stream.read(8 * count))
            values = struct.unpack(f"={6 * count}d", stream.read(stride * count))
            records.extend(zip(times, zip(*(values[i::6] for i in range(6)))))
    assert metadata["numComponents"] == 6
    assert metadata["units"] == "1/s"
    assert metadata["component0"] == "omega_H"
    assert metadata["component1"] == "transit"
    assert [metadata[f"component{i}"] for i in range(6)] == [
        "omega_H", "transit", "transit_max", "transit_min", "omega_d", "nu"]
    assert all(metadata[f"component{i}Definition"] for i in range(6))
    assert records
    # This input has identity geometry, v_parallel,max=4*sqrt(Te/me), Lz=4 m.
    expected_transit = math.sqrt(40 * abs(metadata["charge"]) / metadata["mass"])
    charge, mass = abs(metadata["charge"]), metadata["mass"]
    expected_nu = (1 + math.sqrt(2)) * 0.1 * (6.6 - 0.5*math.log(7e18/1e20)
        + 1.5*math.log(40)) * charge**4 * 7e18 / (
        6*math.sqrt(2)*math.pi**1.5*(8.85418781762038985053656303171075e-12)**2
        * math.sqrt(mass)*(40*charge)**1.5)
    for time, (omega_h, transit, transit_max, transit_min, drift, nu) in records:
        assert math.isfinite(time) and time >= 0
        assert math.isfinite(omega_h) and omega_h > 0
        assert math.isclose(transit_max, expected_transit, rel_tol=1e-12)
        assert 0 < transit < transit_max
        assert transit_min == 0 and drift == 0
        assert math.isclose(nu, expected_nu, rel_tol=1e-8)
        if time == 0:
            # Finite velocity truncation gives a small thermal-moment error.
            assert math.isclose(transit, expected_transit/4, rel_tol=0.02)
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path,
                        default=Path("build/gyrokinetic/creg/rt_gk_sheath_1x2v_p1"))
    parser.add_argument("--mpirun")
    args = parser.parse_args()
    executable = args.executable.resolve()
    env = os.environ.copy()
    build = executable.parents[2]
    env["LD_LIBRARY_PATH"] = ":".join(
        [str(build / solver) for solver in ("core", "moments", "vlasov", "gyrokinetic")]
        + [env.get("LD_LIBRARY_PATH", "")]
    )

    def run(directory, *options, mpi=False):
        command = [str(executable), "-s10000", *options]
        if mpi:
            command = [args.mpirun, "-np", "2", *command, "-M", "-c", "2"]
        result = subprocess.run(command, cwd=directory, env=env, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        assert result.returncode == 0, result.stdout

    filename = executable.name + "-elc_frequencies_max.gkyl"
    with tempfile.TemporaryDirectory() as directory:
        run(directory)
        path = Path(directory) / filename
        serial = read_dynvec(path)
        assert not list(Path(directory).glob("*-ion_frequencies_max.gkyl"))
        original = path.read_bytes()
        run(directory, "-r", "0")
        assert path.read_bytes().startswith(original)
        assert len(read_dynvec(path)) > len(serial)
        # Enabling diagnostics on a restart must create the metadata too.
        path.unlink()
        run(directory, "-r", "0")
        read_dynvec(path)
    if args.mpirun:
        with tempfile.TemporaryDirectory() as directory:
            run(directory, mpi=True)
            parallel = read_dynvec(Path(directory) / filename)
            assert len(serial) == len(parallel)
            for (st, sv), (pt, pv) in zip(serial, parallel):
                assert math.isclose(st, pt, rel_tol=1e-12, abs_tol=1e-20)
                assert all(math.isclose(s, p, rel_tol=1e-12) for s, p in zip(sv, pv))
    print("Frequency metadata, values, opt-out, and restart checks passed."
          + (" Serial/MPI agreement passed." if args.mpirun else ""))


if __name__ == "__main__":
    main()
