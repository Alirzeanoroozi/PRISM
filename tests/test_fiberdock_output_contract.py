
from pathlib import Path

from benchmark.scripts.diagnose_fiberdock_output_contract import (
    declared_value,
    parse_energy,
)


def test_declared_fiberdock_energy_stem_and_solution(tmp_path: Path) -> None:
    params = tmp_path / "fd_params.txt"
    ref = tmp_path / "fiberdock_energies.ref"
    params.write_text(
        f"energiesOutFileName {tmp_path / 'fiberdock_energies'}\n",
        encoding="utf-8",
    )
    ref.write_text("1 | 0.51 | -15.12 | 11.67 |\n", encoding="utf-8")

    stem = declared_value(params, "energiesOutFileName")
    assert stem == str(tmp_path / "fiberdock_energies")
    assert not (tmp_path / "fd_params.ref").exists()

    solution = parse_energy(Path(stem + ".ref"))
    assert solution is not None
    assert solution["solution"] == 1
    assert solution["global_energy"] == 0.51
