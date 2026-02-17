import math
from dataclasses import dataclass

R_KCAL_PER_MOL_K = 1.987204258e-3  # kcal/(mol*K)

@dataclass(frozen=True)
class KdResult:
    kd_M: float
    kd_nM: float
    kd_uM: float
    kd_mM: float

def dg_to_kd(dg_kcal_mol: float, T_K: float = 298.15) -> KdResult:
    """
    Convert standard binding free energy ΔG° (kcal/mol) to Kd (M).
    Uses: ΔG° = RT ln(Kd)  =>  Kd = exp(ΔG° / (RT))
    """
    RT = R_KCAL_PER_MOL_K * T_K
    kd_M = math.exp(dg_kcal_mol / RT)
    return KdResult(
        kd_M=kd_M,
        kd_nM=kd_M * 1e9,
        kd_uM=kd_M * 1e6,
        kd_mM=kd_M * 1e3,
    )

def kd_to_dg(kd_M: float, T_K: float = 298.15) -> float:
    """
    Convert Kd (M) to standard binding free energy ΔG° (kcal/mol).
    Uses: ΔG° = RT ln(Kd)
    """
    if kd_M <= 0:
        raise ValueError("Kd must be > 0")
    RT = R_KCAL_PER_MOL_K * T_K
    return RT * math.log(kd_M)

def pretty_kd(kd_M: float) -> str:
    """
    Human-friendly formatting for Kd.
    """
    if kd_M < 1e-12:
        return f"{kd_M:.3e} M"
    if kd_M < 1e-9:
        return f"{kd_M*1e12:.3g} pM"
    if kd_M < 1e-6:
        return f"{kd_M*1e9:.3g} nM"
    if kd_M < 1e-3:
        return f"{kd_M*1e6:.3g} µM"
    if kd_M < 1:
        return f"{kd_M*1e3:.3g} mM"
    return f"{kd_M:.3g} M"

if __name__ == "__main__":
    # Example: ΔG = -9.5 kcal/mol at 298 K ~ ~100 nM-ish
    dg = -9.5
    res = dg_to_kd(dg, T_K=298.15)
    print(f"ΔG° = {dg:.2f} kcal/mol -> Kd = {res.kd_M:.3e} M ({pretty_kd(res.kd_M)})")

    # Reverse conversion
    kd = 1e-7
    dg_back = kd_to_dg(kd, T_K=298.15)
    print(f"Kd = {kd:.1e} M -> ΔG° = {dg_back:.2f} kcal/mol")