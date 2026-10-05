"""Explicit transformation-filter configuration with legacy-compatible defaults."""

from dataclasses import asdict, dataclass, replace
import os


ALIGNMENT_GATE_MODES = frozenset({"native", "common_match_coverage"})


def _value(environ, name, default, converter):
    raw = environ.get(name)
    return default if raw is None else converter(raw)


@dataclass(frozen=True)
class TransformationThresholds:
    """All adjustable gates applied by the transformation stage.

    Defaults mirror the current module constants.  The MultiProt gates are
    separate because MultiProt records use native match/coverage semantics,
    not the TMalign TM-score contract.
    """

    minimum_residue_match_count: int = 15
    minimum_residue_match_percentage: float = 50.0
    minimum_hotspot_match_number: int = 1
    diff_percentage: float = 20.0
    template_residue_count: float = 50.0
    contact_count_threshold: int = 5
    clashing_distance: float = 3.0
    max_clashing_count: int = 5
    scaffold_threshold: float = 5.0
    tm_score_threshold: float = 0.5
    multiprot_minimum_residue_match_count: int = 10
    multiprot_minimum_residue_match_percentage: float = 30.0
    alignment_gate_mode: str = "native"

    def __post_init__(self):
        mode = str(self.alignment_gate_mode).strip().lower()
        if mode not in ALIGNMENT_GATE_MODES:
            allowed = ", ".join(sorted(ALIGNMENT_GATE_MODES))
            raise ValueError(
                f"alignment_gate_mode must be one of {allowed}; got {self.alignment_gate_mode!r}"
            )
        object.__setattr__(self, "alignment_gate_mode", mode)

    @classmethod
    def from_environment(cls, environ=None):
        """Build thresholds from PRISM environment variables."""

        if environ is None:
            environ = os.environ
        contact_default = _value(environ, "PRISM_CONTACT_COUNT", 5, int)
        return cls(
            minimum_residue_match_count=_value(
                environ, "PRISM_MINIMUM_RESIDUE_MATCH_COUNT", 15, int
            ),
            minimum_residue_match_percentage=_value(
                environ, "PRISM_MINIMUM_RESIDUE_MATCH_PERCENTAGE", 50.0, float
            ),
            minimum_hotspot_match_number=_value(
                environ, "PRISM_MINIMUM_HOTSPOT_MATCH_NUMBER", 1, int
            ),
            diff_percentage=_value(environ, "PRISM_DIFF_PERCENTAGE", 20.0, float),
            template_residue_count=_value(
                environ, "PRISM_TEMPLATE_RESIDUE_COUNT", 50.0, float
            ),
            contact_count_threshold=_value(
                environ, "PRISM_CONTACT_COUNT_THRESHOLD", contact_default, int
            ),
            clashing_distance=_value(
                environ, "PRISM_CLASHING_DISTANCE", 3.0, float
            ),
            max_clashing_count=_value(
                environ, "PRISM_MAX_CLASHING_COUNT", 5, int
            ),
            scaffold_threshold=_value(
                environ, "PRISM_SCFF_THRESHOLD", 5.0, float
            ),
            tm_score_threshold=_value(
                environ, "PRISM_TM_SCORE_THRESHOLD", 0.5, float
            ),
            multiprot_minimum_residue_match_count=_value(
                environ, "PRISM_MULTIPROT_MINIMUM_RESIDUE_MATCH_COUNT", 10, int
            ),
            multiprot_minimum_residue_match_percentage=_value(
                environ,
                "PRISM_MULTIPROT_MINIMUM_RESIDUE_MATCH_PERCENTAGE",
                30.0,
                float,
            ),
            alignment_gate_mode=str(
                environ.get("PRISM_ALIGNMENT_GATE_MODE", "native")
            ).strip().lower(),
        )

    def with_overrides(self, overrides=None):
        """Return a copy with named fields replaced."""

        if not overrides:
            return self
        if isinstance(overrides, type(self)):
            return overrides
        unknown = set(overrides) - set(self.__dataclass_fields__)
        if unknown:
            names = ", ".join(sorted(unknown))
            raise ValueError(f"unknown transformation threshold(s): {names}")
        return replace(self, **dict(overrides))

    def as_dict(self):
        return asdict(self)
