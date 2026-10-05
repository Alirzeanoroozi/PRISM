"""Input and template-panel selection helpers for the current PRISM CLI."""

from pathlib import Path
import re


_TEMPLATE_ID = re.compile(r"^[A-Za-z0-9]{6}$")


def normalize_template_ids(values):
    """Normalize template IDs from CLI tokens or a plain-text manifest.

    Current PRISM template IDs are four-character PDB IDs followed by two
    chain IDs.  Comma-separated CLI values are accepted for notebook and shell
    convenience, while ordering and the first occurrence of each ID are kept.
    """

    if isinstance(values, str):
        values = [values]

    normalized = []
    seen = set()
    for value in values or []:
        for token in str(value).replace(",", " ").split():
            if not _TEMPLATE_ID.fullmatch(token):
                raise ValueError(
                    f"invalid template ID {token!r}; expected a six-character "
                    "PDB-plus-chain ID such as 1abcAB"
                )
            if token not in seen:
                normalized.append(token)
                seen.add(token)
    return normalized


def read_template_list(path):
    """Read a six-character template ID per non-comment line."""

    path = Path(path)
    values = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        values.append(line.split()[0])
    return normalize_template_ids(values)


def select_templates(
    default_templates,
    *,
    explicit_templates=None,
    template_list_path=None,
    template_limit=None,
):
    """Select the template panel while preserving the current default.

    ``explicit_templates`` takes precedence over a list path; callers should
    normally enforce mutual exclusion at the CLI layer.  A limit is applied
    after selection, so it behaves consistently for generated and precomputed
    panels.
    """

    if explicit_templates is not None and template_list_path is not None:
        raise ValueError("use either explicit templates or --template-list, not both")
    if template_limit is not None and template_limit <= 0:
        raise ValueError("--template-limit must be positive")

    if explicit_templates is not None:
        selected = normalize_template_ids(explicit_templates)
    elif template_list_path is not None:
        selected = read_template_list(template_list_path)
    else:
        selected = normalize_template_ids(default_templates)

    if template_limit is not None:
        selected = selected[:template_limit]
    return selected
