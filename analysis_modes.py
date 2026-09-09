"""Analysis modes and their output disclosure policies."""

from __future__ import annotations

from enum import Enum


class AnalysisMode(str, Enum):
    """Where exact-match queries originate."""

    USER_SUPPLIED = "user_supplied"
    PROTECTED_CATALOGUE = "protected_catalogue"


class OutputPolicy(str, Enum):
    """The result fields a completed analysis may disclose."""

    FULL = "full"
    RESTRICTED = "restricted"


ANALYSIS_MODE_CHOICES = [
    ("Upload my epitope list (detailed results)", AnalysisMode.USER_SUPPLIED.value),
    (
        "Scan the protected antibody catalogue (restricted results)",
        AnalysisMode.PROTECTED_CATALOGUE.value,
    ),
]


def parse_analysis_mode(value: str) -> AnalysisMode:
    """Parse a UI value without falling back to a less restrictive mode."""

    try:
        return AnalysisMode(value)
    except ValueError as exc:
        raise ValueError("Unknown analysis mode.") from exc


def output_policy_for(mode: AnalysisMode) -> OutputPolicy:
    """Keep disclosure policy selection in one auditable place."""

    if mode is AnalysisMode.USER_SUPPLIED:
        return OutputPolicy.FULL
    if mode is AnalysisMode.PROTECTED_CATALOGUE:
        return OutputPolicy.RESTRICTED
    raise ValueError("Analysis mode has no output policy.")
