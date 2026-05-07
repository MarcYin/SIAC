"""Shared 6S native output surface definitions.

The :data:`SIXS_CORE_OUTPUT_SPECS` table maps Python-side output names to
the Fortran *expressions* that compute them. Each entry is a 3-tuple:

    (python_name, fortran_expression, optional_guard)

The ``fortran_expression`` is interpolated verbatim into the generated
Fortran source as ``output_values_out(N) = <expression>`` (see
``siac.algorithms.rt.direct.sixs_build._core_output_assignment_lines``).
Most entries simply name a Fortran variable (``"xa"`` -> ``output_values_out(N) = xa``),
but a handful encode small arithmetic expressions — e.g.
``("sttotr", "sdtotr*sutotr", None)`` writes
``output_values_out(N) = sdtotr*sutotr`` so that ``sttotr`` is the product of
the downward and upward Rayleigh transmittances. This is intentional;
REVIEW.md §3.1 sixs_outputs.py:69 incorrectly flagged it as a name parse
hazard — the consumer is an expression emitter, not a name lookup.

The optional third element is a Fortran condition such as
``"irapp.ge.0"``: when set, the assignment is gated by an
``if`` so the output stays at its initialised value otherwise.
"""

from __future__ import annotations

from typing import Final

# (python_name, fortran_expression, optional_fortran_guard)
SixSCoreOutputSpec = tuple[str, str, str | None]

SIXS_BASE_OUTPUTS: Final[tuple[str, ...]] = ("xap", "xbp", "xcp")

SIXS_CORE_OUTPUT_SPECS: Final[tuple[SixSCoreOutputSpec, ...]] = (
    ("xap", "xap", "irapp.ge.0"),
    ("xa", "xa", "irapp.ge.0"),
    ("xb", "xb", "irapp.ge.0"),
    ("xc", "xc", "irapp.ge.0"),
    ("rapp", "rapp", "irapp.ge.0"),
    ("xrad", "xrad", "irapp.ge.0"),
    ("rog", "rog", "irapp.ge.0"),
    ("rogbrdf", "rogbrdf", "irapp.eq.1"),
    ("refet", "refet", None),
    ("alumet", "alumet", None),
    ("refet1", "refet1", None),
    ("refet2", "refet2", None),
    ("refet3", "refet3", None),
    ("rpfet", "rpfet_report", "ipol.ne.0"),
    ("plumet", "plumet_report", "ipol.ne.0"),
    ("xpol", "xpol_report", "ipol.ne.0"),
    ("rpfet_over_refet", "rpfet_over_refet_report", "ipol.ne.0"),
    ("aini_1_1", "aini(1,1)", None),
    ("aini_1_2", "aini(1,2)", None),
    ("aini_1_3", "aini(1,3)", None),
    ("aini_2_1", "aini(2,1)", None),
    ("aini_2_2", "aini(2,2)", None),
    ("aini_2_3", "aini(2,3)", None),
    ("ainr_1_1", "ainr(1,1)", None),
    ("ainr_1_2", "ainr(1,2)", None),
    ("ainr_1_3", "ainr(1,3)", None),
    ("ainr_2_1", "ainr(2,1)", None),
    ("ainr_2_2", "ainr(2,2)", None),
    ("ainr_2_3", "ainr(2,3)", None),
    ("sb", "sb", None),
    ("seb", "seb", None),
    ("dgasm", "dgasm", None),
    ("ugasm", "ugasm", None),
    ("tgasm", "tgasm", None),
    ("sdwava", "sdwava", None),
    ("suwava", "suwava", None),
    ("stwava", "stwava", None),
    ("sdozon", "sdozon", None),
    ("suozon", "suozon", None),
    ("stozon", "stozon", None),
    ("sddica", "sddica", None),
    ("sudica", "sudica", None),
    ("stdica", "stdica", None),
    ("sdoxyg", "sdoxyg", None),
    ("suoxyg", "suoxyg", None),
    ("stoxyg", "stoxyg", None),
    ("sdniox", "sdniox", None),
    ("suniox", "suniox", None),
    ("stniox", "stniox", None),
    ("sdmeth", "sdmeth", None),
    ("sumeth", "sumeth", None),
    ("stmeth", "stmeth", None),
    ("sdmoca", "sdmoca", None),
    ("sumoca", "sumoca", None),
    ("stmoca", "stmoca", None),
    ("sdtotr", "sdtotr", None),
    ("sutotr", "sutotr", None),
    ("sttotr", "sdtotr*sutotr", None),
    ("sdtota", "sdtota", None),
    ("sutota", "sutota", None),
    ("sttota", "sdtota*sutota", None),
    ("sdtott", "sdtott", None),
    ("sutott", "sutott", None),
    ("sttott", "sdtott*sutott", None),
    ("sasr", "sasr", None),
    ("sasa", "sasa", None),
    ("sast", "sast", None),
    ("sodray", "sodray", None),
    ("sodaer", "sodaer", None),
    ("sodtot", "sodtot", None),
    ("sodrayp", "sodrayp", None),
    ("sodaerp", "sodaerp", None),
    ("sodtotp", "sodtotp", None),
    ("sroray", "sroray", None),
    ("sroaer", "sroaer", None),
    ("srotot", "srotot", None),
    ("srqray", "srqray", "ipol.ne.0"),
    ("srqaer", "srqaer", "ipol.ne.0"),
    ("srqtot", "srqtot", "ipol.ne.0"),
    ("sruray", "sruray", "ipol.ne.0"),
    ("sruaer", "sruaer", "ipol.ne.0"),
    ("srutot", "srutot", "ipol.ne.0"),
    ("srpray", "srpray", "ipol.ne.0"),
    ("srpaer", "srpaer", "ipol.ne.0"),
    ("srptot", "srptot", "ipol.ne.0"),
    ("sdpray", "sdpray", "ipol.ne.0"),
    ("sdpaer", "sdpaer", "ipol.ne.0"),
    ("sdptot", "sdptot", "ipol.ne.0"),
    ("sdppray", "sdppray", "ipol.ne.0"),
    ("sdppaer", "sdppaer", "ipol.ne.0"),
    ("sdpptot", "sdpptot", "ipol.ne.0"),
    ("fophsr", "fophsr", None),
    ("fophsa", "fophsa", None),
    ("fophst", "fophst", None),
    ("foqhsr", "foqhsr", "ipol.ne.0"),
    ("foqhsa", "foqhsa", "ipol.ne.0"),
    ("foqhst", "foqhst", "ipol.ne.0"),
    ("fouhsr", "fouhsr", "ipol.ne.0"),
    ("fouhsa", "fouhsa", "ipol.ne.0"),
    ("fouhst", "fouhst", "ipol.ne.0"),
    ("spdpray", "spdpray", "ipol.ne.0"),
    ("spdpaer", "spdpaer", "ipol.ne.0"),
    ("spdptot", "spdptot", "ipol.ne.0"),
    ("pizerr", "pizerr", None),
    ("pizera", "pizera", None),
    ("pizert", "pizert", None),
    ("rocave", "rocave_report", "idirec.eq.1"),
    ("robar1_over_xnorm1", "robar1_over_xnorm1_report", "idirec.eq.1"),
    ("robar2_over_xnorm2", "robar2_over_xnorm2_report", "idirec.eq.1"),
    ("rbard", "rbard_report", "idirec.eq.1"),
    ("albbrdf", "albbrdf_report", "idirec.eq.1"),
    ("rfoamave", "rfoamave_report", "idirec.eq.1.and.ibrdf.eq.6"),
    ("rwatave", "rwatave_report", "idirec.eq.1.and.ibrdf.eq.6"),
    ("rglitave", "rglitave_report", "idirec.eq.1.and.ibrdf.eq.6"),
    ("rooceaw", "rooceaw", "irapp.ge.0.and.idirec.eq.1.and.ibrdf.eq.6"),
)

SIXS_CORE_OUTPUT_NAMES: Final[tuple[str, ...]] = tuple(
    name for name, _, _ in SIXS_CORE_OUTPUT_SPECS
)

_all_names: list[str] = list(SIXS_BASE_OUTPUTS)
for _name in SIXS_CORE_OUTPUT_NAMES:
    if _name not in _all_names:
        _all_names.append(_name)

SIXS_OUTPUT_VARIABLE_CHOICES: Final[tuple[str, ...]] = tuple(_all_names)
SIXS_DEFAULT_OUTPUT_VARIABLES: Final[tuple[str, ...]] = SIXS_OUTPUT_VARIABLE_CHOICES
SIXS_NATIVE_OUTPUT_NAMES: Final[tuple[str, ...]] = tuple(
    name for name in SIXS_OUTPUT_VARIABLE_CHOICES if name not in {"xbp", "xcp"}
)
