from __future__ import annotations

import math

from siac.sixs_upstream_parity import (
    default_parity_cases,
    parse_original_sixs_stdout,
    render_original_sixs_input,
)


def test_parse_original_sixs_stdout_extracts_lambertian_values() -> None:
    stdout = """
*       apparent reflectance  0.0242139  appar. rad.(w/m2/sr/mic)    9.911    *
*                   total gaseous transmittance  0.960                        *
*           wv above aerosol :   0.024     wv mixed with aerosol :   0.024    *
*                       wv under aerosol :   0.024                            *
*       app. polarized refl.  0.0007    app. pol. rad. (w/m2/sr/mic)    0.246 *
*             direction of the plane of polarization-43.00                    *
*                   total polarization ratio     0.027                        *
*     % of direct  irr.    % of diffuse irr.    % of enviro. irr              *
*               0.875               0.125               0.000                 *
*     atm. intrin. ref.   background  ref.  pixel  reflectance                *
*               0.024               0.000               0.000                 *
*     direct solar irr.    atm. diffuse irr.    environment  irr              *
*            1036.755             148.566               0.000                 *
*     atm. intrin. rad.    background  rad.    pixel  radiance                *
*               9.911               0.000               0.000                 *
*          int. funct filter (in mic)              int. sol. spect (in w/m2)  *
*             1.0000000                                1484.817               *
*      global gas. trans. :     0.97777        0.98059        0.95954         *
*      spherical albedo   :     0.04063        0.03943        0.07191         *
*      reflectance I      :     0.01773        0.00703        0.02504         *
*      degree of polar.   :       10.22          18.02           2.67         *
*       input apparent reflectance            :    0.100                      *
*       measured radiance [w/m2/sr/mic]       :   40.931                      *
*       Lambertian case :      0.08746                                        *
*       BRDF       case :      0.08746                                        *
*       coefficients xa xb xc                 :  0.00284  0.02812  0.07191    *
*       coefficients xap xb xc                :  1.161395  0.028122  0.071913 *
"""

    parsed = parse_original_sixs_stdout(stdout)

    assert parsed["refet"].value == 0.0242139
    assert parsed["alumet"].value == 9.911
    assert parsed["tgasm"].value == 0.95954
    assert parsed["refet3"].value == 0.024
    assert parsed["rpfet"].value == 0.0007
    assert parsed["xpol"].value == -43.0
    assert parsed["aini_1_1"].value == 0.875
    assert parsed["ainr_2_1"].value == 9.911
    assert parsed["sb"].value == 1.0
    assert parsed["seb"].value == 1484.817
    assert parsed["dgasm"].value == 0.97777
    assert parsed["sast"].value == 0.07191
    assert parsed["sdptot"].value == 2.67
    assert parsed["rapp"].value == 0.1
    assert parsed["xrad"].value == 40.931
    assert parsed["rog"].value == 0.08746
    assert parsed["xap"].value == 1.161395
    assert parsed["xb"].value == 0.028122
    assert parsed["xc"].value == 0.071913


def test_parse_original_sixs_stdout_extracts_brdf_and_ocean_values() -> None:
    stdout = """
*       Foam:        0.00019  Water:       0.00086  Glint:       0.00881      *
*                    rodir    robar    ropbar    robarbar  albedo             *
*                 0.00986   0.01806   0.02118       NaN   0.06395             *
  roocean water    7.40198269E-02
"""

    parsed = parse_original_sixs_stdout(stdout)

    assert parsed["rfoamave"].value == 0.00019
    assert parsed["rwatave"].value == 0.00086
    assert parsed["rglitave"].value == 0.00881
    assert parsed["rocave"].value == 0.00986
    assert parsed["robar1_over_xnorm1"].value == 0.01806
    assert parsed["robar2_over_xnorm2"].value == 0.02118
    assert math.isnan(parsed["rbard"].value)
    assert parsed["albbrdf"].value == 0.06395
    assert parsed["rooceaw"].value == 7.40198269e-02


def test_render_original_sixs_input_uses_same_basic_structure() -> None:
    case = default_parity_cases()[0]

    deck = render_original_sixs_input(case)
    lines = deck.strip().splitlines()

    assert lines[0] == "0"
    assert lines[2] == "8"
    assert lines[3] == "2 0.3"
    assert lines[4] == "1"
    assert lines[-2] == "0"
    assert lines[-1] == "-0.1"
