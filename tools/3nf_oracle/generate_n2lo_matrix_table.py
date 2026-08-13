#!/usr/bin/env python3
"""Generate a broad, machine-readable spectator-1 N2LO matrix-element table.

Each signed component is evaluated three ways:

* the production C++ ``chiral_N2LO_3NF_factorized`` implementation;
* the separately implemented Python finite-rank partial-wave projector; and
* a direct five-angle projection of the full Cartesian spin-isospin operator.

The cases deliberately include diagonal and off-diagonal entries, S/P/D pair
waves, spectator D waves, both parities, J=1/2 and 3/2, T=1/2 and 3/2, three
momentum regimes, and explicit Hermitian reverse matrix elements.  The direct
five-angle result is stored at two quadrature orders so its finite-order error
is visible rather than hidden in a single quoted number.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import re
import subprocess
import sys
from typing import Iterable


REPO = Path(__file__).resolve().parents[2]
ORACLE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(ORACLE_DIR))

import factorized_n2lo_pwd as FACTORIZED  # noqa: E402
import full_vector_five_angle_pwd as DIRECT  # noqa: E402


COMPONENTS = ("c1", "c3", "c4", "cD", "cE")
LECS = DIRECT._OP.N2LOLECs(
    c1_gev_inverse=-0.81,
    c3_gev_inverse=-3.2,
    c4_gev_inverse=5.4,
    c_d=-0.2,
    c_e=-0.205,
)
COMPONENT_DRIVER_LECS = {
    "c1": (0.0, 0.0, LECS.c1_gev_inverse, 0.0, 0.0),
    "c3": (0.0, 0.0, 0.0, LECS.c3_gev_inverse, 0.0),
    "c4": (0.0, 0.0, 0.0, 0.0, LECS.c4_gev_inverse),
    "cD": (0.0, LECS.c_d, 0.0, 0.0, 0.0),
    "cE": (LECS.c_e, 0.0, 0.0, 0.0, 0.0),
}


@dataclass(frozen=True)
class SpaceSpec:
    name: str
    pair_j_max: int
    two_total_j_max: int
    isospin_breaking_1s0: bool

    def driver_args(self) -> tuple[str, ...]:
        if (
            self.pair_j_max == 1
            and self.two_total_j_max == 1
            and not self.isospin_breaking_1s0
        ):
            return ()
        return (
            "--space",
            str(self.pair_j_max),
            str(self.two_total_j_max),
            "1" if self.isospin_breaking_1s0 else "0",
        )


@dataclass(frozen=True)
class ChannelKey:
    pair_l: int
    pair_s: int
    pair_j: int
    pair_t: int
    spectator_l: int
    two_spectator_j: int
    two_total_j: int
    two_total_t: int
    parity: int


@dataclass(frozen=True)
class Channel:
    alpha: int
    key: ChannelKey

    def jj(self) -> DIRECT.JjChannel:
        key = self.key
        return DIRECT.JjChannel(
            key.pair_l,
            key.pair_s,
            key.pair_j,
            key.spectator_l,
            key.two_spectator_j,
            key.two_total_j,
            key.pair_t,
            key.two_total_t,
        )


@dataclass(frozen=True)
class BaseCase:
    name: str
    space: str
    bra: ChannelKey
    ket: ChannelKey
    bra_momenta: tuple[float, float]
    ket_momenta: tuple[float, float]
    coverage: tuple[str, ...]


@dataclass(frozen=True)
class MatrixCase:
    name: str
    hermitian_pair: str
    direction: str
    space: str
    bra: ChannelKey
    ket: ChannelKey
    bra_momenta: tuple[float, float]
    ket_momenta: tuple[float, float]
    coverage: tuple[str, ...]


DEFAULT_SPACE = SpaceSpec("J12_T12", 1, 1, False)
EXTENDED_SPACE = SpaceSpec("J12_J32_T12_T32", 2, 3, True)
SPACES = {space.name: space for space in (DEFAULT_SPACE, EXTENDED_SPACE)}


def key(
    pair_l: int,
    pair_s: int,
    pair_j: int,
    pair_t: int,
    spectator_l: int,
    two_spectator_j: int,
    two_total_j: int = 1,
    two_total_t: int = 1,
    parity: int = 1,
) -> ChannelKey:
    return ChannelKey(
        pair_l,
        pair_s,
        pair_j,
        pair_t,
        spectator_l,
        two_spectator_j,
        two_total_j,
        two_total_t,
        parity,
    )


def base_cases() -> tuple[BaseCase, ...]:
    low = ((0.46, 0.39), (0.24, 0.31))
    medium = ((0.83, 0.57), (0.35, 0.72))
    high = ((1.25, 0.92), (0.68, 1.10))
    extended = ((0.71, 0.52), (0.33, 0.64))
    isospin = ((0.92, 0.44), (0.51, 0.76))

    s_singlet = key(0, 0, 0, 1, 0, 1)
    s_triplet = key(0, 1, 1, 0, 0, 1)
    pair_p_singlet = key(1, 0, 1, 0, 1, 1)
    pair_d_triplet = key(2, 1, 1, 0, 0, 1)
    pair_p_triplet = key(1, 1, 1, 1, 1, 1)
    spectator_d_triplet = key(0, 1, 1, 0, 2, 3)
    negative_s_singlet = key(0, 0, 0, 1, 1, 1, parity=-1)
    negative_s_triplet = key(0, 1, 1, 0, 1, 1, parity=-1)
    negative_p_triplet = key(1, 1, 0, 1, 0, 1, parity=-1)

    j32_s_triplet = key(0, 1, 1, 0, 0, 1, two_total_j=3)
    j32_pair_d_triplet = key(2, 1, 1, 0, 0, 1, two_total_j=3)
    j32_singlet_spectator_d = key(0, 0, 0, 1, 2, 3, two_total_j=3)
    t32_j12_s_singlet = key(0, 0, 0, 1, 0, 1, two_total_t=3)
    t32_j32_singlet_spectator_d = key(
        0, 0, 0, 1, 2, 3, two_total_j=3, two_total_t=3
    )

    return (
        BaseCase("s_singlet_diagonal", DEFAULT_SPACE.name, s_singlet, s_singlet,
                 *low, ("diagonal", "S-wave", "c1/c3/cE")),
        BaseCase("s_triplet_diagonal", DEFAULT_SPACE.name, s_triplet, s_triplet,
                 *medium, ("diagonal", "S-wave", "c1/c3/cE")),
        BaseCase("pair_d_from_s", DEFAULT_SPACE.name, pair_d_triplet, s_triplet,
                 *medium, ("off-diagonal", "pair tensor S-D", "c1/c3")),
        BaseCase("pair_p_from_s", DEFAULT_SPACE.name, pair_p_singlet, s_triplet,
                 *low, ("off-diagonal", "P-wave", "c1/c3")),
        BaseCase("p_triplet_from_s", DEFAULT_SPACE.name, pair_p_triplet, s_triplet,
                 *medium, ("off-diagonal", "P-wave", "c4")),
        BaseCase("triplet_s_from_singlet_s", DEFAULT_SPACE.name, s_triplet, s_singlet,
                 *low, ("off-diagonal", "spin/isospin transition", "c4/cD")),
        BaseCase("spectator_d_from_singlet_s", DEFAULT_SPACE.name,
                 spectator_d_triplet, s_singlet, *high,
                 ("off-diagonal", "spectator D-wave", "c4/cD")),
        BaseCase("negative_parity_s_transition", DEFAULT_SPACE.name,
                 negative_s_triplet, negative_s_singlet, *low,
                 ("off-diagonal", "negative parity", "spectator P-wave", "c4/cD")),
        BaseCase("negative_parity_p_diagonal", DEFAULT_SPACE.name,
                 negative_p_triplet, negative_p_triplet, *medium,
                 ("diagonal", "negative parity", "pair P-wave", "signed c1/c3")),
        BaseCase("j32_pair_d_from_s", EXTENDED_SPACE.name,
                 j32_pair_d_triplet, j32_s_triplet, *extended,
                 ("off-diagonal", "J=3/2", "pair tensor S-D", "c1/c3")),
        BaseCase("j32_s_from_singlet_spectator_d", EXTENDED_SPACE.name,
                 j32_s_triplet, j32_singlet_spectator_d, *extended,
                 ("off-diagonal", "J=3/2", "spectator D-wave", "c4/cD")),
        BaseCase("t32_j12_s_diagonal", EXTENDED_SPACE.name,
                 t32_j12_s_singlet, t32_j12_s_singlet, *isospin,
                 ("diagonal", "T=3/2", "J=1/2", "c1/c3/cE")),
        BaseCase("t32_j32_spectator_d_diagonal", EXTENDED_SPACE.name,
                 t32_j32_singlet_spectator_d, t32_j32_singlet_spectator_d,
                 *isospin, ("diagonal", "T=3/2", "J=3/2", "spectator D-wave")),
    )


def matrix_cases() -> tuple[MatrixCase, ...]:
    cases: list[MatrixCase] = []
    for base in base_cases():
        for direction in ("forward", "reverse"):
            reverse = direction == "reverse"
            cases.append(MatrixCase(
                name=f"{base.name}:{direction}",
                hermitian_pair=base.name,
                direction=direction,
                space=base.space,
                bra=base.ket if reverse else base.bra,
                ket=base.bra if reverse else base.ket,
                bra_momenta=base.ket_momenta if reverse else base.bra_momenta,
                ket_momenta=base.bra_momenta if reverse else base.ket_momenta,
                coverage=base.coverage,
            ))
    return tuple(cases)


CHANNEL_PATTERN = re.compile(
    r"^# alpha (\d+)  L2=(\d+) S2=(\d+) J2=(\d+) T2=(\d+) "
    r"l1=(\d+) 2j1=(\d+) 2J3=(\d+) 2T3=(\d+) P3=(-?\d+)$",
    re.MULTILINE,
)


def driver_command(
    driver: Path,
    component_lecs: tuple[float, float, float, float, float],
    cutoff_mev: float,
    transfer_order: int,
    space: SpaceSpec,
) -> list[str]:
    c_e, c_d, c1, c3, c4 = component_lecs
    return [
        str(driver), "--factorized", str(c_e), str(c_d), str(c1), str(c3),
        str(c4), str(cutoff_mev), str(transfer_order), *space.driver_args(),
    ]


def invoke_driver(command: list[str], stdin: str = "") -> str:
    result = subprocess.run(
        command,
        input=stdin,
        text=True,
        capture_output=True,
        timeout=900,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"production driver failed with status {result.returncode}:\n"
            f"{' '.join(command)}\n{result.stderr}"
        )
    return result.stdout


def discover_channels(
    driver: Path,
    cutoff_mev: float,
    transfer_order: int,
    space: SpaceSpec,
) -> dict[ChannelKey, Channel]:
    stdout = invoke_driver(
        driver_command(driver, (0.0, 0.0, 0.0, 0.0, 0.0), cutoff_mev,
                       transfer_order, space)
    )
    channels: dict[ChannelKey, Channel] = {}
    for match in CHANNEL_PATTERN.finditer(stdout):
        alpha, pair_l, pair_s, pair_j, pair_t, spectator_l, two_spectator_j, \
            two_total_j, two_total_t, parity = map(int, match.groups())
        channel_key = ChannelKey(
            pair_l, pair_s, pair_j, pair_t, spectator_l, two_spectator_j,
            two_total_j, two_total_t, parity,
        )
        if channel_key in channels:
            raise RuntimeError(f"duplicate production channel {channel_key}")
        channels[channel_key] = Channel(alpha, channel_key)
    if not channels:
        raise RuntimeError("production driver emitted no parseable channel table")
    return channels


def run_production_batch(
    driver: Path,
    cutoff_mev: float,
    transfer_order: int,
    space: SpaceSpec,
    component: str,
    cases: Iterable[MatrixCase],
    channels: dict[ChannelKey, Channel],
) -> dict[str, float]:
    case_list = tuple(cases)
    lines = []
    for case in case_list:
        bra = channels[case.bra]
        ket = channels[case.ket]
        p_bra, q_bra = case.bra_momenta
        p_ket, q_ket = case.ket_momenta
        lines.append(
            f"{bra.alpha} {ket.alpha} {p_bra:.17g} {q_bra:.17g} "
            f"{p_ket:.17g} {q_ket:.17g}\n"
        )
    stdout = invoke_driver(
        driver_command(driver, COMPONENT_DRIVER_LECS[component], cutoff_mev,
                       transfer_order, space),
        "".join(lines),
    )
    values = [
        float(line.split()[-1])
        for line in stdout.splitlines()
        if line.startswith("W1 ")
    ]
    if len(values) != len(case_list):
        raise RuntimeError(
            f"driver returned {len(values)} values for {len(case_list)} "
            f"{space.name}/{component} cases"
        )
    return {case.name: value for case, value in zip(case_list, values)}


def complex_json(value: complex | float) -> dict[str, float]:
    z = complex(value)
    return {"real": float(z.real), "imag": float(z.imag)}


def regulator_product(
    momenta: tuple[float, float, float, float],
    cutoff_mev: float,
    constants: DIRECT._OP.N2LOConstants,
) -> float:
    cutoff = cutoff_mev / constants.hbar_c_mev_fm
    p, q, p_prime, q_prime = momenta

    def regulator(pair: float, spectator: float) -> float:
        invariant = pair * pair + 0.75 * spectator * spectator
        return math.exp(-((invariant / (cutoff * cutoff)) ** 2))

    return regulator(p, q) * regulator(p_prime, q_prime)


def oracle_values(
    case: MatrixCase,
    channels: dict[ChannelKey, Channel],
    cutoff_mev: float,
    factorized_order: int,
    direct_order_low: int,
    direct_order_high: int,
    constants: DIRECT._OP.N2LOConstants,
    projector: DIRECT.FiveAngleProjector,
) -> tuple[dict[str, complex], dict[str, complex], dict[str, complex]]:
    bra = channels[case.bra].jj()
    ket = channels[case.ket].jj()
    p_bra, q_bra = case.bra_momenta
    p_ket, q_ket = case.ket_momenta
    # Both Python projectors use (ket p,q, bra p',q').
    momenta = (p_ket, q_ket, p_bra, q_bra)
    factorized = FACTORIZED.project_n2lo_jj_recoupled(
        bra, ket, momenta, constants, LECS, factorized_order, cutoff_mev
    )
    regulator = regulator_product(momenta, cutoff_mev, constants)

    def direct_at(order: int) -> dict[str, complex]:
        raw = projector.project_jj_direct(bra, ket, momenta, LECS, order)
        normalized = projector.to_tictac_normalization(raw)
        return {name: value * regulator for name, value in normalized.items()}

    return factorized, direct_at(direct_order_low), direct_at(direct_order_high)


def component_record(
    production_low: float,
    production_high: float,
    factorized: complex,
    direct_low: complex,
    direct_high: complex,
) -> dict[str, object]:
    production_factorized_error = float(abs(production_high - factorized))
    production_direct_error = float(abs(production_high - direct_high))
    production_step = float(abs(production_high - production_low))
    direct_step = float(abs(direct_high - direct_low))
    magnitude = float(max(abs(production_high), abs(factorized), abs(direct_high)))
    classification = "selection-rule zero" if magnitude < 1.0e-10 else "nonzero"

    checks = {
        "production_matches_factorized_oracle": bool(
            production_factorized_error <= 5.0e-11 + 5.0e-10 * abs(factorized)
        ),
        "production_transfer_quadrature_stable": bool(
            production_step <= 1.0e-8 + 1.0e-5 * magnitude
        ),
        # The direct projector is deliberately non-factorized but finite-order.
        # Consistency is judged against twice its last observed quadrature step.
        "independent_full_vector_consistent": bool(
            production_direct_error <= 5.0e-7 + 2.0 * direct_step
        ),
        "independent_full_vector_quadrature_resolved": bool(
            direct_step <= 2.0e-5 + 5.0e-3 * magnitude
        ),
        "imaginary_leakage_small": bool(
            max(abs(factorized.imag), abs(direct_high.imag)) <= 5.0e-10
        ),
    }
    return {
        "classification": classification,
        "production_low_order": complex_json(production_low),
        "production_high_order": complex_json(production_high),
        "factorized_python_oracle": complex_json(factorized),
        "independent_full_vector_low_order": complex_json(direct_low),
        "independent_full_vector_high_order": complex_json(direct_high),
        "absolute_errors": {
            "production_vs_factorized": production_factorized_error,
            "production_vs_independent_full_vector": production_direct_error,
            "production_last_order_step": production_step,
            "independent_full_vector_last_order_step": direct_step,
        },
        "checks": checks,
        "passed": all(checks.values()),
    }


def generate_report(
    driver: Path,
    cutoff_mev: float,
    production_order_low: int,
    production_order_high: int,
    direct_order_low: int,
    direct_order_high: int,
) -> dict[str, object]:
    cases = matrix_cases()
    constants = DIRECT._OP.N2LOConstants.tictac()
    projector = DIRECT.FiveAngleProjector(constants)
    channel_maps = {
        name: discover_channels(driver, cutoff_mev, production_order_high, space)
        for name, space in SPACES.items()
    }
    for case in cases:
        if case.bra not in channel_maps[case.space]:
            raise RuntimeError(f"missing bra channel for {case.name}: {case.bra}")
        if case.ket not in channel_maps[case.space]:
            raise RuntimeError(f"missing ket channel for {case.name}: {case.ket}")

    production: dict[tuple[int, str], dict[str, float]] = {}
    for order in (production_order_low, production_order_high):
        for component in COMPONENTS:
            merged: dict[str, float] = {}
            for space_name, space in SPACES.items():
                subset = tuple(case for case in cases if case.space == space_name)
                merged.update(run_production_batch(
                    driver, cutoff_mev, order, space, component, subset,
                    channel_maps[space_name],
                ))
            production[(order, component)] = merged
            print(f"production order={order} component={component}: {len(merged)} values",
                  flush=True)

    rows: list[dict[str, object]] = []
    for index, case in enumerate(cases, start=1):
        factorized, direct_low, direct_high = oracle_values(
            case,
            channel_maps[case.space],
            cutoff_mev,
            production_order_high,
            direct_order_low,
            direct_order_high,
            constants,
            projector,
        )
        components = {
            component: component_record(
                production[(production_order_low, component)][case.name],
                production[(production_order_high, component)][case.name],
                factorized[component],
                direct_low[component],
                direct_high[component],
            )
            for component in COMPONENTS
        }
        channels = channel_maps[case.space]
        rows.append({
            "id": case.name,
            "hermitian_pair": case.hermitian_pair,
            "direction": case.direction,
            "space": case.space,
            "coverage": list(case.coverage),
            "bra": {"alpha": channels[case.bra].alpha, **asdict(case.bra)},
            "ket": {"alpha": channels[case.ket].alpha, **asdict(case.ket)},
            "momenta_fm_inverse": {
                "bra": {"p": case.bra_momenta[0], "q": case.bra_momenta[1]},
                "ket": {"p": case.ket_momenta[0], "q": case.ket_momenta[1]},
            },
            "components": components,
            "passed_before_hermiticity": all(
                component["passed"] for component in components.values()
            ),
        })
        print(f"independent projection {index}/{len(cases)}: {case.name}", flush=True)

    hermiticity: list[dict[str, object]] = []
    rows_by_pair: dict[str, dict[str, dict[str, object]]] = {}
    for row in rows:
        rows_by_pair.setdefault(str(row["hermitian_pair"]), {})[
            str(row["direction"])
        ] = row
    for pair_name, directions in rows_by_pair.items():
        forward = directions["forward"]
        reverse = directions["reverse"]
        component_checks: dict[str, object] = {}
        for component in COMPONENTS:
            f_comp = forward["components"][component]
            r_comp = reverse["components"][component]
            source_errors = {}
            for source in (
                "production_high_order",
                "factorized_python_oracle",
                "independent_full_vector_high_order",
            ):
                f_value = complex(**f_comp[source])
                r_value = complex(**r_comp[source])
                source_errors[source] = float(abs(f_value - r_value.conjugate()))
            direct_tolerance = 5.0e-7 + 2.0 * (
                f_comp["absolute_errors"]["independent_full_vector_last_order_step"]
                + r_comp["absolute_errors"]["independent_full_vector_last_order_step"]
            )
            passed = (
                source_errors["production_high_order"] <= 5.0e-11
                and source_errors["factorized_python_oracle"] <= 5.0e-11
                and source_errors["independent_full_vector_high_order"]
                <= direct_tolerance
            )
            component_checks[component] = {
                "absolute_errors": source_errors,
                "independent_full_vector_tolerance": direct_tolerance,
                "passed": passed,
            }
        pair_passed = all(item["passed"] for item in component_checks.values())
        hermiticity.append({
            "pair": pair_name,
            "components": component_checks,
            "passed": pair_passed,
        })
        forward["hermiticity_passed"] = pair_passed
        reverse["hermiticity_passed"] = pair_passed
        forward["passed"] = forward["passed_before_hermiticity"] and pair_passed
        reverse["passed"] = reverse["passed_before_hermiticity"] and pair_passed

    component_records = [
        component
        for row in rows
        for component in row["components"].values()
    ]
    summary = {
        "matrix_elements": len(rows),
        "component_comparisons": len(component_records),
        "nonzero_component_comparisons": sum(
            item["classification"] == "nonzero" for item in component_records
        ),
        "selection_rule_zero_comparisons": sum(
            item["classification"] == "selection-rule zero"
            for item in component_records
        ),
        "hermitian_pairs": len(hermiticity),
        "passed_matrix_elements": sum(row["passed"] for row in rows),
        "passed_hermitian_pairs": sum(item["passed"] for item in hermiticity),
    }
    summary["passed"] = (
        summary["passed_matrix_elements"] == summary["matrix_elements"]
        and summary["passed_hermitian_pairs"] == summary["hermitian_pairs"]
    )

    return {
        "schema": "tictac.n2lo_3nf_matrix_validation.v1",
        "description": (
            "Signed spectator-1 W^(1) partial-wave matrix elements in Tic-tac "
            "normalization, including the nonlocal regulator on bra and ket."
        ),
        "settings": {
            "lecs": asdict(LECS),
            "cutoff_mev": cutoff_mev,
            "production_transfer_orders": [
                production_order_low, production_order_high,
            ],
            "factorized_python_order": production_order_high,
            "independent_full_vector_orders": [direct_order_low, direct_order_high],
            "fourier_normalization": "(2*pi)^-6",
            "momentum_units": "fm^-1",
            "matrix_element_units": "fm^5",
            "driver": str(driver.relative_to(REPO) if driver.is_relative_to(REPO) else driver),
        },
        "threshold_policy": {
            "production_vs_factorized": "abs <= 5e-11 + 5e-10*|oracle|",
            "production_transfer_step": "abs <= 1e-8 + 1e-5*|element|",
            "production_vs_independent": (
                "abs <= 5e-7 + 2*|independent(high)-independent(low)|"
            ),
            "independent_quadrature_step": "abs <= 2e-5 + 5e-3*|element|",
            "hermiticity": (
                "production/factorized abs <= 5e-11; independent abs <= "
                "5e-7 + 2*(forward last step + reverse last step)"
            ),
        },
        "spaces": [asdict(space) for space in SPACES.values()],
        "summary": summary,
        "hermiticity": hermiticity,
        "matrix_elements": rows,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--driver",
        type=Path,
        default=REPO / "build" / "tools" / "3nf_oracle" / "print_w1_element",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO / "output" / "validation" / "n2lo_3nf_matrix_elements.json",
    )
    parser.add_argument("--cutoff-mev", type=float, default=450.0)
    parser.add_argument("--production-order-low", type=int, default=8)
    parser.add_argument("--production-order", type=int, default=10)
    parser.add_argument("--direct-order-low", type=int, default=6)
    parser.add_argument("--direct-order", type=int, default=8)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    driver = args.driver.resolve()
    if not driver.is_file():
        raise SystemExit(
            f"driver not found: {driver}\n"
            "build it with: cmake --build build --target print_w1_element"
        )
    orders = (
        args.production_order_low,
        args.production_order,
        args.direct_order_low,
        args.direct_order,
    )
    if any(order < 1 for order in orders):
        raise SystemExit("all quadrature orders must be positive")
    if args.production_order_low >= args.production_order:
        raise SystemExit("production low order must be smaller than high order")
    if args.direct_order_low >= args.direct_order:
        raise SystemExit("direct low order must be smaller than high order")

    report = generate_report(
        driver,
        args.cutoff_mev,
        args.production_order_low,
        args.production_order,
        args.direct_order_low,
        args.direct_order,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["summary"], indent=2))
    print(f"wrote {args.output}")
    return 0 if report["summary"]["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
