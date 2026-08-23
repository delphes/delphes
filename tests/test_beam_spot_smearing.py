import pytest
from conftest import C_LIGHT_MM_PER_NS, assert_deterministic, make_candidate, make_config


def beam_spot_config(**extra):
    return make_config("BeamSpotSmearing", **extra)


def run_beam_spot(load_delphes, config, particles):
    module, factory, _ = load_delphes(config)
    input_array = module.ExportArray("inputParticles")
    for pos, decay in particles:
        c = make_candidate(factory, 50.0, 0.5)
        c.Position.SetXYZT(*pos)
        c.DecayPosition.SetXYZT(*decay)
        input_array.Add(c)
    module.Init()
    module.Process()
    return [
        (
            c.Position.X(),
            c.Position.Y(),
            c.Position.Z(),
            c.Position.T(),
            c.DecayPosition.X(),
            c.DecayPosition.Y(),
            c.DecayPosition.Z(),
            c.DecayPosition.T(),
        )
        for c in (input_array.At(i) for i in range(input_array.GetEntries()))
    ]


PARTICLES = [
    ((0.0, 0.0, 0.0, 0.0), (10.0, 20.0, 30.0, 40.0)),
    ((1.0, 2.0, 3.0, 4.0), (11.0, 22.0, 33.0, 44.0)),
]


def test_all_sigmas_zero_leaves_positions_unchanged(load_delphes):
    config = beam_spot_config(SigmaX=0.0, SigmaY=0.0, SigmaZ=0.0, SigmaT=0.0)
    result = run_beam_spot(load_delphes, config, PARTICLES)
    for pos, decay in PARTICLES:
        out_pos = result[PARTICLES.index((pos, decay))][:4]
        out_decay = result[PARTICLES.index((pos, decay))][4:]
        assert out_pos == pytest.approx(pos, abs=1e-12)
        assert out_decay == pytest.approx(decay, abs=1e-12)


def test_rigid_shift_same_for_all_candidates(load_delphes):
    config = beam_spot_config(SigmaX=0.001, SigmaY=0.002, SigmaZ=0.003, SigmaT=0.0)
    result = run_beam_spot(load_delphes, config, PARTICLES)

    shifts = []
    for (pos, decay), out in zip(PARTICLES, result):
        out_pos = out[:4]
        out_decay = out[4:]

        shift_pos = tuple(b - a for a, b in zip(pos, out_pos))
        shift_decay = tuple(b - a for a, b in zip(decay, out_decay))
        assert shift_pos == pytest.approx(shift_decay, abs=1e-9)
        shifts.append(shift_pos)

    for shift in shifts[1:]:
        assert shift == pytest.approx(shifts[0], abs=1e-9)

    dx, dy, dz, dt = shifts[0]
    assert abs(dx) <= 5.0 * 1.0
    assert abs(dy) <= 5.0 * 2.0
    assert abs(dz) <= 5.0 * 3.0
    assert dt == pytest.approx(0.0, abs=1e-9)


def test_sigma_t_unit_conversion_and_scaling(load_delphes):
    result_1ns = run_beam_spot(load_delphes, beam_spot_config(SigmaT=1.0e-9), PARTICLES)
    result_2ns = run_beam_spot(load_delphes, beam_spot_config(SigmaT=2.0e-9), PARTICLES)

    for (pos, decay), out1, out2 in zip(PARTICLES, result_1ns, result_2ns):
        dt1 = out1[3] - pos[3]
        dt2 = out2[3] - pos[3]

        assert dt2 == pytest.approx(2.0 * dt1, rel=1e-6)
        assert abs(dt1) <= 5.0 * C_LIGHT_MM_PER_NS

        assert out1[:3] == pytest.approx(pos[:3], abs=1e-12)
        assert out1[4:7] == pytest.approx(decay[:3], abs=1e-12)


def test_deterministic_with_fixed_seed(load_delphes):
    config = beam_spot_config(SigmaX=0.001, SigmaY=0.001, SigmaZ=0.001, SigmaT=1.0e-9)
    assert_deterministic(
        lambda: run_beam_spot(load_delphes, config, PARTICLES),
        extract=lambda results: [value for row in results for value in row],
        abs_tol=1e-12,
    )


def test_empty_input(load_delphes):
    config = beam_spot_config(SigmaX=0.001)
    result = run_beam_spot(load_delphes, config, [])
    assert result == []
