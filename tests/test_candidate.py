import math

import cppyy
import pytest
import ROOT


def test_new_candidate_unique_ids(load_delphes):
    _, factory, _ = load_delphes({})
    c1 = factory.NewCandidate()
    c2 = factory.NewCandidate()
    assert c1.GetUniqueID() != c2.GetUniqueID()
    assert c1.GetUniqueID() > 0


def test_candidate_copy(load_delphes):
    _, factory, _ = load_delphes({})
    source = factory.NewCandidate()
    source.PID = 211
    source.Charge = 1
    source.Status = 3
    source.Momentum.SetPtEtaPhiE(50.0, 0.5, 0.0, 56.0)
    target = factory.NewCandidate()
    source.Copy(target)
    assert target.PID == 211
    assert target.Charge == 1
    assert target.Status == 3
    assert target.Momentum.Pt() == pytest.approx(50.0, rel=1e-6)
    assert target.Momentum.Eta() == pytest.approx(0.5, rel=1e-6)


COPY_SCALAR_FIELDS = [
    ("PID", 211),
    ("Status", 3),
    ("M1", 1),
    ("M2", 2),
    ("D1", 3),
    ("D2", 4),
    ("Charge", 1),
    ("IsPU", 1),
    ("IsRecoPU", 1),
    ("IsConstituent", 1),
    ("IsFromConversion", 1),
    ("ClusterIndex", 7),
    ("ClusterNDF", 5),
    ("NTimeHits", 9),
    ("NCharged", 2),
    ("NNeutrals", 3),
    ("NSubJetsTrimmed", 1),
    ("NSubJetsPruned", 2),
    ("NSubJetsSoftDropped", 3),
    ("Flavor", 5),
    ("FlavorAlgo", 6),
    ("FlavorPhys", 7),
    ("TauFlavor", 8),
    ("BTag", 0x010),
    ("BTagAlgo", 0x020),
    ("BTagPhys", 0x040),
    ("TauTag", 0x100),
    ("BoostedTag", 0x200),
    ("Mass", 0.135),
    ("TauWeight", 0.42),
    ("Eem", 11.0),
    ("Ehad", 22.0),
    ("Etrk", 33.0),
    ("DeltaEta", 0.05),
    ("DeltaPhi", 0.06),
    ("L", 1.23),
    ("DZ", 0.045),
    ("ErrorDZ", 0.001),
    ("ErrorT", 0.002),
    ("D0", 0.012),
    ("ErrorD0", 0.003),
    ("C", 0.02),
    ("ErrorC", 0.004),
    ("P", 50.5),
    ("ErrorP", 0.5),
    ("PT", 51.5),
    ("ErrorPT", 0.6),
    ("CtgTheta", 1.1),
    ("ErrorCtgTheta", 0.01),
    ("Phi", 1.23),
    ("ErrorPhi", 0.02),
    ("Nclusters", 12.5),
    ("dNdx", 2.5),
    ("Xd", 0.1),
    ("Yd", 0.2),
    ("Zd", 0.3),
    ("XFirstHit", 0.4),
    ("YFirstHit", 0.5),
    ("ZFirstHit", 0.6),
    ("TrackResolution", 0.05),
    ("Beta", 0.25),
    ("BetaStar", 0.35),
    ("MeanSqDeltaR", 0.45),
    ("PTD", 0.55),
    ("NeutralEnergyFraction", 0.65),
    ("ChargedEnergyFraction", 0.75),
    ("IsolationVar", 0.85),
    ("IsolationVarRhoCorr", 0.95),
    ("SumPtCharged", 1.05),
    ("SumPtNeutral", 1.15),
    ("SumPtChargedPU", 1.25),
    ("SumPt", 1.35),
    ("ParticleDensity", 7.75),
    ("ClusterSigma", 2.5),
    ("SumPT2", 3.5),
    ("BTVSumPT2", 4.5),
    ("GenDeltaZ", 5.5),
    ("GenSumPT2", 6.5),
    ("ExclYmerge12", 11.0),
    ("ExclYmerge23", 12.0),
    ("ExclYmerge34", 13.0),
    ("ExclYmerge45", 14.0),
    ("ExclYmerge56", 15.0),
]

COPY_VECTOR_FIELDS = [
    "Momentum",
    "Position",
    "InitialPosition",
    "DecayPosition",
    "PositionError",
    "Area",
    "SoftDroppedJet",
    "SoftDroppedSubJet1",
    "SoftDroppedSubJet2",
]

COPY_VECTOR_ARRAY_FIELDS = ["TrimmedP4", "PrunedP4", "SoftDroppedP4"]


COPY_VECTOR_ARRAY_OFFSETS = {"TrimmedP4": 3.0, "PrunedP4": 4.0, "SoftDroppedP4": 5.0}


def set_xyzt(vec, i, offset):
    vec.SetXYZT(float(i) + offset, float(i) + offset + 0.1, float(i) + offset + 0.2, 100.0 + i + offset)


def test_candidate_copy_full_fields(load_delphes):
    _, factory, _ = load_delphes({})
    source = factory.NewCandidate()

    for name, value in COPY_SCALAR_FIELDS:
        setattr(source, name, value)

    source.Edges[0] = 0.5
    source.Edges[1] = -0.5
    source.Edges[2] = 1.5
    source.Edges[3] = -1.5
    source.FracPt[0] = 0.1
    source.FracPt[1] = 0.2
    source.FracPt[2] = 0.3
    source.FracPt[3] = 0.4
    source.FracPt[4] = 0.5
    source.Tau[0] = 0.01
    source.Tau[1] = 0.02
    source.Tau[2] = 0.03
    source.Tau[3] = 0.04
    source.Tau[4] = 0.05

    for i, name in enumerate(COPY_VECTOR_FIELDS):
        vec = getattr(source, name)
        vec.SetXYZT(1.0 * (i + 1), 2.0 * (i + 1), 3.0 * (i + 1), 10.0 * (i + 1))

    for name in COPY_VECTOR_ARRAY_FIELDS:
        for i in range(5):
            set_xyzt(getattr(source, name)[i], i, COPY_VECTOR_ARRAY_OFFSETS[name])

    cov = source.TrackCovariance
    cov.ResizeTo(6, 6)
    for i in range(6):
        cov[i][i] = float(i + 1)
    cov[1][0] = 0.1
    cov[2][0] = 0.2
    cov[4][1] = 0.3

    source.ECalEnergyTimePairs.push_back((1.5, 2.5))
    source.ECalEnergyTimePairs.push_back((3.5, 4.5))

    child1 = factory.NewCandidate()
    child1.PID = 91
    child2 = factory.NewCandidate()
    child2.PID = 111
    source.AddCandidate(child1)
    source.AddCandidate(child2)

    target = factory.NewCandidate()
    source.Copy(target)

    for name, value in COPY_SCALAR_FIELDS:
        got = getattr(target, name)
        if isinstance(value, int):
            assert got == value, f"field {name} not copied"
        else:
            assert got == pytest.approx(value, rel=1e-6, abs=1e-9), f"field {name} not copied"

    for i in range(4):
        assert target.Edges[i] == source.Edges[i], f"Edges[{i}] not copied"
        assert target.FracPt[i] == source.FracPt[i], f"FracPt[{i}] not copied"
        assert target.Tau[i] == source.Tau[i], f"Tau[{i}] not copied"

    for name in COPY_VECTOR_FIELDS:
        s, t = getattr(source, name), getattr(target, name)
        assert (t.X(), t.Y(), t.Z(), t.T()) == pytest.approx(
            (s.X(), s.Y(), s.Z(), s.T()), rel=1e-9, abs=1e-12
        ), f"{name} not copied"

    for name in COPY_VECTOR_ARRAY_FIELDS:
        for i in range(5):
            s, t = getattr(source, name)[i], getattr(target, name)[i]
            assert (t.X(), t.Y(), t.Z(), t.T()) == pytest.approx(
                (s.X(), s.Y(), s.Z(), s.T()), rel=1e-9, abs=1e-12
            ), f"{name}[{i}] not copied"

    tcov = target.TrackCovariance
    assert tcov.GetNrows() == 6 and tcov.GetNcols() == 6
    for i in range(6):
        assert tcov(i, i) == pytest.approx(float(i + 1), rel=1e-12), f"cov ({i},{i}) not copied"
    for r, c, v in ((1, 0, 0.1), (2, 0, 0.2), (4, 1, 0.3)):
        assert tcov(r, c) == pytest.approx(v, rel=1e-12), f"cov ({r},{c}) not copied"

    tp = target.ECalEnergyTimePairs
    assert tp.size() == 2
    for i, (e, t) in enumerate(((1.5, 2.5), (3.5, 4.5))):
        assert tp.at(i).first == pytest.approx(e, rel=1e-6)
        assert tp.at(i).second == pytest.approx(t, rel=1e-6)

    tchildren = target.GetCandidates()
    assert tchildren.GetEntries() == 2
    assert tchildren.At(0).GetUniqueID() == child1.GetUniqueID()
    assert tchildren.At(1).GetUniqueID() == child2.GetUniqueID()

    extra = factory.NewCandidate()
    extra.PID = 321
    target.AddCandidate(extra)
    tchildren = target.GetCandidates()
    assert tchildren.GetEntries() == 3
    assert tchildren.At(2).GetUniqueID() == extra.GetUniqueID()


def test_add_and_get_candidates(load_delphes):
    _, factory, _ = load_delphes({})
    parent = factory.NewCandidate()
    child = factory.NewCandidate()
    child.PID = 999
    parent.AddCandidate(child)
    candidates = parent.GetCandidates()
    assert candidates.GetEntries() == 1
    assert candidates.At(0).PID == 999

    assert parent.GetCandidates() is candidates


def test_factory_clear_reuses_memory_and_ids(load_delphes):
    _, factory, _ = load_delphes({})
    first = factory.NewCandidate()
    first.PID = 211
    factory.Clear()

    second = factory.NewCandidate()
    third = factory.NewCandidate()
    assert first.PID == 0

    assert second.GetUniqueID() != third.GetUniqueID()


def test_sort_candidates_by_pt(load_delphes):
    _, factory, _ = load_delphes({})
    array = factory.NewArray()
    for pt in (30.0, 10.0, 50.0):
        c = factory.NewCandidate()
        c.Momentum.SetPtEtaPhiE(pt, 0.0, 0.0, pt)
        array.Add(c)

    assert array.At(0).IsSortable()
    array.Sort()
    assert [array.At(i).Momentum.Pt() for i in range(3)] == [50.0, 30.0, 10.0]


def test_sort_candidates_by_sumpt2_comparator(load_delphes):
    _, factory, _ = load_delphes({})
    array = factory.NewArray()
    for pt, sumpt2 in ((30.0, 3.0), (10.0, 9.0), (50.0, 1.0)):
        c = factory.NewCandidate()
        c.Momentum.SetPtEtaPhiE(pt, 0.0, 0.0, pt)
        c.SumPT2 = sumpt2
        array.Add(c)

    comp = ROOT.CompSumPT2[ROOT.Candidate].Instance()
    ROOT.Candidate.fgCompare = comp
    array.Sort()
    assert [array.At(i).SumPT2 for i in range(3)] == [9.0, 3.0, 1.0]

    ROOT.Candidate.fgCompare = ROOT.CompMomentumPt[ROOT.Candidate].Instance()


def test_candidate_kinematics_helpers_massless(load_delphes):
    _, factory, _ = load_delphes({})
    c = factory.NewCandidate()
    pt, eta, phi = 40.0, 1.0, 0.5
    e = pt * math.cosh(eta)
    c.Momentum.SetPtEtaPhiE(pt, eta, phi, e)
    assert c.Momentum.Pt() == pytest.approx(pt, rel=1e-12)
    assert c.Momentum.Eta() == pytest.approx(eta, rel=1e-12)
    assert c.Momentum.Phi() == pytest.approx(phi, rel=1e-12)
    assert c.Momentum.E() == pytest.approx(e, rel=1e-12)

    assert c.Momentum.M() == pytest.approx(0.0, abs=1e-9)
    assert c.Momentum.M2() == pytest.approx(0.0, abs=1e-9)
    assert c.Momentum.P() == pytest.approx(e, rel=1e-12)


def test_candidate_kinematics_helpers_massive(load_delphes):
    _, factory, _ = load_delphes({})
    c = factory.NewCandidate()
    px, py, pz, e = 3.0, 4.0, 12.0, 13.0
    c.Momentum.SetXYZT(px, py, pz, e)
    p = (px * px + py * py + pz * pz) ** 0.5
    m2 = e * e - p * p
    assert c.Momentum.Pt() == pytest.approx(5.0, rel=1e-12)
    assert c.Momentum.E() == pytest.approx(e, rel=1e-12)
    assert c.Momentum.M2() == pytest.approx(m2, rel=1e-12)
    assert c.Momentum.M() == pytest.approx(m2**0.5, rel=1e-12)
    assert c.Momentum.Eta() == pytest.approx(0.5 * math.log((p + pz) / (p - pz)), rel=1e-12)


def test_candidate_zero_momentum_no_nan(load_delphes):
    _, factory, _ = load_delphes({})
    c = factory.NewCandidate()

    assert math.isfinite(c.Momentum.Eta())
    assert math.isfinite(c.Momentum.Phi())
    assert math.isfinite(c.Momentum.P())
    c.Momentum.SetPtEtaPhiE(0.0, 0.0, 0.0, 0.0)
    assert math.isfinite(c.Momentum.Eta())
    assert math.isfinite(c.Momentum.Phi())
    assert c.Momentum.Pt() == pytest.approx(0.0, abs=1e-12)
    assert c.Momentum.P() == pytest.approx(0.0, abs=1e-12)
    assert c.Momentum.M2() == pytest.approx(0.0, abs=1e-12)


def test_candidate_children_reference_resolution(load_delphes):
    _, factory, _ = load_delphes({})
    parent = factory.NewCandidate()
    child = factory.NewCandidate()
    child.PID = 999
    parent.AddCandidate(child)
    assert parent.GetCandidates().GetEntries() == 1
    entry = parent.GetCandidates().At(0)
    assert entry.GetUniqueID() == child.GetUniqueID()
    assert entry.PID == 999


def test_candidate_pool_reuse_after_clear(load_delphes):
    _, factory, _ = load_delphes({})
    first = factory.NewCandidate()
    second_pre = factory.NewCandidate()
    first_addr = cppyy.addressof(first)
    second_addr = cppyy.addressof(second_pre)
    factory.Clear()
    new_first = factory.NewCandidate()
    new_second = factory.NewCandidate()

    assert cppyy.addressof(new_first) == first_addr
    assert cppyy.addressof(new_second) == second_addr

    new_first.PID = 211
    assert new_second.PID == 0
    assert new_first.GetUniqueID() != new_second.GetUniqueID()


@pytest.mark.parametrize(
    "class_name, field",
    [
        ("Photon", "PT"),
        ("Electron", "PT"),
        ("Muon", "PT"),
        ("Jet", "PT"),
        ("Track", "PT"),
        ("Tower", "E"),
        ("Vertex", "SumPT2"),
    ],
)
def test_default_comparator_per_subclass(load_delphes, class_name, field):
    _, factory, _ = load_delphes({})
    klass = getattr(ROOT, class_name)
    values = [30.0, 10.0, 50.0]
    array = factory.NewArray()
    for value in values:
        obj = factory.New(klass.Class())
        setattr(obj, field, value)
        assert obj.IsSortable(), f"{class_name} should be sortable by default"
        array.Add(obj)
    array.Sort()
    assert [array.At(i).__getattribute__(field) for i in range(3)] == [50.0, 30.0, 10.0]


def test_gen_particle_not_sortable_by_default(load_delphes):
    _, factory, _ = load_delphes({})
    obj = factory.New(ROOT.GenParticle.Class())
    assert not obj.IsSortable()


def make_tower(factory, et, e):
    tower = factory.New(ROOT.Tower.Class())
    tower.ET = et
    tower.E = e
    return tower


def test_comp_et_orders_towers_by_transverse_energy(load_delphes):
    _, factory, _ = load_delphes({})
    comp = ROOT.CompET[ROOT.Tower].Instance()
    low = make_tower(factory, 10.0, 30.0)
    high = make_tower(factory, 30.0, 10.0)
    same = make_tower(factory, 10.0, 20.0)
    assert comp.Compare(low, high) == 1
    assert comp.Compare(high, low) == -1
    assert comp.Compare(low, same) == 0


def test_comp_e_orders_towers_by_energy(load_delphes):
    _, factory, _ = load_delphes({})
    comp = ROOT.CompE[ROOT.Tower].Instance()
    low = make_tower(factory, 30.0, 5.0)
    high = make_tower(factory, 10.0, 7.0)
    same = make_tower(factory, 20.0, 5.0)
    assert comp.Compare(low, high) == 1
    assert comp.Compare(high, low) == -1
    assert comp.Compare(low, same) == 0


def make_jet_at(factory, eta, phi):
    jet = factory.New(ROOT.Jet.Class())
    jet.Eta = eta
    jet.Phi = phi
    return jet


def test_comp_delta_r_orders_by_reference_distance(load_delphes):
    _, factory, _ = load_delphes({})
    comp = ROOT.CompDeltaR[ROOT.Jet, ROOT.Jet].Instance()
    reference = make_jet_at(factory, 0.0, 0.0)
    comp.SetObject(reference)

    near = make_jet_at(factory, 0.0, 3.0)
    wrapped = make_jet_at(factory, 0.0, -3.1)
    further = make_jet_at(factory, 1.0, 0.0)

    assert comp.Compare(near, wrapped) == -1
    assert comp.Compare(wrapped, near) == 1

    assert comp.Compare(near, further) == 1
    assert comp.Compare(further, near) == -1
    assert comp.Compare(near, near) == 0


def test_sort_jets_by_delta_r_to_reference(load_delphes):
    _, factory, _ = load_delphes({})
    comp = ROOT.CompDeltaR[ROOT.Jet, ROOT.Jet].Instance()
    reference = make_jet_at(factory, 0.0, 0.0)
    comp.SetObject(reference)
    array = factory.NewArray()
    for eta, phi in ((1.0, 0.0), (0.0, -3.1), (0.0, 3.0)):
        array.Add(make_jet_at(factory, eta, phi))

    ROOT.Jet.fgCompare = comp
    array.Sort()
    ROOT.Jet.fgCompare = ROOT.CompPT[ROOT.Jet].Instance()

    phis = [array.At(i).Phi for i in range(3)]
    assert phis[0] == pytest.approx(0.0, abs=1e-6)
    assert phis[1] == pytest.approx(3.0, abs=1e-6)
    assert phis[2] == pytest.approx(-3.1, abs=1e-6)
