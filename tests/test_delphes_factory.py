import itertools

import cppyy


def test_permanent_array_survives_clear(load_delphes):
    _, factory, _ = load_delphes({})
    array = factory.NewPermanentArray()
    member = factory.NewCandidate()
    member.PID = 211
    array.Add(member)
    address = cppyy.addressof(array)

    factory.Clear()

    assert cppyy.addressof(array) == address
    assert array.GetEntries() == 0

    fresh_member = factory.NewCandidate()
    fresh_member.PID = 13

    array.Add(fresh_member)

    assert array.GetEntries() == 1
    assert array.At(0).PID == 13


def test_non_permanent_array_slot_is_recycled(load_delphes):
    _, factory, _ = load_delphes({})
    old = factory.NewArray()
    old_addr = cppyy.addressof(old)

    factory.Clear()

    new = factory.NewArray()

    assert cppyy.addressof(new) == old_addr
    assert new.GetEntries() == 0


def test_clear_resets_candidate_fields_to_defaults(load_delphes):
    _, factory, _ = load_delphes({})
    first = factory.NewCandidate()
    first.PID = 211
    first.Charge = 2
    first.Status = 3
    first.Momentum.SetPtEtaPhiE(50.0, 0.5, 0.3, 60.0)
    first.Position.SetXYZT(1.0, 2.0, 3.0, 4.0)
    first.D0 = 0.1
    first.DZ = 0.2
    first.CtgTheta = 0.5
    first.PT = 42.0
    first.Eta = 0.7
    first.Rapidity = 1.3
    address = cppyy.addressof(first)

    factory.Clear()

    fresh = factory.NewCandidate()
    assert cppyy.addressof(fresh) == address
    assert fresh.PID == 0
    assert fresh.Charge == 0
    assert fresh.Status == 0
    assert fresh.Momentum.Px() == 0.0
    assert fresh.Momentum.Py() == 0.0
    assert fresh.Momentum.Pz() == 0.0
    assert fresh.Momentum.E() == 0.0
    assert fresh.Position.X() == 0.0
    assert fresh.Position.Y() == 0.0
    assert fresh.Position.Z() == 0.0
    assert fresh.Position.T() == 0.0
    assert fresh.D0 == 0.0
    assert fresh.DZ == 0.0
    assert fresh.CtgTheta == 0.0
    assert fresh.PT == 0.0
    assert fresh.Eta == 0.7
    assert fresh.Rapidity == 1.3


def test_candidate_ids_monotonic_within_event(load_delphes):
    _, factory, _ = load_delphes({})
    factory.Clear()
    first_id = factory.NewCandidate().GetUniqueID()
    ids = [first_id] + [factory.NewCandidate().GetUniqueID() for _ in range(200)]

    assert len(set(ids)) == 201
    assert all(a < b for a, b in itertools.pairwise(ids))
    assert all(b - a == 1 for a, b in itertools.pairwise(ids))
    assert first_id == 1


def test_candidate_ids_restart_after_clear(load_delphes):
    _, factory, _ = load_delphes({})
    factory.Clear()
    before = [factory.NewCandidate().GetUniqueID() for _ in range(3)]

    factory.Clear()

    after = [factory.NewCandidate().GetUniqueID() for _ in range(3)]

    assert before == [1, 2, 3]
    assert after == [1, 2, 3]


def test_many_allocations_do_not_break_pool_invariants(load_delphes):
    _, factory, _ = load_delphes({})
    factory.Clear()
    ids = set()
    for i in range(100):
        c = factory.NewCandidate()
        assert c.GetUniqueID() not in ids
        ids.add(c.GetUniqueID())
        if i % 10 == 0:
            factory.NewArray()
        if i % 25 == 0:
            factory.NewPermanentArray()

    factory.Clear()

    new_ids = [factory.NewCandidate().GetUniqueID() for _ in range(10)]

    assert len(set(new_ids)) == 10
    assert new_ids[0] == 1
