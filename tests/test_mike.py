import json
from pathlib import Path

import pytest

from mike.mike import Mike, PARAMS_I, PARAMS_III, PARAMS_V

KAT_PATH = Path(__file__).parent / "test_data.json"
KAT_DATA = json.loads(KAT_PATH.read_text())

PARAMS = {
    "PARAMS_I": PARAMS_I,
    "PARAMS_III": PARAMS_III,
    "PARAMS_V": PARAMS_V,
}


def load_party(entry):
    return {
        "seed": bytes.fromhex(entry["seed"]),
        "pk": bytes.fromhex(entry["pk"]),
        "sk": bytes.fromhex(entry["sk"]),
        "secret": bytes.fromhex(entry["secret"]),
    }


TEST_CASES = {
    name: {
        "params": PARAMS[name],
        "alice": load_party(party["alice"]),
        "bob": load_party(party["bob"]),
    }
    for name, party in KAT_DATA.items()
}


@pytest.mark.parametrize("name", TEST_CASES, ids=str.lower)
def test_keygen_known_answers(name):
    case = TEST_CASES[name]
    params = case["params"]

    for party in ("alice", "bob"):
        data = case[party]
        m = Mike(params, data["seed"])
        pk, sk = m.keygen()
        assert pk == data["pk"]
        assert sk == data["sk"]


@pytest.mark.parametrize("name", TEST_CASES, ids=str.lower)
def test_shared_secret_known_answer(name):
    case = TEST_CASES[name]
    params = case["params"]
    alice_data = case["alice"]
    bob_data = case["bob"]

    alice = Mike(params, alice_data["seed"])
    bob = Mike(params, bob_data["seed"])

    alice_pk, slice_sk = alice.keygen()
    bob_pk, bob_sk = bob.keygen()

    alice_secret = alice.shared_secret(slice_sk, bob_pk)
    bob_secret = bob.shared_secret(bob_sk, alice_pk)
    assert alice_secret == bob_secret
    assert alice_secret == alice_data["secret"]
