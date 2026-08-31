from backend.project_adaptation_v2 import QualifiedEvidencePair, fit_project_adapter


def events(n, *, duplicate=False):
    return [QualifiedEvidencePair(str(i), i, "c1ccccc1" + "C" * i, "sol", 1.0, {"A": 1.8, "B": 1.0}, duplicate_status="SAME_MEASUREMENT" if duplicate else "DISTINCT_MEASUREMENT") for i in range(n)]


def test_under_five_is_base_only_and_weights_are_simplex():
    result = fit_project_adapter("sol", events(4), {"A": .8, "B": .2})
    assert result.status == "BASE_ONLY" and result.project_weights == result.global_weights
    assert sum(result.project_weights.values()) == 1 and min(result.project_weights.values()) >= 0


def test_sequential_project_series_shifts_only_after_gate_and_loo():
    early = fit_project_adapter("sol", events(5), {"A": .8, "B": .2})
    assert early.status == "LIGHT_PROJECT_ADAPTATION"
    assert early.project_weights["B"] > early.global_weights["B"]
    result = fit_project_adapter("sol", events(12), {"A": .8, "B": .2})
    assert result.status == "REGULARIZED_PROJECT_ENSEMBLE"
    assert result.project_weights["B"] > result.global_weights["B"]
    assert result.adapted_validation_error <= result.base_validation_error


def test_duplicates_and_external_quality_do_not_inflate_effective_n():
    result = fit_project_adapter("sol", events(10, duplicate=True), {"A": .5, "B": .5})
    assert result.raw_n == 0 and result.status == "BASE_ONLY"
    external = QualifiedEvidencePair("x", 1, "CC", "sol", 1, {"A": 1, "B": 1}, origin="EXPERIMENTAL_EXTERNAL", source_quality="D")
    assert fit_project_adapter("sol", [external], {"A": .5, "B": .5}).raw_n == 0


def test_same_compound_is_excluded_for_query_similarity_fit():
    source = events(6)
    result = fit_project_adapter("sol", source, {"A": .5, "B": .5}, query_smiles=source[0].smiles)
    assert result.raw_n == 5
