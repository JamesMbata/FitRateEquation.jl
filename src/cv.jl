# Cross-validation folds. Every fit is cross-validated under TWO fold schemes, both written
# to the run dir (goodness_of_fit.csv summary + cv_folds.csv per-fold losses):
#
#   :article — leave-one-ARTICLE-out (the HEADLINE estimate). Folds by article (the part of
#              the Article|Fig group key before '|'), so sibling figures from one paper —
#              which share a lab, an enzyme prep and an assay — never leak across train/test.
#   :figure  — leave-one-FIGURE-out (one fold per full Article|Fig group key). Optimistic
#              relative to :article (a held-out figure's sibling figures stay in training);
#              it exists for a like-for-like comparison with EnzymeRates.jl, whose LOOCV holds
#              out one `group` (= Article_Fig, the same partition) at a time, includes 1-row
#              groups, walks groups in data order, and reports SE = std/sqrt(n_groups).
#              `_group_folds` reproduces each of those choices.
#
# Held-out rows are scored with the same per-group centered log-ratio loss at the train-fit
# params (÷ held-out row count); a 1-row group carries no centered signal and scores 0 unless
# it takes the sign/finite penalty. Absolute mode fits a single article, so its :article
# scheme has no non-degenerate fold and reports NaN with n_folds = 0; :figure is the genuine
# held-out CV there.

"The fold schemes every fit runs, in report order (headline first)."
const CV_SCHEMES = (:article, :figure)

_article(group::AbstractString) = String(split(group, '|')[1])

"Train/test row-index folds, one per article (in first-appearance order), degenerate folds dropped."
function _article_folds(d::Dataset)
    arts = _article.(d.group)
    folds = NamedTuple[]
    for a in unique(arts)
        test  = findall(==(a), arts)
        train = findall(!=(a), arts)
        (isempty(test) || isempty(train)) && continue
        push!(folds, (held_out=a, train=train, test=test))
    end
    folds
end

# Carry Et through the subset so absolute-mode CV folds keep each row's enzyme concentration
# (byte-identical for relative mode: Et is all-NaN there, so d.Et[idx] == fill(NaN,...)).
_subset(d::Dataset, idx) =
    Dataset(d.concs[idx], d.rate[idx], d.group[idx], d.keq[idx], d.Et[idx])

# Leave-one-GROUP-out folds (by the full Article|Fig group key), one per unique group in
# first-appearance order. Unlike `_article_folds` this returns EVERY group (incl. a degenerate
# empty-train fold when the corpus has a single group); `_cv_folds` filters those.
function _group_folds(d::Dataset)
    groups  = unique(d.group)
    all_idx = collect(1:length(d.group))
    [ (held_out = g,
       train    = [i for i in all_idx if d.group[i] != g],
       test     = [i for i in all_idx if d.group[i] == g]) for g in groups ]
end

"Non-degenerate (non-empty train AND test) folds of `scheme` ∈ `CV_SCHEMES`."
function _cv_folds(d::Dataset, scheme::Symbol)
    folds = scheme === :article ? _article_folds(d) :
            scheme === :figure  ? _group_folds(d)   :
            error("_cv_folds: unknown CV scheme :$scheme (expected one of $(CV_SCHEMES))")
    [f for f in folds if !isempty(f.train) && !isempty(f.test)]
end

# Summary of one scheme's per-fold results (`per` = NamedTuples carrying `loss`). The NaN /
# 0.0 guards cover a scheme with no non-degenerate fold (absolute mode's :article).
function _cv_summary(per)
    losses = Float64[x.loss for x in per]
    (per_fold = per,
     mean_cv  = isempty(losses) ? NaN : mean(losses),
     se       = isempty(losses) ? 0.0 : std(losses) / sqrt(length(losses)),
     n_folds  = length(losses))
end

# The live drivers are `_cha_loocv` (serial, run.jl) and the pmap task path
# (`_build_tasks` → `_run_fit_task` → `_reduce_cells`, run.jl); both iterate CV_SCHEMES ×
# `_cv_folds` in the same order and reduce through `_cv_summary`, so they agree bit-for-bit.
