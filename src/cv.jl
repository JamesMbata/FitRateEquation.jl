# Leave-one-ARTICLE-out CV (goodness-of-fit). Folds by article (the part of the
# Article|Fig group key before '|'), so sibling figures from one paper never leak
# across train/test. Held-out rows are scored with the same per-group centered
# log-ratio loss at the train-fit params; singleton groups carry no centered signal.

_article(group::AbstractString) = String(split(group, '|')[1])

"Train/test row-index folds, one per article."
function _article_folds(d::Dataset)
    arts = _article.(d.group)
    folds = NamedTuple[]
    for a in unique(arts)
        test  = findall(==(a), arts)
        train = findall(!=(a), arts)
        (isempty(test) || isempty(train)) && continue
        push!(folds, (article=a, train=train, test=test))
    end
    folds
end

# Carry Et through the subset so absolute-mode CV folds keep each row's enzyme concentration
# (byte-identical for relative mode: Et is all-NaN there, so d.Et[idx] == fill(NaN,...)).
_subset(d::Dataset, idx) =
    Dataset(d.concs[idx], d.rate[idx], d.group[idx], d.keq[idx], d.Et[idx])

# Leave-one-GROUP-out folds (by the full Article|Fig group key). Absolute mode fits a single
# article, where leave-one-article-out is degenerate; folding on Fig group instead gives a
# genuine held-out CV. One fold per unique group; `article` carries the held-out group key so
# the fold shares `_article_folds`' NamedTuple shape and run.jl's task-build / reduce / _cha_loocv
# treat both uniformly. Unlike `_article_folds` this returns EVERY group (incl. any degenerate
# empty-train singleton) — callers filter degenerate folds at the use site.
function _group_folds(d::Dataset)
    groups  = unique(d.group)
    all_idx = collect(1:length(d.group))
    [ (article = g,
       train   = [i for i in all_idx if d.group[i] != g],
       test    = [i for i in all_idx if d.group[i] == g]) for g in groups ]
end

# Leave-one-article-out CV runs in Cha macro-coordinate space; the live driver is
# `_cha_loocv` (run.jl), which reuses the `_article`/`_article_folds`/`_subset` helpers
# above. The retired coordinate-space `loocv_by_article` twin was removed with the
# coefficient-space path.
