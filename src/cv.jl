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

_subset(d::Dataset, idx) =
    Dataset(d.concs[idx], d.rate[idx], d.group[idx], d.keq[idx])

# Leave-one-article-out CV runs in Cha macro-coordinate space; the live driver is
# `_cha_loocv` (run.jl), which reuses the `_article`/`_article_folds`/`_subset` helpers
# above. The retired coordinate-space `loocv_by_article` twin was removed with the
# coefficient-space path.
