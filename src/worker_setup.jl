# =========================================================================================
#              worker_setup.jl — shared local/Savio worker dispatch setup
# =========================================================================================
#
# `setup_workers(nprocs)` centralizes the env-detect + addprocs block used by every
# FitRateEquation runner (run_g6pd/run_pgd/run_hk1). Library-safe: nothing here runs at
# include (module-load) time — a library must not mutate the active project or spawn
# processes just because it was `using`'d. Everything is wrapped in `setup_workers`, called
# explicitly by the runners. Workers load the module with `@everywhere using
# FitRateEquation` (they already have it as an installed dependency via the active
# project), not `@everywhere include(...)` — a library doesn't know its own source path
# the way a standalone script did.
#
#   * Savio (inside an sbatch alloc): SLURM_JOB_CPUS_PER_NODE present -> SlurmManager
#     sized to the allocation (unaffected by `nprocs`).
#   * Local PC: addprocs(n) workers; n defaults to the guarded value below when `nprocs`
#     is nothing.
#
# GUARDED DEFAULT: max(1, min(3, Sys.CPU_THREADS - 1)). Capped at 3 to avoid the RAM
# oversubscription that OOM-crashes full-budget fits — each worker holds its own copy
# of the model/data, so memory (not CPU) is the binding constraint. Floored at 1 so a
# 1-2 logical-core machine yields n == 1, which the `n > 1` guard below turns into
# serial master-only execution (never addprocs(0)). Override with FRE_NPROCS.

"""
    setup_workers(nprocs=nothing) -> Vector{Int}

Ensure worker processes are available for a FitRateEquation run and load the module
(`FitRateEquation`, `LinearAlgebra` with single-threaded BLAS) on each of them. Returns
`Distributed.workers()`. Idempotent: if worker processes already exist, no new ones are
spawned and the module load is skipped.

  * Inside a Savio/SLURM allocation (`SLURM_JOB_CPUS_PER_NODE` set), spawns a
    `SlurmManager` sized to the allocation; `nprocs` is ignored in that branch.
  * Otherwise (local), spawns `addprocs(n)` where `n = nprocs` if given, else the
    `FRE_NPROCS` env var if set, else `max(1, min(3, Sys.CPU_THREADS - 1))`. `n <= 1`
    runs serially on the master (no `addprocs` call).
"""
function setup_workers(nprocs=nothing)
    if nworkers() == 1 && workers() == [1]   # no worker processes spawned yet
        if haskey(ENV, "SLURM_JOB_CPUS_PER_NODE")
            # Savio: submitter set --project; only instantiate (idempotent). No resolve
            # against stale cluster registries. Called here (at setup_workers time), not
            # at include time, so `using FitRateEquation` never touches the active project.
            Pkg.instantiate()
            ENV["JULIA_WORKER_TIMEOUT"] = "600"
            let subs = Dict("x" => "*", "(" => "", ")" => "")
                np = sum(eval(Meta.parse(replace(
                    ENV["SLURM_JOB_CPUS_PER_NODE"], r"x|\(|\)" => s -> subs[s],
                ))))
                addprocs(SlurmManager(np); exeflags = "--project=$(Base.active_project())")
            end
        else
            n = nprocs !== nothing ? nprocs :
                parse(Int, get(ENV, "FRE_NPROCS", string(max(1, min(3, Sys.CPU_THREADS - 1)))))
            n > 1 && addprocs(n; exeflags = "--project=$(Base.active_project())")
        end
    end

    new_workers = filter(!=(1), workers())
    if !isempty(new_workers)
        # `@everywhere` macroexpands to a `:toplevel` Expr, which Julia only allows at file
        # top level — not inside a function body (this function). Call the function it
        # forwards to directly instead; same effect (`using` runs under `Main` on each pid).
        Distributed.remotecall_eval(Main, new_workers, quote
            using FitRateEquation
            using LinearAlgebra
            LinearAlgebra.BLAS.set_num_threads(1)   # avoid BLAS/CMA-ES oversubscription
        end)
    end
    LinearAlgebra.BLAS.set_num_threads(1)   # master too
    workers()
end
