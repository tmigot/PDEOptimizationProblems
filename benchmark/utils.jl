using DataFrames, Dates, JLD2, SolverTools, SolverBenchmark
using Plots; # pgfplots()

function make_figures(name, file_prefix::String = "results")
  legend = Dict(
    :neval_obj => "number of f evals",
    :neval_cons => "number of c evals",
    :neval_grad => "number of ∇f evals",
    :neval_jac => "number of ∇c evals",
    :neval_jprod => "number of ∇c*v evals",
    :neval_jtprod => "number of ∇cᵀ*v evals",
    :neval_hess => "number of ∇²f evals",
    :elapsed_time => "elapsed time",
  )
  styles = [:solid, :dash, :dot, :dashdot] #[:auto, :solid, :dash, :dot, :dashdot, :dashdotdot]
  perf_title(col) = "Performance profile on CUTEst w.r.t. $(string(legend[col]))"

  @load string(name, ".jld2") stats

  for col in keys(legend)
    empty = false
    for df in values(stats)
      if all(df[!, col] .== 0)
        empty = true
      end
    end

    if !empty
      ϵ = minimum(minimum(filter(x -> x > 0, df[!, col])) for df in values(stats))
      first_order(df) = df.status .== :first_order
      unbounded(df) = df.status .== :unbounded
      solved(df) = first_order(df) .| unbounded(df)
      cost(df) = (max.(df[!, col], ϵ) + .!solved(df) .* Inf)
      #p = performance_profile(stats, cost)
      p = performance_profile(
        stats,
        cost,
        title = perf_title(col),
        legend = :bottomright,
        linestyles = styles,
      )
      # savefig("$(file_prefix)_perf-$col.tex")
      png("$(file_prefix)_perf-$col")
      #profile_solvers(stats, [cost], ["$(col)"])
      costs = [cost]
      solvers = collect(keys(stats))
      nsolvers = length(solvers)
      npairs = div(nsolvers * (nsolvers - 1), 2)
      colors = get_color_palette(:auto, nsolvers)
      if nsolvers > 2
        ipairs = 0
        # combinations of solvers 2 by 2
        for i = 2:nsolvers
          for j = 1:(i - 1)
            ipairs += 1
            pair = [solvers[i], solvers[j]]
            dfs = (stats[solver] for solver in pair)
            Ps = [hcat([cost(df) for df in dfs]...) for cost in costs]

            clrs = [colors[i], colors[j]]
            stls = [styles[i], styles[j]]
            p = performance_profile(
              Ps[1],
              string.(pair),
              palette = clrs,
              legend = :bottomright,
              styles = stls,
            )
            ipairs < npairs && xlabel!(p, "")
            # savefig("$(tod)_$(solvers[i])_$(solvers[j])_perf-$col.tex")
            png("$(file_prefix)_$(solvers[i])_$(solvers[j])_perf-$col")
          end
        end
      else
      end
    end
  end
end

using JLD2, DataFrames, Dates, SolverCore, NLPModels, SolverBenchmark

function make_md(name, opt_val, file_prefix::String = "results")
  @load string(name, ".jld2") stats

  open("$(file_prefix).md", "w") do io
    for solver in collect(keys(stats))
      println(io, "## $solver")
      df = stats[solver]
      insertcols!(df, 1, :opt_val => opt_val)
      pretty_stats(
        io,
        df[
          !,
          [
            :name,
            :nvar,
            :ncon,
            :status,
            :objective,
            :opt_val,
            :elapsed_time,
            :neval_obj,
            :neval_cons,
            :dual_feas,
            :primal_feas,
          ],
        ],
        tf = tf_markdown,
        hdr_override = Dict(
          :neval_obj => "#f",
          :neval_cons => "#c",
          :elapsed_time => "time (s)",
          :primal_feas => "feas",
        ),
      )
    end
  end
end
