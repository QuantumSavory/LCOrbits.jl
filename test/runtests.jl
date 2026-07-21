JET_flag = ARGS == ["jet"]

if JET_flag
  @info "Running JET tests in their dedicated test environment."
  using Pkg
  Pkg.activate(joinpath(@__DIR__, "projects", "jet"))
  Pkg.instantiate()
else
  @info "Skipping JET tests -- pass `test_args=[\"jet\"]` to Pkg.test to enable them."
end

using LCOrbits
using TestItemRunner

# filter for the test
testfilter = ti -> begin
  exclude = Symbol[]
  if !JET_flag
    push!(exclude, :jet)
  else
    return :jet in ti.tags
  end
  if !(VERSION >= v"1.10")
    push!(exclude, :doctests)
    push!(exclude, :aqua)
  end

  return all(!in(exclude), ti.tags)
end

println("Starting tests with $(Threads.nthreads()) threads out of `Sys.CPU_THREADS = $(Sys.CPU_THREADS)`...")

@run_package_tests filter=testfilter
