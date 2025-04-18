using LCOrbits
using TestItemRunner

# filter for the test
testfilter = ti -> begin
  exclude = Symbol[]
  if get(ENV,"JET_TEST","")!="true"
    push!(exclude, :jet)
  end
  if !(VERSION >= v"1.10")
    push!(exclude, :doctests)
    push!(exclude, :aqua)
  end

  return all(!in(exclude), ti.tags)
end

println("Starting tests with $(Threads.nthreads()) threads out of `Sys.CPU_THREADS = $(Sys.CPU_THREADS)`...")

@run_package_tests filter=testfilter
