@testitem "JET analysis" tags=[:jet] begin

using JET
using Test
using LCOrbits
using CSV, DataFrames

rep = report_package("LCOrbits";
    ignored_modules=(
        LastFrameModule(Base),
        AnyFrameModule(CSV),
        LastFrameModule(DataFrames),
    )
)
@show rep
@test_broken length(JET.get_reports(rep)) == 0
@test length(JET.get_reports(rep)) <= 1

end
