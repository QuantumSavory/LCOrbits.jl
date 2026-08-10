@testitem "Aqua analysis" tags=[:aqua] begin

using Aqua, LCOrbits

Aqua.test_all(LCOrbits)

end
