using RigorousInvariantMeasures, Test

@testset "Progress bars" begin
    RIM = RigorousInvariantMeasures
    for (v, expected) in (("0", false), ("false", false), ("off", false), ("no", false),
                          ("1", true), ("true", true), ("ON", true), (" yes ", true))
        @test withenv(() -> RIM._progress_bars_default(), "PROGRESS_BARS" => v) == expected
    end
    @test withenv(() -> RIM._progress_bars_default(), "PROGRESS_BARS" => nothing) ==
          isinteractive()
    @test withenv(() -> RIM._progress_bars_default(), "PROGRESS_BARS" => "maybe") ==
          isinteractive()
    old = RIM.SHOW_PROGRESS_BARS[]
    @test show_progress_bars!(true) === true && RIM.SHOW_PROGRESS_BARS[]
    @test show_progress_bars!(false) === false && !RIM.SHOW_PROGRESS_BARS[]
    show_progress_bars!(old)
end
