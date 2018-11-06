using Test
using BARONTools

@testset "sum" begin
    data = parse_sum(joinpath(@__DIR__, "data", "sum.lst"))
    @test data.iteration == [1, 66, 172, 275, 379, 482, 583, 683, 784, 882, 982, 1078, 1141]
    @test data.subproblem_type == fill(BARONTools.RELAXATION, length(data.iteration))
    @test data.open_nodes == [1, 39, 104, 163, 230, 293, 349, 412, 473, 531, 586, 645, 685]
    @test data.time == [1.00, 6.00, 11.00, 16.00, 21.00, 26.00, 32.00, 37.00, 42.00, 47.00, 52.00, 57.00, 60.00]
    @test data.lower_bound == [0.306712, 0.310216, 0.311283, 0.312151, 0.312805, 0.313522, 0.313815, 0.314224, 0.314815, 0.315142, 0.315705, 0.316093, 0.316249]
    @test data.upper_bound == [0.417344, 0.417344, 0.417344, 0.417344, 0.417344, 0.417344, 0.417344, 0.417344, 0.417344, 0.417344, 0.417344, 0.417344, 0.417344]
end
