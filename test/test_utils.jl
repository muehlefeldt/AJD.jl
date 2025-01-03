using AJD


@testset "Generate commuting matrices" begin
    matrices = AJD.random_normal_commuting_matrices(4, 10)
    for index in 1:length(matrices)-1
        @test isapprox(matrices[index] * matrices[index+1], matrices[index+1] * matrices[index])
    end

    matrices = AJD.random_normal_commuting_matrices(4, 10; complex=true)
    for index in 1:length(matrices)-1
        @test isapprox(matrices[index] * matrices[index+1], matrices[index+1] * matrices[index])
    end
end