using ITensors
using ITensorUnicodePlots
using ReferenceTests
using Test

@testset "Basic test for UnicodePlots backend" begin
  extension = "txt"

  N = 10
  s(n) = Index([QN("Sz", 0) => 1, QN("Sz", 1) => 1]; tags="S=1/2,Site,n=$n")
  l(n) = Index([QN("Sz", 0) => 10, QN("Sz", 1) => 10]; tags="Link,l=$n")
  h(n) = Index([QN("Sz", 0) => 5, QN("Sz", 1) => 5]; tags="ham,Link,l=$n")
  s⃗ = [s(n) for n in 1:N]
  l⃗ = [l(n) for n in 1:(N - 1)]
  h⃗ = [h(n) for n in 1:(N - 1)]

  # Add some more indices between two of the tensors
  x = Index([QN("Sz", 0) => 2]; tags="X")
  y = Index([QN("Sz", 0) => 2]; tags="Y")

  n = 2
  ψn1n2 = randomITensor(l⃗[n - 1], s⃗[n], s⃗[n + 1], l⃗[n + 1], dag(x), dag(y))
  hn1 = randomITensor(dag(h⃗[n - 1]), s⃗[n]', dag(s⃗[n]), h⃗[n], x, y)
  hn2 = randomITensor(dag(h⃗[n]), s⃗[n + 1]', dag(s⃗[n + 1]), h⃗[n + 1])
  ELn0 = randomITensor(l⃗[n - 1]', h⃗[n - 1], dag(l⃗[n - 1]))
  ERn2 = randomITensor(l⃗[n + 1]', dag(h⃗[n + 1]), dag(l⃗[n + 1]))

  tn = [ELn0, ψn1n2, hn1, hn2, ERn2]

  R = @visualize ELn0 * ψn1n2 * hn1 * hn2 * ERn2
  R1 = @visualize ELn0 * ψn1n2 * hn1
  R2 = @visualize R1 * hn2 * ERn2 vertex_labels = ["T1", "T2", "T3"]
  tn2 = @visualize tn

  @test R ≈ ELn0 * ψn1n2 * hn1 * hn2 * ERn2
  @test R1 ≈ ELn0 * ψn1n2 * hn1
  @test R2 ≈ ELn0 * ψn1n2 * hn1 * hn2 * ERn2
  @test all(tn .== tn2)

  R = @visualize figR ELn0 * ψn1n2 * hn1 * hn2 * ERn2
  R_tags = @visualize figR_tags ELn0 * ψn1n2 * hn1 * hn2 * ERn2 edge_labels = (tags=true,)
  R1 = @visualize figR1 ELn0 * ψn1n2 * hn1
  R2 = @visualize figR2 R1 * hn2 * ERn2 vertex_labels = ["T1", "T2", "T3"]

  fig_tn = @visualize_noeval tn

  by = isequal

  @test_reference "references/R.$extension" figR by = by
  @test_reference "references/R_tags.$extension" figR_tags by = by
  @test_reference "references/R1.$extension" figR1 by = by
  @test_reference "references/R2.$extension" figR2 by = by
  @test_reference "references/tn.$extension" fig_tn by = by

  @test_throws DimensionMismatch @visualize fig R1 * hn2 * ERn2 vertex_labels = ["T1", "T2"]
end
