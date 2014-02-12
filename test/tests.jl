println("--------- statring LinearExpressions tests ----------")

using Base.Test

using LinearExpressions
@test isdefined(:LinearExpressions)
@test typeof(LinearExpressions) === Module


a, b, c, d, w, x, y, z = map(RealVariable, [:a, :b, :c, :d, :w, :x, :y, :z])
ix = TypedVariable{Integer}(:x)

@test x == x
@test x != y
@test x != ix
@test 1x == x
@test 1x == 1.0x
@test -x == -1x
@test x+2x == 3x
@test x-x == 0
@test 1+x-x == 1
@test 2*(x-x) == 0
@test x+y-x == y
@test 2.2x+4.4y == 2*(1.1x+2.2y)

@test_throws (1+2x) * (3+4y) # no polynomials

@test conj(x) == x
@test conj(2x) == 2x


A = rand(Int, 2, 2)
B = rand(Int, 2, 2)

@test A + A*x + B + B*x == (A+B) + (A+B)*x


le1 = y + 2.3 * (w + z) - 3.4 * (a + b) - c - d + 1.23
le2 = 1.23 - d - 3.4b - 3.4a + y + 2.3w + 2.3z - c

@test le1 == le2
@test differentiate(le1, x) == 0
@test differentiate(le1, w) == 2.3


function test_eq()
  a1 = LinExpr(-1.2, [x=>-1., y=>2.])
  a2 = LinExpr(-1.2, [x=>-1., y=>2.])
  b = LinExpr(-3.4, [x=>-3., y=>4.])
  @test a1 == a2
  @test a1 != b
end

test_eq()


function test_sum()
  a = LinExpr(-1.2, [x=>-1., y=>2.])
  b = LinExpr(-3.4, [x=>-3., y=>4.])
  c = LinExpr(-3.4, [x=>1., z=>4.])
  r_a_b = LinExpr(-4.6, [x=>-4., y=>6.])
  r_a_c = LinExpr(-4.6, [y=>2., z=>4.])
  @test a+b == r_a_b
  @test a+c == r_a_c
end

test_sum()


@test convert(RealVariable, LinExpr{Float64, RealVariable}([x=>1.0])) == x

@test_throws convert(RealVariable, LinExpr(0.0, [x=>1.0, y=>1.0]))
@test_throws convert(RealVariable, LinExpr{Float64, RealVariable}(1.23))

@test convert(Float64, LinExpr{Float64, RealVariable}(1.23)) == 1.23
@test_throws convert(Float64, LinExpr(1.23, [x=>4.56]))


println("✓ all LinearExpressions tests passed")

