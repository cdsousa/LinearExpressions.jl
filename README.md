LinearExpressions
=================

**Linear symbolic expressions for the Julia language**


[![Build Status](https://travis-ci.org/cdsousa/LinearExpressions.jl.png)](https://travis-ci.org/cdsousa/LinearExpressions.jl)

Examples
--------

```Julia
julia> using LinearExpressions

julia> a, b, c, d, w, x, y, z = map(RealVariable, [:a, :b, :c, :d, :w, :x, :y, :z])
8-element Array{TypedVariable{Real},1}:
 a
 b
 c
 d
 w
 x
 y
 z

julia> 2x
2x

julia> 2x + 1
1 + 2x

julia> 4 * (3 + 2x)
12 + 8x

julia> x - x == 0
true

julia> le = x + y + 2.3 * (w + z) - 3.4 * (a + b) - c - d + 1.23
1.23 - d - 3.4b - 3.4a + y + 2.3w + 2.3z - c + x

julia> le + 3*le
4.92 - 4.0d - 13.6b - 13.6a + 4.0y + 9.2w + 9.2z - 4.0c + 4.0x

julia> [x y; y 1] + [1.5 2; 3 4] * [x 0; z 0]
2x2 Array{LinExpr{Float64,TypedVariable{Real}},2}:
 2.0z+2.5x    y
 y+4.0z+3.0x  1.0

julia> [1 2; 2 1] + [0.0 1.2; 1.2 0.0]*x + 4*eye(2)*y
1	2
2	1

+

4	0
0	4
y

+

0	1.2
1.2	0
x
```


Author
------

[Cristóvão Duarte Sousa](https://github.com/cdsousa)

Install
-------

Within Julia:

    Pkg.clone("git://github.com/cdsousa/LinearExpressions.jl.git")

License
-------

MIT "Expat" License. See [License File](LICENSE.md)
