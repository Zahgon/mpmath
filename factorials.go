// Copyright (c) 2005-2026 Fredrik Johansson and mpmath contributors.
// Go port distributed under the same BSD 3-clause licence (see LICENSE).

package mpmath

import "github.com/mpmath/mpmath/internal/libmp"

// Gamma products and factorial-type functions (functions/factorials.py).

func (c *Context) convertAll(xs []any) []Number {
	out := make([]Number, len(xs))
	for i, x := range xs {
		out[i] = c.convert(x, true)
	}
	return out
}

// epsShift is `x and x*(1+eps) or x+eps`, nudging a pole off to the right.
func (c *Context) epsShift(x Number) Number {
	eps := c.constant("eps")
	if c.truthy(x) {
		if v := c.binop("*", x, c.binop("+", 1, eps)); c.truthy(v) {
			return v
		}
	}
	return c.binop("+", x, eps)
}

func (c *Context) gammaprod(a, b []Number, infsign bool) Number {
	var polesNum, polesDen, regNum, regDen []Number
	for _, x := range a {
		if c.isnpint(x) {
			polesNum = append(polesNum, x)
		} else {
			regNum = append(regNum, x)
		}
	}
	for _, x := range b {
		if c.isnpint(x) {
			polesDen = append(polesDen, x)
		} else {
			regDen = append(regDen, x)
		}
	}
	if len(polesNum) < len(polesDen) {
		return c.zero()
	}
	if len(polesNum) > len(polesDen) {
		if infsign {
			var na, nb []Number
			for _, x := range polesNum {
				na = append(na, c.epsShift(x))
			}
			for _, x := range polesDen {
				nb = append(nb, c.epsShift(x))
			}
			g := c.gammaprod(append(na, regNum...), append(nb, regDen...), false)
			s, err := c.Sign(g)
			if err != nil {
				panic(err)
			}
			return c.binop("*", s, c.mkf(libmp.Finf))
		}
		return c.mkf(libmp.Finf)
	}
	var p Number = c.one()
	c.withPrec(c.prec+15, func() {
		for len(polesNum) > 0 {
			i := polesNum[len(polesNum)-1]
			polesNum = polesNum[:len(polesNum)-1]
			j := polesDen[len(polesDen)-1]
			polesDen = polesDen[:len(polesDen)-1]
			t := c.binop("*", c.binop("**", -1, c.binop("+", i, j)), c.gamma(c.binop("-", 1, j)))
			p = c.binop("*", p, c.binop("/", t, c.gamma(c.binop("-", 1, i))))
		}
		for _, x := range regNum {
			p = c.binop("*", p, c.gamma(x))
		}
		for _, x := range regDen {
			p = c.binop("/", p, c.gamma(x))
		}
	})
	return c.pos(p)
}

// Gammaprod is gamma(a[0])*gamma(a[1])*... / (gamma(b[0])*gamma(b[1])*...)
// with poles cancelled in the limit sense.
func (c *Context) Gammaprod(a, b []any) (Number, error) {
	return call(func() Number { return c.gammaprod(c.convertAll(a), c.convertAll(b), false) })
}

// Beta is the beta function B(x, y).
func (c *Context) Beta(x, y any) (Number, error) {
	return call(func() Number {
		xn, yn := c.convert(x, true), c.convert(y, true)
		if c.isinf(yn) {
			xn, yn = yn, xn
		}
		if c.isinf(xn) {
			if c.cmpop("==", xn, c.mkf(libmp.Finf)) && !c.truthy(c.im(yn)) {
				if c.cmpop("==", yn, c.ninf()) {
					return c.mkf(libmp.Fnan)
				}
				if c.cmpop(">", yn, 0) {
					return c.zero()
				}
				if c.isint(yn, false) {
					return c.mkf(libmp.Fnan)
				}
				if c.cmpop("<", yn, 0) {
					s, err := c.Sign(c.gamma(yn))
					if err != nil {
						panic(err)
					}
					return c.binop("*", s, c.mkf(libmp.Finf))
				}
			}
			return c.mkf(libmp.Fnan)
		}
		p := 2 * c.prec
		xy := c.farith("+", xn, yn, PrecOptions{Prec: &p})
		return c.gammaprod([]Number{xn, yn}, []Number{xy}, false)
	})
}

// Binomial is the binomial coefficient (n choose k) for any n, k.
func (c *Context) Binomial(n, k any) (Number, error) {
	return call(func() Number {
		p := 2 * c.prec
		o := PrecOptions{Prec: &p}
		n1 := c.farith("+", n, 1, o)
		k1 := c.farith("+", k, 1, o)
		nk1 := c.farith("-", n1, k, o)
		return c.gammaprod([]Number{n1}, []Number{k1, nk1}, false)
	})
}

// Rf is the rising factorial (x)_n.
func (c *Context) Rf(x, n any) (Number, error) {
	return call(func() Number {
		p := 2 * c.prec
		xn := c.farith("+", x, n, PrecOptions{Prec: &p})
		return c.gammaprod([]Number{xn}, []Number{c.convert(x, true)}, false)
	})
}

// Ff is the falling factorial x(x-1)...(x-n+1).
func (c *Context) Ff(x, n any) (Number, error) {
	return call(func() Number {
		p := 2 * c.prec
		o := PrecOptions{Prec: &p}
		x1 := c.farith("+", x, 1, o)
		xn1 := c.farith("+", c.farith("-", x, n, o), 1, o)
		return c.gammaprod([]Number{x1}, []Number{xn1}, false)
	})
}

// Fac2 is the double factorial x!!.
func (c *Context) Fac2(x any) (Number, error) {
	return c.wrapped(w1(func(x Number) Number {
		if c.isinf(x) {
			if c.cmpop("==", x, c.mkf(libmp.Finf)) {
				return x
			}
			return c.mkf(libmp.Fnan)
		}
		a := c.binop("**", 2, c.binop("/", x, 2))
		e := c.binop("/", c.binop("-", c.b1(x, libmp.MpfCosPi, libmp.MpcCosPi), 1), 4)
		b := c.binop("**", c.binop("/", c.constant("pi"), 2), e)
		g := c.gamma(c.binop("+", c.binop("/", x, 2), 1))
		return c.binop("*", c.binop("*", a, b), g)
	}), x)
}
