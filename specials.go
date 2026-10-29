// Copyright (c) 2005-2026 Fredrik Johansson and mpmath contributors.
// Go port distributed under the same BSD 3-clause licence (see LICENSE).

package mpmath

import "github.com/mpmath/mpmath/internal/libmp"

func altzetaF(x libmp.MPF, p int, r libmp.Rnd) libmp.MPF { return libmp.MpfAltzeta(x, p, r) }
func altzetaC(z libmp.MPC, p int, r libmp.Rnd) libmp.MPC { return libmp.MpcAltzeta(z, p, r) }

func must[T any](v T, err error) T {
	if err != nil {
		if e, ok := err.(*Error); ok {
			libmp.Raise(e.Class, "%s", e.Msg)
		}
		panic(err)
	}
	return v
}

// zeta is zeta(s) with a = 1 and no derivative: the Borwein builtin, then the
// original's shortcuts before its Hurwitz code, which this port does not carry.
func (c *Context) zeta(s any) Number {
	var r Number
	pe := libmp.Try(func() { r = c.b1(s, zetaF, zetaC) })
	if pe == nil {
		return r
	}
	if !libmp.Is(pe, "NotImplementedError") {
		panic(pe)
	}
	n := c.convert(s, true)
	prec := c.prec
	if !c.truthy(n) {
		return c.binop("-", must(c.Mpf(0.5)), 1)
	}
	im, re := c.absN(c.im(n)), c.absN(c.re(n))
	if must(c.Gt(c.absN(im), 500*prec)) && must(c.Lt(c.binop("*", 10, re), prec)) {
		libmp.Raise("NotImplementedError", "zeta: the Riemann-Siegel algorithm is not part of this port")
	}
	if must(c.Eq(n, 1)) {
		return c.mkf(libmp.Finf)
	}
	abss := c.absN(n)
	if must(c.Eq(abss, c.mkf(libmp.Finf))) {
		if must(c.Eq(c.re(n), c.mkf(libmp.Finf))) {
			return c.one()
		}
		return c.binop("*", n, 0)
	} else if must(c.IsNan(abss)) {
		return c.binop("/", 1, n)
	}
	if must(c.Gt(c.re(n), 2*c.prec)) {
		return c.binop("+", c.one(), must(c.Power(2, c.neg(n))))
	}
	libmp.Raise("NotImplementedError", "zeta: the Hurwitz zeta code is not part of this port")
	return nil
}

// Zeta is zeta(s), the Riemann zeta function (a = 1, no derivative).
func (c *Context) Zeta(s any) (Number, error) { return call(func() Number { return c.zeta(s) }) }

// Altzeta is altzeta(s), the Dirichlet eta function.
func (c *Context) Altzeta(s any) (Number, error) {
	var r Number
	pe := libmp.Try(func() { r = c.b1(s, altzetaF, altzetaC) })
	if pe == nil {
		return r, nil
	}
	if !libmp.Is(pe, "NotImplementedError") {
		return nil, &Error{Class: pe.Class, Msg: pe.Msg}
	}
	return c.wrapped(func(a []Number) Number {
		x := a[0]
		if must(c.Eq(x, 1)) {
			return c.binop("+", c.constant("ln2"), c.binop("*", 0, x))
		}
		return c.binop("*", c.neg(must(c.Powm1(2, c.binop("-", 1, x)))), c.zeta(x))
	}, s)
}

// Ellipe is ellipe(m), the complete elliptic integral of the second kind.
func (c *Context) Ellipe(m any) (Number, error) {
	return c.wrapped(func(a []Number) Number {
		return c.builtin(a[0], c.prec, c.rnd, libmp.MpfEllipe, libmp.MpcEllipe)
	}, m)
}

// Ci is ci(z), the cosine integral.
func (c *Context) Ci(z any) (Number, error) { return c.pub(z, libmp.MpfCi, libmp.MpcCi) }

// Si is si(z), the sine integral.
func (c *Context) Si(z any) (Number, error) { return c.pub(z, libmp.MpfSi, libmp.MpcSi) }

func (c *Context) erfLike(z any, ff mpfFn) Number {
	n := c.convert(z, true)
	if x, ok := n.(Mpf); ok {
		return c.builtin(x, c.prec, c.rnd, ff, nil)
	}
	w := n.(Mpc)
	if w.v.Im.Eq(libmp.Fzero) {
		r := c.builtin(c.mkf(w.v.Re), c.prec, c.rnd, ff, nil)
		return c.newMpc(r, 0)
	}
	libmp.Raise("NotImplementedError", "erf and erfc of a complex number are not part of this port")
	return nil
}

// Erf is erf(z) for a real z (or a complex z with zero imaginary part).
func (c *Context) Erf(z any) (Number, error) {
	return call(func() Number { return c.erfLike(z, libmp.MpfErf) })
}

// Erfc is erfc(z) for a real z (or a complex z with zero imaginary part).
func (c *Context) Erfc(z any) (Number, error) {
	return call(func() Number { return c.erfLike(z, libmp.MpfErfc) })
}

func (c *Context) besselj(n int, z any) Number {
	if n < 0 {
		sign := 1
		if n%2 != 0 {
			sign = -1
		}
		return c.binop("*", sign, c.besselj(-n, z))
	}
	w := c.convert(z, true)
	M := c.mag(w)
	if m, ok := M.(int); ok && m > -10 && m < 10 && n < 20 {
		if x, ok := w.(Mpf); ok {
			return c.mkf(libmp.MpfBesseljn(n, x.v, c.prec, c.rnd))
		}
		return c.mkc(libmp.MpcBesseljn(n, w.(Mpc).v, c.prec, c.rnd))
	}
	if !c.truthy(w) {
		if n == 0 {
			return c.binop("+", c.binop("+", c.one(), n), w)
		}
		return c.binop("*", n, w)
	}
	libmp.Raise("NotImplementedError", "besselj outside the fixed-point range is not part of this port")
	return nil
}

// Besselj is besselj(n, z) for an integer order n.
func (c *Context) Besselj(n int, z any) (Number, error) {
	return call(func() Number { return c.besselj(n, z) })
}

// Format is format(x, spec) of an mpf or mpc.
func (c *Context) Format(x Number, spec string) (string, error) {
	return call(func() string {
		switch t := x.(type) {
		case Mpf:
			return libmp.FormatMpf(t.v, spec, c.prec, c.rnd, c.prettyReprDps)
		case Mpc:
			return libmp.FormatMpc(t.v.Re, t.v.Im, spec, c.prec, c.rnd, c.prettyReprDps)
		}
		libmp.Raise("TypeError", "unsupported format string passed to %s.__format__", typeName(x))
		return ""
	})
}
