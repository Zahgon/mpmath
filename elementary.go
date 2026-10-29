// Copyright (c) 2005-2026 Fredrik Johansson and mpmath contributors.
// Go port distributed under the same BSD 3-clause licence (see LICENSE).

package mpmath

import (
	"math/big"
	"strings"

	"github.com/mpmath/mpmath/internal/libmp"
)

// Elementary functions: the libmp builtins of ctx_mp.py, the constants, and
// the high-level functions of functions/functions.py.

type mpfFn func(libmp.MPF, int, libmp.Rnd) libmp.MPF
type mpcFn func(libmp.MPC, int, libmp.Rnd) libmp.MPC

// withPrec runs f at precision p. Setting the precision back afterwards also
// recomputes dps from it, exactly as `ctx.prec = prec` does in the original.
func (c *Context) withPrec(p int, f func()) {
	saved := c.prec
	c.SetPrec(p)
	defer c.SetPrec(saved)
	f()
}

func (c *Context) pos(x Number) Number {
	if t, ok := x.(Mpc); ok {
		return c.mkc(libmp.MpcPos(t.v, c.prec, c.rnd))
	}
	return c.mkf(libmp.MpfPos(x.(Mpf).v, c.prec, c.rnd))
}

func (c *Context) neg(x Number) Number {
	if t, ok := x.(Mpc); ok {
		return c.mkc(libmp.MpcNeg(t.v, c.prec, c.rnd))
	}
	return c.mkf(libmp.MpfNeg(x.(Mpf).v, c.prec, c.rnd))
}

func (c *Context) one() Mpf  { return c.mkf(libmp.Fone) }
func (c *Context) zero() Mpf { return c.mkf(libmp.Fzero) }
func (c *Context) jj() Mpc   { return c.mkc(libmp.MPC{Re: libmp.Fzero, Im: libmp.Fone}) }

func (c *Context) re(x Number) Number {
	if t, ok := x.(Mpc); ok {
		return t.Real()
	}
	return x
}

func (c *Context) im(x Number) Number {
	if t, ok := x.(Mpc); ok {
		return t.Imag()
	}
	return c.zero()
}

// wrapped is _wrap_specfun: the arguments are converted at the caller's
// precision, the body runs 10 bits higher, and the result is rounded back.
func (c *Context) wrapped(f func(a []Number) Number, args ...any) (Number, error) {
	return call(func() Number {
		conv := make([]Number, len(args))
		for i, a := range args {
			conv[i] = c.convert(a, true)
		}
		var r Number
		c.withPrec(c.prec+10, func() { r = f(conv) })
		return c.pos(r)
	})
}

// builtin is _wrap_libmp_function at an explicit precision and rounding.
func (c *Context) builtin(x any, prec int, rnd libmp.Rnd, ff mpfFn, cf mpcFn) Number {
	var n Number
	switch t := x.(type) {
	case Mpf:
		n = t
	case Mpc:
		n = t
	default:
		n = c.convert(x, true)
	}
	if t, ok := n.(Mpf); ok {
		var r libmp.MPF
		pe := libmp.Try(func() { r = ff(t.v, prec, rnd) })
		if pe == nil {
			return c.mkf(r)
		}
		if !libmp.Is(pe, "ComplexResult") || c.trapComplex {
			panic(pe)
		}
		if cf == nil {
			libmp.Raise("TypeError", "'NoneType' object is not callable")
		}
		return c.mkc(cf(libmp.MPC{Re: t.v, Im: libmp.Fzero}, prec, rnd))
	}
	if cf == nil {
		libmp.Raise("TypeError", "'NoneType' object is not callable")
	}
	return c.mkc(cf(n.(Mpc).v, prec, rnd))
}

func (c *Context) b1(x any, ff mpfFn, cf mpcFn) Number { return c.builtin(x, c.prec, c.rnd, ff, cf) }

func (c *Context) pub(x any, ff mpfFn, cf mpcFn) (Number, error) {
	return call(func() Number { return c.b1(x, ff, cf) })
}

func zetaF(x libmp.MPF, p int, r libmp.Rnd) libmp.MPF { return libmp.MpfZeta(x, p, r, false) }
func zetaC(z libmp.MPC, p int, r libmp.Rnd) libmp.MPC { return libmp.MpcZeta(z, p, r, false, false) }

func (c *Context) sqrt(x any) Number  { return c.b1(x, libmp.MpfSqrt, libmp.MpcSqrt) }
func (c *Context) exp(x any) Number   { return c.b1(x, libmp.MpfExp, libmp.MpcExp) }
func (c *Context) ln(x any) Number    { return c.b1(x, libmp.MpfLn, libmp.MpcLn) }
func (c *Context) sin(x any) Number   { return c.b1(x, libmp.MpfSin, libmp.MpcSin) }
func (c *Context) cos(x any) Number   { return c.b1(x, libmp.MpfCos, libmp.MpcCos) }
func (c *Context) tan(x any) Number   { return c.b1(x, libmp.MpfTan, libmp.MpcTan) }
func (c *Context) sinh(x any) Number  { return c.b1(x, libmp.MpfSinh, libmp.MpcSinh) }
func (c *Context) cosh(x any) Number  { return c.b1(x, libmp.MpfCosh, libmp.MpcCosh) }
func (c *Context) tanh(x any) Number  { return c.b1(x, libmp.MpfTanh, libmp.MpcTanh) }
func (c *Context) asin(x any) Number  { return c.b1(x, libmp.MpfAsin, libmp.MpcAsin) }
func (c *Context) acos(x any) Number  { return c.b1(x, libmp.MpfAcos, libmp.MpcAcos) }
func (c *Context) atan(x any) Number  { return c.b1(x, libmp.MpfAtan, libmp.MpcAtan) }
func (c *Context) asinh(x any) Number { return c.b1(x, libmp.MpfAsinh, libmp.MpcAsinh) }
func (c *Context) acosh(x any) Number { return c.b1(x, libmp.MpfAcosh, libmp.MpcAcosh) }
func (c *Context) atanh(x any) Number { return c.b1(x, libmp.MpfAtanh, libmp.MpcAtanh) }
func (c *Context) sinpi(x any) Number { return c.b1(x, libmp.MpfSinPi, libmp.MpcSinPi) }
func (c *Context) gamma(x any) Number { return c.b1(x, libmp.MpfGamma, libmp.MpcGamma) }
func (c *Context) expjpi(x any) Number {
	return c.b1(x, libmp.MpfExpjpi, libmp.MpcExpjpi)
}

// The libmp builtins.
func (c *Context) Sqrt(x any) (Number, error) { return c.pub(x, libmp.MpfSqrt, libmp.MpcSqrt) }
func (c *Context) Cbrt(x any) (Number, error) { return c.pub(x, libmp.MpfCbrt, libmp.MpcCbrt) }
func (c *Context) Ln(x any) (Number, error)   { return c.pub(x, libmp.MpfLn, libmp.MpcLn) }
func (c *Context) Exp(x any) (Number, error)  { return c.pub(x, libmp.MpfExp, libmp.MpcExp) }
func (c *Context) Expj(x any) (Number, error) { return c.pub(x, libmp.MpfExpj, libmp.MpcExpj) }
func (c *Context) Expjpi(x any) (Number, error) {
	return c.pub(x, libmp.MpfExpjpi, libmp.MpcExpjpi)
}
func (c *Context) Sin(x any) (Number, error)   { return c.pub(x, libmp.MpfSin, libmp.MpcSin) }
func (c *Context) Cos(x any) (Number, error)   { return c.pub(x, libmp.MpfCos, libmp.MpcCos) }
func (c *Context) Tan(x any) (Number, error)   { return c.pub(x, libmp.MpfTan, libmp.MpcTan) }
func (c *Context) Sinh(x any) (Number, error)  { return c.pub(x, libmp.MpfSinh, libmp.MpcSinh) }
func (c *Context) Cosh(x any) (Number, error)  { return c.pub(x, libmp.MpfCosh, libmp.MpcCosh) }
func (c *Context) Tanh(x any) (Number, error)  { return c.pub(x, libmp.MpfTanh, libmp.MpcTanh) }
func (c *Context) Asin(x any) (Number, error)  { return c.pub(x, libmp.MpfAsin, libmp.MpcAsin) }
func (c *Context) Acos(x any) (Number, error)  { return c.pub(x, libmp.MpfAcos, libmp.MpcAcos) }
func (c *Context) Atan(x any) (Number, error)  { return c.pub(x, libmp.MpfAtan, libmp.MpcAtan) }
func (c *Context) Asinh(x any) (Number, error) { return c.pub(x, libmp.MpfAsinh, libmp.MpcAsinh) }
func (c *Context) Acosh(x any) (Number, error) { return c.pub(x, libmp.MpfAcosh, libmp.MpcAcosh) }
func (c *Context) Atanh(x any) (Number, error) { return c.pub(x, libmp.MpfAtanh, libmp.MpcAtanh) }
func (c *Context) Sinpi(x any) (Number, error) { return c.pub(x, libmp.MpfSinPi, libmp.MpcSinPi) }
func (c *Context) Cospi(x any) (Number, error) { return c.pub(x, libmp.MpfCosPi, libmp.MpcCosPi) }
func (c *Context) Floor(x any) (Number, error) { return c.pub(x, libmp.MpfFloor, libmp.MpcFloor) }
func (c *Context) Ceil(x any) (Number, error)  { return c.pub(x, libmp.MpfCeil, libmp.MpcCeil) }
func (c *Context) Nint(x any) (Number, error)  { return c.pub(x, libmp.MpfNint, libmp.MpcNint) }
func (c *Context) Frac(x any) (Number, error)  { return c.pub(x, libmp.MpfFrac, libmp.MpcFrac) }
func (c *Context) Fib(x any) (Number, error) {
	return c.pub(x, libmp.MpfFibonacci, libmp.MpcFibonacci)
}
func (c *Context) Gamma(x any) (Number, error)  { return c.pub(x, libmp.MpfGamma, libmp.MpcGamma) }
func (c *Context) Rgamma(x any) (Number, error) { return c.pub(x, libmp.MpfRgamma, libmp.MpcRgamma) }
func (c *Context) Loggamma(x any) (Number, error) {
	return c.pub(x, libmp.MpfLoggamma, libmp.MpcLoggamma)
}
func (c *Context) Factorial(x any) (Number, error) {
	return c.pub(x, libmp.MpfFactorial, libmp.MpcFactorial)
}
func (c *Context) Digamma(x any) (Number, error) { return c.pub(x, libmp.MpfPsi0, libmp.MpcPsi0) }
func (c *Context) Harmonic(x any) (Number, error) {
	return c.pub(x, libmp.MpfHarmonic, libmp.MpcHarmonic)
}
func (c *Context) Ei(x any) (Number, error)     { return c.pub(x, libmp.MpfEi, libmp.MpcEi) }
func (c *Context) E1(x any) (Number, error)     { return c.pub(x, libmp.MpfE1, libmp.MpcE1) }
func (c *Context) Ellipk(x any) (Number, error) { return c.pub(x, libmp.MpfEllipk, libmp.MpcEllipk) }
func (c *Context) Agm1(x any) (Number, error)   { return c.pub(x, libmp.MpfAgm1, libmp.MpcAgm1) }

// Psi is the polygamma function psi(m, z).
func (c *Context) Psi(m int, z any) (Number, error) {
	return call(func() Number {
		n := c.convert(z, true)
		if t, ok := n.(Mpf); ok {
			return c.mkf(libmp.MpfPsi(m, t.v, c.prec, c.rnd))
		}
		return c.mkc(libmp.MpcPsi(m, n.(Mpc).v, c.prec, c.rnd))
	})
}

// Bernoulli is the n-th Bernoulli number (B_1 = +1/2 when plus is set).
func (c *Context) Bernoulli(n int, plus bool) (Mpf, error) {
	return call(func() Mpf { return c.mkf(libmp.MpfBernoulli(n, c.prec, c.rnd, plus)) })
}

// Bernfrac is B_n as an exact fraction p/q.
func (c *Context) Bernfrac(n int, plus bool) (*big.Int, *big.Int, error) {
	type pq struct{ p, q *big.Int }
	r, err := call(func() pq { p, q := libmp.Bernfrac(n, plus); return pq{p, q} })
	return r.p, r.q, err
}

var constants = map[string]func(int, libmp.Rnd) libmp.MPF{
	"pi": libmp.MpfPi, "e": libmp.MpfE, "ln2": libmp.MpfLn2, "ln10": libmp.MpfLn10,
	"phi": libmp.MpfPhi, "euler": libmp.MpfEuler, "catalan": libmp.MpfCatalan,
	"khinchin": libmp.MpfKhinchin, "glaisher": libmp.MpfGlaisher, "apery": libmp.MpfApery,
	"degree": libmp.MpfDegree, "twinprime": libmp.MpfTwinprime, "mertens": libmp.MpfMertens,
	"eps": func(prec int, _ libmp.Rnd) libmp.MPF {
		return libmp.MPF{Sign: 0, Man: big.NewInt(1), Exp: 1 - prec, BC: 1}
	},
}

func (c *Context) constant(name string) Mpf { return c.mkf(constants[name](c.prec, c.rnd)) }

// Constant is +mp.<name> for pi, e, ln2, ln10, phi, euler, catalan, khinchin,
// glaisher, apery, degree, twinprime, mertens and eps: the constant at the
// working precision and rounding. o.Prec/o.Dps/o.Rounding stand for the
// keywords of calling the constant, as in mp.pi(dps=50).
func (c *Context) Constant(name string, o PrecOptions) (Mpf, error) {
	return call(func() Mpf {
		f, ok := constants[name]
		if !ok {
			libmp.Raise("AttributeError", "'MPContext' object has no attribute %s", libmp.PyRepr(name))
		}
		prec, rnd := c.prec, c.rnd
		if o.Prec != nil && *o.Prec != 0 {
			prec = *o.Prec
		}
		if o.Rounding != nil && *o.Rounding != "" {
			rnd = libmp.Rnd((*o.Rounding)[0])
		}
		if o.Dps != nil && *o.Dps != 0 {
			prec = libmp.DpsToPrec(*o.Dps)
		}
		return c.mkf(f(prec, rnd))
	})
}

// iop is Python arithmetic on magnitudes (an int, or an mpf for the special
// magnitudes).
func (c *Context) iop(op string, a, b any) any {
	if x, ok := a.(int); ok {
		if y, ok := b.(int); ok {
			switch op {
			case "+":
				return x + y
			case "-":
				return x - y
			}
		}
	}
	return c.binop(op, a, b)
}

func (c *Context) icmp(op string, a, b any) bool {
	if x, ok := a.(int); ok {
		if y, ok := b.(int); ok {
			switch op {
			case "<":
				return x < y
			case "<=":
				return x <= y
			case ">":
				return x > y
			case ">=":
				return x >= y
			case "==":
				return x == y
			case "!=":
				return x != y
			}
		}
	}
	return c.cmpop(op, a, b)
}

// sumAccurately is ctx.sum_accurately.
func (c *Context) sumAccurately(terms func() []any, checkStep int) Number {
	prec := c.prec
	defer c.SetPrec(prec)
	extraprec := 10
	for {
		c.SetPrec(prec + extraprec + 5)
		var maxMag any = c.ninf()
		var sumMag any = c.zero()
		var s Number = c.zero()
		for k, term := range terms() {
			s = c.binop("+", s, term)
			if k%checkStep == 0 && c.truthy(term) {
				termMag := c.mag(term)
				maxMag = c.pyMax(maxMag, termMag)
				sumMag = c.mag(s)
				if c.icmp(">", c.iop("-", sumMag, termMag), c.prec) {
					break
				}
			}
		}
		cancellation := c.iop("-", maxMag, sumMag)
		if c.icmp("!=", cancellation, cancellation) {
			return s
		}
		if c.icmp("<", cancellation, extraprec) {
			return s
		}
		var m any = c.prec
		if c.icmp("<", cancellation, c.prec) {
			m = cancellation
		}
		extraprec = c.iop("+", extraprec, m).(int)
	}
}

func (c *Context) recip(x Number) Number { return c.binop("/", c.one(), x) }

func (c *Context) invOrZero(z Number) Number {
	if c.isinf(z) {
		return c.zero()
	}
	return c.recip(z)
}

func w1(f func(z Number) Number) func([]Number) Number {
	return func(a []Number) Number { return f(a[0]) }
}

// Cot and the other reciprocal trigonometric and hyperbolic functions.
func (c *Context) Cot(z any) (Number, error) {
	return c.wrapped(w1(func(z Number) Number { return c.recip(c.tan(z)) }), z)
}
func (c *Context) Sec(z any) (Number, error) {
	return c.wrapped(w1(func(z Number) Number { return c.recip(c.cos(z)) }), z)
}
func (c *Context) Csc(z any) (Number, error) {
	return c.wrapped(w1(func(z Number) Number { return c.recip(c.sin(z)) }), z)
}
func (c *Context) Coth(z any) (Number, error) {
	return c.wrapped(w1(func(z Number) Number { return c.recip(c.tanh(z)) }), z)
}
func (c *Context) Sech(z any) (Number, error) {
	return c.wrapped(w1(func(z Number) Number { return c.recip(c.cosh(z)) }), z)
}
func (c *Context) Csch(z any) (Number, error) {
	return c.wrapped(w1(func(z Number) Number { return c.recip(c.sinh(z)) }), z)
}
func (c *Context) Acot(z any) (Number, error) {
	return c.wrapped(w1(func(z Number) Number {
		if !c.truthy(z) {
			return c.binop("*", c.constant("pi"), 0.5)
		}
		return c.atan(c.invOrZero(z))
	}), z)
}
func (c *Context) Asec(z any) (Number, error) {
	return c.wrapped(w1(func(z Number) Number { return c.acos(c.invOrZero(z)) }), z)
}
func (c *Context) Acsc(z any) (Number, error) {
	return c.wrapped(w1(func(z Number) Number { return c.asin(c.invOrZero(z)) }), z)
}
func (c *Context) Acoth(z any) (Number, error) {
	return c.wrapped(w1(func(z Number) Number {
		if !c.truthy(z) {
			return c.binop("*", c.constant("pi"), complex(0, 0.5))
		}
		return c.atanh(c.invOrZero(z))
	}), z)
}
func (c *Context) Asech(z any) (Number, error) {
	return c.wrapped(w1(func(z Number) Number { return c.acosh(c.recip(z)) }), z)
}
func (c *Context) Acsch(z any) (Number, error) {
	return c.wrapped(w1(func(z Number) Number { return c.asinh(c.recip(z)) }), z)
}

// Sign is sign(x): -1, 0 or 1 for a real x, x/|x| for a complex one.
func (c *Context) Sign(x any) (Number, error) {
	return call(func() Number {
		n := c.convert(x, true)
		if !c.truthy(n) || c.isnan(n) {
			return n
		}
		if _, ok := n.(Mpf); ok {
			if c.cmpop(">", n, 0) {
				return c.one()
			}
			return c.neg(c.one())
		}
		return c.binop("/", n, c.absN(n))
	})
}

func (c *Context) absN(x Number) Number {
	if t, ok := x.(Mpc); ok {
		return c.mkf(libmp.MpcAbs(t.v, c.prec, c.rnd))
	}
	return c.mkf(libmp.MpfAbs(x.(Mpf).v, c.prec, c.rnd))
}

// Agm is agm(a, b); a nil b means the default b = 1.
func (c *Context) Agm(a, b any) (Number, error) {
	return call(func() Number {
		if b == nil || c.eqOne(b) {
			return c.b1(a, libmp.MpfAgm1, libmp.MpcAgm1)
		}
		x, y := c.convert(a, true), c.convert(b, true)
		xf, okx := x.(Mpf)
		yf, oky := y.(Mpf)
		if okx && oky {
			var v libmp.MPF
			pe := libmp.Try(func() { v = libmp.MpfAgm(xf.v, yf.v, c.prec, c.rnd) })
			if pe == nil {
				return c.mkf(v)
			}
			if !libmp.Is(pe, "ComplexResult") {
				panic(pe)
			}
		}
		return c.mkc(libmp.MpcAgm(asMPC(x), asMPC(y), c.prec, c.rnd))
	})
}

// eqOne is Python's `b == 1` for any argument the API accepts.
func (c *Context) eqOne(b any) bool {
	switch t := b.(type) {
	case Number:
		return c.cmpop("==", t, 1)
	case float64:
		return t == 1
	case complex128:
		return t == 1
	case *big.Rat:
		return t.Cmp(big.NewRat(1, 1)) == 0
	case string, nil:
		return false
	}
	if v, ok := toBig(b); ok {
		return v.Cmp(big.NewInt(1)) == 0
	}
	return false
}

func asMPC(x Number) libmp.MPC {
	if t, ok := x.(Mpc); ok {
		return t.v
	}
	return libmp.MPC{Re: x.(Mpf).v, Im: libmp.Fzero}
}

// Sinc is sin(x)/x.
func (c *Context) Sinc(x any) (Number, error) {
	return c.wrapped(w1(func(x Number) Number {
		if c.isinf(x) {
			return c.binop("/", 1, x)
		}
		if !c.truthy(x) {
			return c.binop("+", x, 1)
		}
		return c.binop("/", c.sin(x), x)
	}), x)
}

// Sincpi is sin(pi*x)/(pi*x).
func (c *Context) Sincpi(x any) (Number, error) {
	return c.wrapped(w1(func(x Number) Number {
		if c.isinf(x) {
			return c.binop("/", 1, x)
		}
		if !c.truthy(x) {
			return c.binop("+", x, 1)
		}
		return c.binop("/", c.sinpi(x), c.binop("*", c.constant("pi"), x))
	}), x)
}

// Expm1 is exp(x) - 1.
func (c *Context) Expm1(x any) (Number, error) {
	return c.wrapped(w1(func(x Number) Number {
		if !c.truthy(x) {
			return c.zero()
		}
		if c.icmp("<", c.mag(x), -c.prec) {
			return c.binop("+", x, c.binop("*", 0.5, c.binop("**", x, 2)))
		}
		return c.sumAccurately(func() []any { return []any{c.exp(x), -1} }, 1)
	}), x)
}

func (c *Context) fop(op string, x, y any, prec int) Number {
	return c.farith(op, x, y, PrecOptions{Prec: &prec})
}

// Log1p is log(1+x).
func (c *Context) Log1p(x any) (Number, error) {
	return c.wrapped(w1(func(x Number) Number {
		if !c.truthy(x) {
			return c.zero()
		}
		const log1pExtraprec = 10
		cc := x
		cmag := c.mag(cc)
		a, b := c.re(cc), c.im(cc)
		wp := c.prec
		if c.icmp(">=", cmag, -wp) {
			if c.icmp("<", cmag, 4) {
				if c.icmp(">", c.mag(a), c.mag(b)) {
					wp *= 2
				} else {
					wp *= 3
				}
			}
			arg := c.fop("+", 1.0, cc, wp)
			return c.ln(arg)
		}
		real := c.binop("+", a, c.binop("*", c.binop("*", b, b), 0.5))
		if c.cmpop("<", a, 0.0) && c.icmp("<=", c.mag(real), c.iop("-", c.mag(a), log1pExtraprec)) {
			a2 := c.binop("*", a, a)
			b2 := c.fop("*", b, b, 2*wp)
			diff := c.fop("-", b2, a2, 3*wp)
			real = c.binop("+", a, c.mkf(libmp.MpfShift(mpfOf(diff), -1)))
		}
		if _, ok := x.(Mpf); ok {
			return real
		}
		return c.newMpc(real, b)
	}), x)
}

// Powm1 is x**y - 1.
func (c *Context) Powm1(x, y any) (Number, error) {
	return c.wrapped(func(a []Number) Number {
		x, y := a[0], a[1]
		one := c.one()
		w := c.binop("-", c.binop("**", x, y), one)
		M := c.mag(w)
		if c.icmp(">", M, -8) {
			return w
		}
		if !c.truthy(w) {
			if !c.truthy(y) {
				return w
			}
			unit := false
			for _, u := range []any{1, -1, complex(0, 1), complex(0, -1)} {
				if c.cmpop("==", x, u) {
					unit = true
					break
				}
			}
			if unit && c.isint(y, false) {
				return w
			}
		}
		c.binop("-", x, one)
		magy := c.mag(y)
		lnx := c.ln(x)
		if c.icmp("<", c.iop("+", magy, c.mag(lnx)), -c.prec) {
			t := c.binop("*", lnx, y)
			return c.binop("+", t, c.binop("/", c.binop("**", t, 2), 2))
		}
		return c.sumAccurately(func() []any { return []any{c.binop("**", x, y), -1} }, 1)
	}, x, y)
}

func (c *Context) rootof1(k, n int) Number {
	if n == 0 {
		libmp.Raise("ZeroDivisionError", "integer modulo by zero")
	}
	k %= n
	if k != 0 && (k < 0) != (n < 0) {
		k += n
	}
	switch {
	case k == 0:
		return c.one()
	case 2*k == n:
		return c.neg(c.one())
	case 4*k == n:
		return c.jj()
	case 4*k == 3*n:
		return c.neg(c.jj())
	}
	mk := c.mkf(libmp.MpfPos(libmp.FromInt64(int64(k)), c.prec, c.rnd))
	return c.expjpi(c.binop("/", c.binop("*", 2, mk), n))
}

func (c *Context) nthroot(x Number, n int) Number {
	if t, ok := x.(Mpf); ok {
		var r libmp.MPF
		pe := libmp.Try(func() { r = libmp.MpfNthroot(t.v, n, c.prec, c.rnd) })
		if pe == nil {
			return c.mkf(r)
		}
		if !libmp.Is(pe, "ComplexResult") || c.trapComplex {
			panic(pe)
		}
		return c.mkc(libmp.MpcNthroot(libmp.MPC{Re: t.v, Im: libmp.Fzero}, n, c.prec, c.rnd))
	}
	return c.mkc(libmp.MpcNthroot(x.(Mpc).v, n, c.prec, c.rnd))
}

func (c *Context) root(x any, n, k int) Number {
	xn := c.convert(x, true)
	if k != 0 {
		if n&1 != 0 && 2*k == n-1 && !c.truthy(c.im(xn)) && c.cmpop("<", c.re(xn), 0) {
			return c.neg(c.root(c.neg(xn), n, 0))
		}
		var v Number
		c.withPrec(c.prec+10, func() { v = c.binop("*", c.root(xn, n, 0), c.rootof1(k, n)) })
		return c.pos(v)
	}
	return c.nthroot(xn, n)
}

// Root is the k-th n-th root of x, root(x, n, k).
func (c *Context) Root(x any, n, k int) (Number, error) {
	return call(func() Number { return c.root(x, n, k) })
}

// Unitroots is unitroots(n, primitive): the n-th roots of unity.
func (c *Context) Unitroots(n int, primitive bool) ([]Number, error) {
	return call(func() []Number {
		var v []Number
		c.withPrec(c.prec+10, func() {
			for k := 0; k < n; k++ {
				if primitive && libmp.Gcd(big.NewInt(int64(k)), big.NewInt(int64(n))).Cmp(big.NewInt(1)) != 0 {
					continue
				}
				v = append(v, c.rootof1(k, n))
			}
		})
		out := make([]Number, len(v))
		for i, x := range v {
			out[i] = c.pos(x)
		}
		return out
	})
}

func (c *Context) atan2(y, x any) Number {
	yv, xv := c.realArg(y), c.realArg(x)
	return c.mkf(libmp.MpfAtan2(yv, xv, c.prec, c.rnd))
}

func (c *Context) arg(x any) Number {
	n := c.convert(x, true)
	return c.atan2(c.im(n), c.re(n))
}

// Arg is arg(x), the complex phase.
func (c *Context) Arg(x any) (Number, error) { return call(func() Number { return c.arg(x) }) }

// Fabs is abs(x) as an mpf.
func (c *Context) Fabs(x any) (Number, error) {
	return call(func() Number { return c.absN(c.convert(x, true)) })
}

// Re is the real part of x.
func (c *Context) Re(x any) (Number, error) {
	return call(func() Number { return c.re(c.convert(x, true)) })
}

// Im is the imaginary part of x.
func (c *Context) Im(x any) (Number, error) {
	return call(func() Number { return c.im(c.convert(x, true)) })
}

// Conj is the complex conjugate of x.
func (c *Context) Conj(x any) (Number, error) {
	return call(func() Number {
		n := c.convert(x, true)
		if t, ok := n.(Mpc); ok {
			return c.mkc(libmp.MpcConjugate(t.v, c.prec, c.rnd))
		}
		return n
	})
}

// Polar is (|z|, arg z).
func (c *Context) Polar(z any) (Number, Number, error) {
	type pair struct{ a, b Number }
	p, err := call(func() pair { return pair{c.absN(c.convert(z, true)), c.arg(z)} })
	return p.a, p.b, err
}

// Rect is r*(cos(phi) + i*sin(phi)).
func (c *Context) Rect(r, phi any) (Number, error) {
	return c.wrapped(func(a []Number) Number {
		cs, sn := c.cosSin(a[1], PrecOptions{}, false)
		return c.binop("*", a[0], c.newMpc(cs, sn))
	}, r, phi)
}

func (c *Context) log(x, b any) Number {
	if b == nil {
		return c.ln(x)
	}
	wp := c.prec + 20
	lx := c.builtin(x, wp, c.rnd, libmp.MpfLn, libmp.MpcLn)
	lb := c.builtin(b, wp, c.rnd, libmp.MpfLn, libmp.MpcLn)
	return c.binop("/", lx, lb)
}

// Log is log(x, b); a nil b is the natural logarithm.
func (c *Context) Log(x, b any) (Number, error) { return call(func() Number { return c.log(x, b) }) }

// Log10 is log(x, 10).
func (c *Context) Log10(x any) (Number, error) { return call(func() Number { return c.log(x, 10) }) }

// Log2 is log(x, 2).
func (c *Context) Log2(x any) (Number, error) { return call(func() Number { return c.log(x, 2) }) }

// Exp2 is 2**x.
func (c *Context) Exp2(x any) (Number, error) { return c.Power(2, x) }

// Fmod is x % y.
func (c *Context) Fmod(x, y any) (Number, error) {
	return call(func() Number { return c.binop("%", c.convert(x, true), c.convert(y, true)) })
}

// withConstant is x OP mp.degree: the constant's type is named 'constant' in
// the TypeError Python raises for an operand it cannot convert.
func (c *Context) withConstant(op string, x any, name string) Number {
	var r Number
	pe := libmp.Try(func() { r = c.binop(op, x, c.constant(name)) })
	if pe != nil {
		if pe.Class == "TypeError" {
			pe = &libmp.PyError{Class: pe.Class, Msg: strings.ReplaceAll(pe.Msg, "'mpf'", "'constant'")}
		}
		panic(pe)
	}
	return r
}

// Degrees converts radians to degrees.
func (c *Context) Degrees(x any) (Number, error) {
	return call(func() Number { return c.withConstant("/", x, "degree") })
}

// Radians converts degrees to radians.
func (c *Context) Radians(x any) (Number, error) {
	return call(func() Number { return c.withConstant("*", x, "degree") })
}
