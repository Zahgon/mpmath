// Copyright (c) 2005-2026 Fredrik Johansson and mpmath contributors.
// Go port distributed under the same BSD 3-clause licence (see LICENSE).

package mpmath

import (
	"math/big"

	"github.com/mpmath/mpmath/internal/libmp"
)

// PrecOptions are the precision keywords of fadd, fsub, fmul, fdiv, fneg,
// cos_sin and cospi_sinpi. Each field is passed only when it is set; PrecInf
// and DpsInf stand for prec=inf and dps=inf.
type PrecOptions struct {
	Prec, Dps       *int
	PrecInf, DpsInf bool
	Rounding        *string
	Exact           *bool
}

func (o PrecOptions) given() bool {
	return o.Prec != nil || o.Dps != nil || o.PrecInf || o.DpsInf || o.Rounding != nil || o.Exact != nil
}

var mpfrRounding = map[string]libmp.Rnd{"N": 'n', "D": 'f', "U": 'c', "Y": 'u', "Z": 'd',
	"n": 'n', "f": 'f', "c": 'c', "u": 'u', "d": 'd'}

// parsePrec is ctx._parse_prec.
func (c *Context) parsePrec(o PrecOptions) (int, libmp.Rnd) {
	if !o.given() {
		return c.prec, c.rnd
	}
	if o.Exact != nil && *o.Exact {
		return 0, libmp.RoundFloor
	}
	prec, rnd := c.prec, c.rnd
	if o.Rounding != nil {
		r, ok := mpfrRounding[*o.Rounding]
		if !ok {
			libmp.Raise("KeyError", "%s", libmp.PyRepr(*o.Rounding))
		}
		rnd = r
	}
	switch {
	case o.PrecInf:
		return 0, libmp.RoundFloor
	case o.Prec != nil:
		prec = *o.Prec
	case o.DpsInf:
		return 0, libmp.RoundFloor
	case o.Dps != nil:
		prec = libmp.DpsToPrec(*o.Dps)
	}
	return prec, rnd
}

const exactOverflowMsg = "the exact result does not fit in memory"

func (c *Context) farith(op string, x, y any, o PrecOptions) Number {
	prec, rnd := c.parsePrec(o)
	if op == "/" && prec == 0 {
		libmp.Raise("ValueError", "division is not an exact operation")
	}
	xn := c.convert(x, true)
	yn := c.convert(y, true)
	var r Number
	pe := libmp.Try(func() {
		xf, xIsF := xn.(Mpf)
		yf, yIsF := yn.(Mpf)
		xc, _ := xn.(Mpc)
		yc, _ := yn.(Mpc)
		switch op {
		case "+":
			switch {
			case xIsF && yIsF:
				r = c.mkf(libmp.MpfAdd(xf.v, yf.v, prec, rnd))
			case xIsF:
				r = c.mkc(libmp.MpcAddMpf(yc.v, xf.v, prec, rnd))
			case yIsF:
				r = c.mkc(libmp.MpcAddMpf(xc.v, yf.v, prec, rnd))
			default:
				r = c.mkc(libmp.MpcAdd(xc.v, yc.v, prec, rnd))
			}
		case "-":
			switch {
			case xIsF && yIsF:
				r = c.mkf(libmp.MpfSub(xf.v, yf.v, prec, rnd))
			case xIsF:
				r = c.mkc(libmp.MpcMpfSub(xf.v, yc.v, prec, rnd))
			case yIsF:
				r = c.mkc(libmp.MpcSubMpf(xc.v, yf.v, prec, rnd))
			default:
				r = c.mkc(libmp.MpcSub(xc.v, yc.v, prec, rnd))
			}
		case "*":
			switch {
			case xIsF && yIsF:
				r = c.mkf(libmp.MpfMul(xf.v, yf.v, prec, rnd))
			case xIsF:
				r = c.mkc(libmp.MpcMulMpf(yc.v, xf.v, prec, rnd))
			case yIsF:
				r = c.mkc(libmp.MpcMulMpf(xc.v, yf.v, prec, rnd))
			default:
				r = c.mkc(libmp.MpcMul(xc.v, yc.v, prec, rnd))
			}
		case "/":
			switch {
			case xIsF && yIsF:
				r = c.mkf(libmp.MpfDiv(xf.v, yf.v, prec, rnd))
			case xIsF:
				r = c.mkc(libmp.MpcMpfDiv(xf.v, yc.v, prec, rnd))
			case yIsF:
				r = c.mkc(libmp.MpcDivMpf(xc.v, yf.v, prec, rnd))
			default:
				r = c.mkc(libmp.MpcDiv(xc.v, yc.v, prec, rnd))
			}
		}
	})
	if pe != nil {
		if op != "/" && (libmp.Is(pe, "ValueError") || libmp.Is(pe, "OverflowError")) {
			libmp.Raise("OverflowError", exactOverflowMsg)
		}
		panic(pe)
	}
	return r
}

// Fadd is ctx.fadd(x, y, **kwargs).
func (c *Context) Fadd(x, y any, o PrecOptions) (Number, error) {
	return call(func() Number { return c.farith("+", x, y, o) })
}

// Fsub is ctx.fsub(x, y, **kwargs).
func (c *Context) Fsub(x, y any, o PrecOptions) (Number, error) {
	return call(func() Number { return c.farith("-", x, y, o) })
}

// Fmul is ctx.fmul(x, y, **kwargs).
func (c *Context) Fmul(x, y any, o PrecOptions) (Number, error) {
	return call(func() Number { return c.farith("*", x, y, o) })
}

// Fdiv is ctx.fdiv(x, y, **kwargs).
func (c *Context) Fdiv(x, y any, o PrecOptions) (Number, error) {
	return call(func() Number { return c.farith("/", x, y, o) })
}

// Fneg is ctx.fneg(x, **kwargs).
func (c *Context) Fneg(x any, o PrecOptions) (Number, error) {
	return call(func() Number {
		prec, rnd := c.parsePrec(o)
		switch t := c.convert(x, true).(type) {
		case Mpf:
			return c.mkf(libmp.MpfNeg(t.v, prec, rnd))
		case Mpc:
			return c.mkc(libmp.MpcNeg(t.v, prec, rnd))
		}
		return nil
	})
}

// Fsum is ctx.fsum(terms, absolute, squared).
func (c *Context) Fsum(terms []any, absolute, squared bool) (Number, error) {
	return call(func() Number {
		prec, rnd := c.prec, c.rnd
		var real, imag []libmp.MPF
		for _, term := range terms {
			n, ok := term.(Number)
			if !ok {
				n = c.convert(term, true)
			}
			if f, isF := n.(Mpf); isF {
				re := f.v
				if squared {
					re = libmp.MpfMul(re, re, 0, libmp.RoundFast)
				} else if absolute {
					re = libmp.MpfAbs(re, 0, libmp.RoundFast)
				}
				real = append(real, re)
				continue
			}
			z := n.(Mpc).v
			switch {
			case squared && absolute:
				real = append(real, libmp.MpfMul(z.Re, z.Re, 0, libmp.RoundFast),
					libmp.MpfMul(z.Im, z.Im, 0, libmp.RoundFast))
			case squared:
				w := libmp.MpcPowInt(z, big.NewInt(2), prec+10, libmp.RoundFast)
				real, imag = append(real, w.Re), append(imag, w.Im)
			case absolute:
				real = append(real, libmp.MpcAbs(z, prec, libmp.RoundFast))
			default:
				real, imag = append(real, z.Re), append(imag, z.Im)
			}
		}
		s := libmp.MpfSum(real, prec, rnd, absolute)
		if len(imag) > 0 {
			return c.mkc(libmp.MPC{Re: s, Im: libmp.MpfSum(imag, prec, rnd, false)})
		}
		return c.mkf(s)
	})
}

// Fdot is ctx.fdot(A, B, conjugate) for two equally long lists (zip stops at
// the shorter one).
func (c *Context) Fdot(a, b []any, conjugate bool) (Number, error) {
	return call(func() Number {
		prec, rnd := c.prec, c.rnd
		f := libmp.RoundFast
		mul := func(x, y libmp.MPF) libmp.MPF { return libmp.MpfMul(x, y, 0, f) }
		var real, imag []libmp.MPF
		for i := 0; i < len(a) && i < len(b); i++ {
			x, ok := a[i].(Number)
			if !ok {
				x = c.convert(a[i], true)
			}
			y, ok := b[i].(Number)
			if !ok {
				y = c.convert(b[i], true)
			}
			xf, xr := x.(Mpf)
			yf, yr := y.(Mpf)
			xc, _ := x.(Mpc)
			yc, _ := y.(Mpc)
			switch {
			case xr && yr:
				real = append(real, mul(xf.v, yf.v))
			case xr:
				bim := yc.v.Im
				if conjugate {
					bim = libmp.MpfNeg(bim, 0, f)
				}
				real = append(real, mul(xf.v, yc.v.Re))
				imag = append(imag, mul(xf.v, bim))
			case yr:
				real = append(real, mul(xc.v.Re, yf.v))
				imag = append(imag, mul(xc.v.Im, yf.v))
			default:
				bim := yc.v.Im
				if conjugate {
					bim = libmp.MpfNeg(bim, 0, f)
				}
				real = append(real, mul(xc.v.Re, yc.v.Re), libmp.MpfNeg(mul(xc.v.Im, bim), 0, f))
				imag = append(imag, mul(xc.v.Re, bim), mul(xc.v.Im, yc.v.Re))
			}
		}
		s := libmp.MpfSum(real, prec, rnd, false)
		if len(imag) > 0 {
			return c.mkc(libmp.MPC{Re: s, Im: libmp.MpfSum(imag, prec, rnd, false)})
		}
		return c.mkf(s)
	})
}

// Fprod is ctx.fprod(factors).
func (c *Context) Fprod(factors []any) (Number, error) {
	return call(func() Number {
		var v Number = c.mkf(libmp.Fone)
		for _, p := range factors {
			v = c.binop2("*", v, p, true)
		}
		if t, ok := v.(Mpc); ok {
			return c.mkc(libmp.MpcPos(t.v, c.prec, c.rnd))
		}
		return c.mkf(libmp.MpfPos(v.(Mpf).v, c.prec, c.rnd))
	})
}

// Power is ctx.power(x, y).
func (c *Context) Power(x, y any) (Number, error) {
	return call(func() Number { return c.binop("**", c.convert(x, true), c.convert(y, true)) })
}

func noAttr(n Number, attr string) {
	libmp.Raise("AttributeError", "%s", "'"+typeName(n)+"' object has no attribute '"+attr+"'")
}

func (c *Context) realArg(x any) libmp.MPF { return mpfOf(c.convert(x, true)) }

func mpfOf(n Number) libmp.MPF {
	f, ok := n.(Mpf)
	if !ok {
		noAttr(n, "_mpf_")
	}
	return f.v
}

// Ldexp is ctx.ldexp(x, n).
func (c *Context) Ldexp(x any, n int) (Mpf, error) {
	return call(func() Mpf { return c.mkf(libmp.MpfShift(c.realArg(x), n)) })
}

// Frexp is ctx.frexp(x).
func (c *Context) Frexp(x any) (Mpf, int, error) {
	type pair struct {
		y Mpf
		n int
	}
	p, err := call(func() pair { y, n := libmp.MpfFrexp(c.realArg(x)); return pair{c.mkf(y), n} })
	return p.y, p.n, err
}

// Hypot is ctx.hypot(x, y).
func (c *Context) Hypot(x, y any) (Mpf, error) {
	return call(func() Mpf {
		xn, yn := c.convert(x, true), c.convert(y, true)
		a := mpfOf(xn)
		b := mpfOf(yn)
		return c.mkf(libmp.MpfHypot(a, b, c.prec, c.rnd))
	})
}

// Atan2 is ctx.atan2(y, x).
func (c *Context) Atan2(y, x any) (Mpf, error) {
	return call(func() Mpf {
		xn, yn := c.convert(x, true), c.convert(y, true)
		yv := mpfOf(yn)
		xv := mpfOf(xn)
		return c.mkf(libmp.MpfAtan2(yv, xv, c.prec, c.rnd))
	})
}

func (c *Context) cosSin(x any, o PrecOptions, pi bool) (Number, Number) {
	n, ok := x.(Number)
	if !ok {
		n = c.convert(x, true)
	}
	prec, rnd := c.parsePrec(o)
	switch t := n.(type) {
	case Mpf:
		var cs, sn libmp.MPF
		if pi {
			cs, sn = libmp.MpfCosSinPi(t.v, prec, rnd)
		} else {
			cs, sn = libmp.MpfCosSin(t.v, prec, rnd)
		}
		return c.mkf(cs), c.mkf(sn)
	}
	var cs, sn libmp.MPC
	if pi {
		cs, sn = libmp.MpcCosSinPi(n.(Mpc).v, prec, rnd)
	} else {
		cs, sn = libmp.MpcCosSin(n.(Mpc).v, prec, rnd)
	}
	return c.mkc(cs), c.mkc(sn)
}

// CosSin is ctx.cos_sin(x, **kwargs).
func (c *Context) CosSin(x any, o PrecOptions) (Number, Number, error) {
	type pair struct{ a, b Number }
	p, err := call(func() pair { a, b := c.cosSin(x, o, false); return pair{a, b} })
	return p.a, p.b, err
}

// CospiSinpi is ctx.cospi_sinpi(x, **kwargs).
func (c *Context) CospiSinpi(x any, o PrecOptions) (Number, Number, error) {
	type pair struct{ a, b Number }
	p, err := call(func() pair { a, b := c.cosSin(x, o, true); return pair{a, b} })
	return p.a, p.b, err
}

func (c *Context) ninf() Mpf { return c.mkf(libmp.Fninf) }

// mpfMag is ctx._mpf_mag: an int, or an mpf for zero (-inf), inf and nan.
func (c *Context) mpfMag(x libmp.MPF) any {
	if x.Eq(libmp.Fzero) {
		return c.ninf()
	}
	if x.IsSpecial() {
		return c.mkf(libmp.MpfAbs(x, 0, libmp.RoundFast))
	}
	return x.Exp + x.Man.BitLen()
}

func intOrBig(n int) *big.Int { return big.NewInt(int64(n)) }

func magArg(v any) any {
	if n, ok := v.(int); ok {
		return intOrBig(n)
	}
	return v
}

// pyMax is max(a, b) over ints and mpfs: a unless b > a.
func (c *Context) pyMax(a, b any) any {
	ai, aInt := a.(int)
	bi, bInt := b.(int)
	if aInt && bInt {
		if bi > ai {
			return b
		}
		return a
	}
	if c.cmpop(">", magArg(b), magArg(a)) {
		return b
	}
	return a
}

func (c *Context) mag(x any) any {
	switch t := x.(type) {
	case Mpf:
		return c.mpfMag(t.v)
	case Mpc:
		r, i := t.v.Re, t.v.Im
		if r.Eq(libmp.Fzero) {
			return c.mpfMag(i)
		}
		if i.Eq(libmp.Fzero) {
			return c.mpfMag(r)
		}
		m := c.pyMax(c.mpfMag(r), c.mpfMag(i))
		if n, ok := m.(int); ok {
			return 1 + n
		}
		return c.binop("+", 1, m)
	case *big.Rat:
		if t.Num().Sign() == 0 {
			return c.ninf()
		}
		return 1 + t.Num().BitLen() - t.Denom().BitLen()
	}
	if b, ok := toBig(x); ok {
		if b.Sign() != 0 {
			return b.BitLen()
		}
		return c.ninf()
	}
	return c.mag(c.convert(x, true))
}

// Mag is ctx.mag(x): an int, or an Mpf (-inf, inf or nan) for zero and the
// special values.
func (c *Context) Mag(x any) (any, error) {
	return call(func() any { return c.mag(x) })
}

// NintDistance is ctx.nint_distance(x): the nearest integer n and d, an int
// or the Mpf -inf.
func (c *Context) NintDistance(x any) (*big.Int, any, error) {
	type pair struct {
		n *big.Int
		d any
	}
	p, err := call(func() pair { n, d := c.nintDistance(x); return pair{n, d} })
	return p.n, p.d, err
}

func (c *Context) nintDistance(x any) (*big.Int, any) {
	if b, ok := toBig(x); ok {
		if _, isBool := x.(bool); !isBool {
			return b, c.ninf()
		}
	}
	if q, ok := x.(*big.Rat); ok {
		p, d := q.Num(), q.Denom()
		n, r := new(big.Int).DivMod(p, d, new(big.Int))
		if new(big.Int).Lsh(r, 1).Cmp(d) >= 0 {
			n.Add(n, big.NewInt(1))
		} else if r.Sign() == 0 {
			return n, c.ninf()
		}
		t := new(big.Int).Sub(p, new(big.Int).Mul(n, d))
		return n, t.BitLen() - d.BitLen()
	}
	var re libmp.MPF
	var imDist any = c.ninf()
	switch t := x.(type) {
	case Mpf:
		re = t.v
	case Mpc:
		re = t.v.Re
		iman, iexp := libmp.ToManExp(t.v.Im, true)
		if iman.Sign() != 0 {
			imDist = iexp + iman.BitLen()
		}
	default:
		return c.nintDistance(c.convert(x, true))
	}
	man, exp := libmp.ToManExp(re, true)
	mag := exp + man.BitLen()
	var n *big.Int
	var reDist any
	switch {
	case mag < 0:
		n, reDist = big.NewInt(0), mag
	case man.Sign() != 0:
		sign := man.Sign() < 0
		man = new(big.Int).Abs(man)
		switch {
		case exp >= 0:
			n, reDist = new(big.Int).Lsh(man, uint(exp)), c.ninf()
		case exp == -1:
			n = new(big.Int).Add(new(big.Int).Rsh(man, 1), big.NewInt(1))
			reDist = 0
		default:
			d := uint(-exp - 1)
			t := new(big.Int).Rsh(man, d)
			if t.Bit(0) == 1 {
				t.Add(t, big.NewInt(1))
				man = new(big.Int).Sub(new(big.Int).Lsh(t, d), man)
			} else {
				man = new(big.Int).Sub(man, new(big.Int).Lsh(t, d))
			}
			n = new(big.Int).Rsh(t, 1)
			reDist = exp + man.BitLen()
		}
		if sign {
			n.Neg(n)
		}
	default:
		n, reDist = big.NewInt(0), c.ninf()
	}
	return n, c.pyMax(reDist, imDist)
}

// truthy is Python's bool(x) for the values the API accepts.
func (c *Context) truthy(x any) bool {
	switch t := x.(type) {
	case nil:
		return false
	case string:
		return t != ""
	case float64:
		return t != 0
	case complex128:
		return t != 0
	case *big.Rat:
		return t.Sign() != 0
	case Number:
		return c.Bool(t)
	case RawMpf, ManExp:
		return true
	}
	b, _ := toBig(x)
	return b.Sign() != 0
}

func (c *Context) isinf(x any) bool {
	switch t := x.(type) {
	case Mpf:
		return t.v.IsInf()
	case Mpc:
		return libmp.MpcIsInf(t.v)
	case *big.Rat:
		return false
	}
	if _, ok := toBig(x); ok {
		return false
	}
	return c.isinf(c.convert(x, true))
}

func (c *Context) isnan(x any) bool {
	switch t := x.(type) {
	case Mpf:
		return t.v.Eq(libmp.Fnan)
	case Mpc:
		return t.v.Re.Eq(libmp.Fnan) || t.v.Im.Eq(libmp.Fnan)
	case *big.Rat:
		return false
	}
	if _, ok := toBig(x); ok {
		return false
	}
	return c.isnan(c.convert(x, true))
}

func isIntRaw(v libmp.MPF) bool { return (v.Man.Sign() != 0 && v.Exp >= 0) || v.Eq(libmp.Fzero) }

func (c *Context) isint(x any, gaussian bool) bool {
	if _, ok := toBig(x); ok {
		return true
	}
	switch t := x.(type) {
	case Mpf:
		return !c.isinf(t) && !c.isnan(t) && isIntRaw(t.v)
	case Mpc:
		if c.isinf(t) || c.isnan(t) {
			return false
		}
		if gaussian {
			return isIntRaw(t.v.Re) && isIntRaw(t.v.Im)
		}
		return isIntRaw(t.v.Re) && t.v.Im.Eq(libmp.Fzero)
	case *big.Rat:
		return t.IsInt()
	}
	return c.isint(c.convert(x, true), gaussian)
}

func (c *Context) isnormal(x any) bool {
	switch t := x.(type) {
	case Mpf:
		return !t.v.IsSpecial() && t.v.Man.Sign() != 0
	case Mpc:
		re, im := t.v.Re, t.v.Im
		reN, imN := re.Man.Sign() != 0, im.Man.Sign() != 0
		if re.Eq(libmp.Fzero) {
			return imN
		}
		if im.Eq(libmp.Fzero) {
			return reN
		}
		return reN && imN
	case *big.Rat:
		return t.Sign() != 0
	}
	if b, ok := toBig(x); ok {
		return b.Sign() != 0
	}
	return c.isnormal(c.convert(x, true))
}

func (c *Context) isnpint(x any) bool {
	if !c.truthy(x) {
		return true
	}
	switch t := x.(type) {
	case Mpf:
		if c.isinf(t) || c.isnan(t) {
			return false
		}
		return t.v.Sign == 1 && t.v.Exp >= 0
	case Mpc:
		return t.v.Im.Eq(libmp.Fzero) && c.isnpint(c.mkf(t.v.Re))
	case *big.Rat:
		return t.IsInt() && t.Sign() <= 0
	}
	if b, ok := toBig(x); ok {
		if _, isBool := x.(bool); !isBool {
			return b.Sign() <= 0
		}
	}
	return c.isnpint(c.convert(x, true))
}

func pred(f func() bool) (bool, error) { return call(f) }

// IsInf is ctx.isinf(x).
func (c *Context) IsInf(x any) (bool, error) { return pred(func() bool { return c.isinf(x) }) }

// IsNan is ctx.isnan(x).
func (c *Context) IsNan(x any) (bool, error) { return pred(func() bool { return c.isnan(x) }) }

// IsFinite is ctx.isfinite(x).
func (c *Context) IsFinite(x any) (bool, error) {
	return pred(func() bool { return !(c.isinf(x) || c.isnan(x)) })
}

// IsNormal is ctx.isnormal(x).
func (c *Context) IsNormal(x any) (bool, error) { return pred(func() bool { return c.isnormal(x) }) }

// IsInt is ctx.isint(x, gaussian).
func (c *Context) IsInt(x any, gaussian bool) (bool, error) {
	return pred(func() bool { return c.isint(x, gaussian) })
}

// IsNpint is ctx.isnpint(x).
func (c *Context) IsNpint(x any) (bool, error) { return pred(func() bool { return c.isnpint(x) }) }

// AlmostEq is ctx.almosteq(s, t, rel_eps, abs_eps); nil epsilons are None.
func (c *Context) AlmostEq(s, t, relEps, absEps any) (bool, error) {
	return pred(func() bool {
		sn, tn := c.convert(s, true), c.convert(t, true)
		for _, v := range []Number{sn, tn} {
			if c.isinf(v) || c.isnan(v) {
				return c.cmpop("==", sn, tn)
			}
		}
		if absEps == nil && relEps == nil {
			e := c.mkf(libmp.MpfShift(libmp.Fone, -c.prec+4))
			relEps, absEps = e, e
		}
		if absEps == nil {
			absEps = relEps
		} else if relEps == nil {
			relEps = absEps
		}
		absOf := func(n Number) Number {
			if z, ok := n.(Mpc); ok {
				return c.mkf(libmp.MpcAbs(z.v, c.prec, c.rnd))
			}
			return c.mkf(libmp.MpfAbs(n.(Mpf).v, c.prec, c.rnd))
		}
		diff := absOf(c.binop("-", sn, tn))
		if c.cmpop("<=", diff, absEps) {
			return true
		}
		abss, abst := absOf(sn), absOf(tn)
		var err Number
		if c.cmpop("<", abss, abst) {
			err = c.binop("/", diff, abst)
		} else {
			err = c.binop("/", diff, abss)
		}
		return c.cmpop("<=", err, relEps)
	})
}
