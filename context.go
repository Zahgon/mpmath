// Copyright (c) 2005-2026 Fredrik Johansson and mpmath contributors.
// Go port distributed under the same BSD 3-clause licence (see LICENSE).

// Package mpmath is a Go port of the arbitrary-precision real and complex
// core of mpmath 1.4.1 (pure-Python backend).
package mpmath

import (
	"fmt"
	"math/big"
	"regexp"
	"strings"

	"github.com/mpmath/mpmath/internal/libmp"
)

// Error is an exception the original raises: its class name and str() text.
type Error struct {
	Class string
	Msg   string
}

// Error renders the exception as "Class: message".
func (e *Error) Error() string { return e.Class + ": " + e.Msg }

func call[T any](f func() T) (T, error) {
	var r T
	if pe := libmp.Try(func() { r = f() }); pe != nil {
		var zero T
		return zero, &Error{Class: pe.Class, Msg: pe.Msg}
	}
	return r, nil
}

// Context is an mpmath MPContext: working precision, rounding and printing
// settings shared by every value it creates.
type Context struct {
	prec, dps     int
	rnd           libmp.Rnd
	pretty        bool
	prettyReprDps bool
	trapComplex   bool
}

// NewContext returns a context with the original's defaults (53 bits,
// 15 digits, round to nearest).
func NewContext() *Context {
	return &Context{prec: 53, dps: 15, rnd: libmp.RoundNearest}
}

// Prec is ctx.prec.
func (c *Context) Prec() int { return c.prec }

// SetPrec is ctx.prec = n.
func (c *Context) SetPrec(n int) {
	c.prec = max(1, n)
	c.dps = libmp.PrecToDps(n)
}

// Dps is ctx.dps.
func (c *Context) Dps() int { return c.dps }

// SetDps is ctx.dps = n.
func (c *Context) SetDps(n int) {
	c.prec = libmp.DpsToPrec(n)
	c.dps = max(1, n)
}

var roundingMap = map[string]libmp.Rnd{
	"N": 'n', "D": 'f', "U": 'c', "Y": 'u', "Z": 'd',
	"n": 'n', "f": 'f', "c": 'c', "u": 'u', "d": 'd',
}

// Rounding is ctx.rounding: one of "n", "f", "c", "u", "d".
func (c *Context) Rounding() string { return string(rune(c.rnd)) }

// SetRounding is ctx.rounding = r; MPFR letters N, D, U, Y, Z are accepted.
func (c *Context) SetRounding(r string) error {
	v, ok := roundingMap[r]
	if !ok {
		return &Error{Class: "ValueError", Msg: "invalid rounding mode"}
	}
	c.rnd = v
	return nil
}

// Pretty is ctx.pretty.
func (c *Context) Pretty() bool { return c.pretty }

// SetPretty is ctx.pretty = v.
func (c *Context) SetPretty(v bool) { c.pretty = v }

// PrettyDps is ctx.pretty_dps: "repr" or "str".
func (c *Context) PrettyDps() string {
	if c.prettyReprDps {
		return "repr"
	}
	return "str"
}

// SetPrettyDps is ctx.pretty_dps = v.
func (c *Context) SetPrettyDps(v string) { c.prettyReprDps = v == "repr" }

// TrapComplex is ctx.trap_complex.
func (c *Context) TrapComplex() bool { return c.trapComplex }

// SetTrapComplex is ctx.trap_complex = v.
func (c *Context) SetTrapComplex(v bool) { c.trapComplex = v }

// String is str(ctx).
func (c *Context) String() string {
	pad := func(s string) string { return fmt.Sprintf("%-30s", s) }
	tc := "False"
	if c.trapComplex {
		tc = "True"
	}
	return strings.Join([]string{"Mpmath settings:",
		pad(fmt.Sprintf("  mp.prec = %d", c.prec)) + "[default: 53]",
		pad(fmt.Sprintf("  mp.dps = %d", c.dps)) + "[default: 15]",
		pad(fmt.Sprintf("  mp.rounding = '%c'", c.rnd)) + "[default: 'n']",
		pad("  mp.trap_complex = "+tc) + "[default: False]",
	}, "\n")
}

func (c *Context) reprDigits() int { return libmp.ReprDps(c.prec) }

// Number is an Mpf or an Mpc.
type Number interface {
	String() string
	Repr() string
	isNumber()
}

// Mpf is an immutable real arbitrary-precision value bound to its context.
type Mpf struct {
	ctx *Context
	v   libmp.MPF
}

// Mpc is an immutable complex arbitrary-precision value bound to its context.
type Mpc struct {
	ctx *Context
	v   libmp.MPC
}

func (Mpf) isNumber() {}
func (Mpc) isNumber() {}

func (c *Context) mkf(v libmp.MPF) Mpf { return Mpf{c, v} }
func (c *Context) mkc(v libmp.MPC) Mpc { return Mpc{c, v} }

// Raw is the value's _mpf_ tuple: (sign, mantissa, exponent, bitcount).
func (x Mpf) Raw() (int, *big.Int, int, int) {
	return x.v.Sign, new(big.Int).Set(x.v.Man), x.v.Exp, x.v.BC
}

// Real is z.real.
func (z Mpc) Real() Mpf { return Mpf{z.ctx, z.v.Re} }

// Imag is z.imag.
func (z Mpc) Imag() Mpf { return Mpf{z.ctx, z.v.Im} }

// RawMpf is a raw 4-tuple (sign, man, exp, bc) passed to mpf(), unchecked.
type RawMpf struct {
	Sign int
	Man  *big.Int
	Exp  int
	BC   int
}

// ManExp is a 2-tuple (man, exp) passed to mpf().
type ManExp struct {
	Man *big.Int
	Exp int
}

// MpfOptions are mpf()'s keyword arguments; nil fields are not passed.
type MpfOptions struct {
	Prec     *int
	Dps      *int
	Rounding *string
	Base     *int
}

// NstrOptions are nstr()'s keyword arguments; nil fields are not passed.
type NstrOptions struct {
	StripZeros       *bool
	MinFixed         *float64
	MaxFixed         *float64
	ShowZeroExponent *bool
	Base             *int
	BinaryExp        *bool
	Rounding         *string
}

func toBig(x any) (*big.Int, bool) {
	switch t := x.(type) {
	case int:
		return big.NewInt(int64(t)), true
	case int64:
		return big.NewInt(t), true
	case *big.Int:
		return t, true
	case bool:
		if t {
			return big.NewInt(1), true
		}
		return big.NewInt(0), true
	}
	return nil, false
}

// pyRepr is repr() of a value the API accepts, for error messages.
func pyRepr(x any) string {
	switch t := x.(type) {
	case nil:
		return "None"
	case string:
		return libmp.PyRepr(t)
	case float64:
		return pyFloatRepr(t)
	case complex128:
		return pyComplexRepr(t)
	case Mpf:
		return t.Repr()
	case Mpc:
		return t.Repr()
	case *big.Rat:
		return "Fraction(" + t.Num().String() + ", " + t.Denom().String() + ")"
	case ManExp:
		return fmt.Sprintf("(%s, %d)", t.Man, t.Exp)
	case RawMpf:
		return fmt.Sprintf("(%d, %s, %d, %d)", t.Sign, t.Man, t.Exp, t.BC)
	}
	if b, ok := toBig(x); ok {
		if v, isBool := x.(bool); isBool {
			if v {
				return "True"
			}
			return "False"
		}
		return b.String()
	}
	return fmt.Sprint(x)
}

// convertArg is mpf_convert_arg.
func (c *Context) convertArg(x any, prec int, rnd libmp.Rnd) libmp.MPF {
	if b, ok := toBig(x); ok {
		return libmp.FromInt(b, 0, libmp.RoundFast)
	}
	switch t := x.(type) {
	case float64:
		return libmp.FromFloat(t, 53, libmp.RoundFast)
	case Mpf:
		return t.v
	case *big.Rat:
		return libmp.FromRational(t.Num(), t.Denom(), prec, rnd)
	}
	libmp.Raise("TypeError", "%s", "cannot create mpf from "+pyRepr(x))
	return libmp.MPF{}
}

func (c *Context) newMpf(val any, o MpfOptions) libmp.MPF {
	prec, rnd, base := c.prec, c.rnd, 0
	if o.Prec != nil {
		prec = *o.Prec
	}
	if o.Dps != nil {
		prec = libmp.DpsToPrec(*o.Dps)
	}
	if o.Rounding != nil && *o.Rounding != "" {
		rnd = libmp.Rnd((*o.Rounding)[0])
	}
	if o.Base != nil {
		base = *o.Base
	}
	var v libmp.MPF
	switch t := val.(type) {
	case RawMpf:
		v = libmp.MPF{Sign: t.Sign, Man: t.Man, Exp: t.Exp, BC: t.BC}
	case ManExp:
		return libmp.FromManExp(t.Man, t.Exp, prec, rnd)
	case string:
		v = libmp.FromStr(t, prec, rnd, base)
	default:
		v = c.convertArg(val, prec, rnd)
	}
	return libmp.MpfPos(v, prec, rnd)
}

// Mpf is ctx.mpf(val).
func (c *Context) Mpf(val any) (Mpf, error) { return c.MpfWith(val, MpfOptions{}) }

// MpfWith is ctx.mpf(val, **kwargs).
func (c *Context) MpfWith(val any, o MpfOptions) (Mpf, error) {
	return call(func() Mpf { return c.mkf(c.newMpf(val, o)) })
}

var getComplex = regexp.MustCompile(`(?i)^\(?` +
	`(?P<re>[+-]?(?:\d*(?:\.\d*)?(?:e[+-]?\d+)?|\d+/\d+))??` +
	`(?P<im>[+-]?(?:\d*(?:\.\d*)?(?:e[+-]?\d+)?|\d+/\d+)\*?[ji])?` +
	`\)?\n?$`)

// convert is ctx.convert(x, strings).
func (c *Context) convert(x any, strs bool) Number {
	switch t := x.(type) {
	case Mpf:
		return t
	case Mpc:
		return t
	case float64:
		return c.mkf(libmp.FromFloat(t, 53, libmp.RoundFast))
	case complex128:
		return c.mkc(libmp.MPC{Re: libmp.FromFloat(real(t), 53, libmp.RoundFast),
			Im: libmp.FromFloat(imag(t), 53, libmp.RoundFast)})
	case *big.Rat:
		return c.mkf(libmp.FromRational(t.Num(), t.Denom(), c.prec, c.rnd))
	}
	if b, ok := toBig(x); ok {
		return c.mkf(libmp.FromInt(b, 0, libmp.RoundFast))
	}
	if s, ok := x.(string); ok && strs {
		var v libmp.MPF
		pe := libmp.Try(func() { v = libmp.FromStr(s, c.prec, c.rnd, 0) })
		if pe == nil {
			return c.mkf(v)
		}
		if !libmp.Is(pe, "ValueError") {
			panic(pe)
		}
		m := getComplex.FindStringSubmatchIndex(strings.ReplaceAll(s, " ", ""))
		if m != nil {
			src := strings.ReplaceAll(s, " ", "")
			var re any = 0
			if m[2] >= 0 && m[3] > m[2] {
				re = src[m[2]:m[3]]
			}
			if m[4] < 0 {
				libmp.Raise("AttributeError", "'NoneType' object has no attribute 'rstrip'")
			}
			im := strings.TrimRight(src[m[4]:m[5]], "jiJI*")
			return c.newMpc(c.convert(re, true), c.convert(im, true))
		}
	}
	libmp.Raise("TypeError", "%s", "cannot create mpf from "+pyRepr(x))
	return nil
}

// Convert is ctx.convert(x) (mpmathify).
func (c *Context) Convert(x any) (Number, error) {
	return call(func() Number { return c.convert(x, true) })
}

func parts(x any) (any, any) {
	switch t := x.(type) {
	case complex128:
		return real(t), imag(t)
	case Mpc:
		return Mpf{t.ctx, t.v.Re}, Mpf{t.ctx, t.v.Im}
	}
	return x, 0
}

func (c *Context) newMpc(re, im any) Mpc {
	if s, ok := re.(string); ok {
		re = c.convert(s, true)
	}
	rr, ri := parts(re)
	ir, ii := parts(im)
	vrr := c.newMpf(rr, MpfOptions{})
	vri := c.newMpf(ri, MpfOptions{})
	vir := c.newMpf(ir, MpfOptions{})
	vii := c.newMpf(ii, MpfOptions{})
	return c.mkc(libmp.MPC{
		Re: libmp.MpfSub(vrr, vii, c.prec, c.rnd),
		Im: libmp.MpfAdd(vri, vir, c.prec, c.rnd),
	})
}

// Mpc is ctx.mpc(real, imag).
func (c *Context) Mpc(re, im any) (Mpc, error) {
	return call(func() Mpc { return c.newMpc(re, im) })
}

func strOpts(rnd libmp.Rnd) libmp.StrOpts {
	o := libmp.DefaultStrOpts()
	o.Rnd = rnd
	return o
}

// String is str(x).
func (x Mpf) String() string {
	s, err := call(func() string { return libmp.ToStr(x.v, x.ctx.dps, strOpts(x.ctx.rnd)) })
	if err != nil {
		return err.Error()
	}
	return s
}

// Repr is repr(x).
func (x Mpf) Repr() string {
	c := x.ctx
	s, err := call(func() string {
		if c.pretty {
			n := c.dps
			if c.prettyReprDps {
				n = c.reprDigits()
			}
			return libmp.ToStr(x.v, n, strOpts(c.rnd))
		}
		return "mpf('" + libmp.ToStr(x.v, c.reprDigits(), strOpts(c.rnd)) + "')"
	})
	if err != nil {
		return err.Error()
	}
	return s
}

// String is str(z).
func (z Mpc) String() string {
	return "(" + libmp.MpcToStr(z.v, z.ctx.dps, libmp.DefaultStrOpts()) + ")"
}

// Repr is repr(z).
func (z Mpc) Repr() string {
	c := z.ctx
	if c.pretty {
		n := c.dps
		if c.prettyReprDps {
			n = c.reprDigits()
		}
		return "(" + libmp.MpcToStr(z.v, n, libmp.DefaultStrOpts()) + ")"
	}
	r := z.Real().Repr()
	i := z.Imag().Repr()
	return "mpc(real=" + r[4:len(r)-1] + ", imag=" + i[4:len(i)-1] + ")"
}

func (o NstrOptions) lib() libmp.StrOpts {
	s := libmp.DefaultStrOpts()
	if o.StripZeros != nil {
		s.StripZeros = *o.StripZeros
	}
	s.MinFixed, s.MaxFixed = o.MinFixed, o.MaxFixed
	if o.ShowZeroExponent != nil {
		s.ShowZeroExponent = *o.ShowZeroExponent
	}
	if o.Base != nil {
		s.Base = *o.Base
	}
	if o.BinaryExp != nil {
		s.BinaryExp = *o.BinaryExp
	}
	if o.Rounding != nil {
		s.Rnd = 0
		if *o.Rounding != "" {
			s.Rnd = libmp.Rnd((*o.Rounding)[0])
		}
		if len(*o.Rounding) != 1 {
			s.Rnd = '?'
		}
	}
	return s
}

// Nstr is ctx.nstr(x, n, **kwargs) for an Mpf or an Mpc.
func (c *Context) Nstr(x Number, n int, o NstrOptions) (string, error) {
	return call(func() string {
		switch t := x.(type) {
		case Mpf:
			return libmp.ToStr(t.v, n, o.lib())
		case Mpc:
			return "(" + libmp.MpcToStr(t.v, n, o.lib()) + ")"
		}
		return fmt.Sprint(x)
	})
}
