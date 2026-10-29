// Copyright (c) 2005-2026 Fredrik Johansson and mpmath contributors.
// Go port distributed under the same BSD 3-clause licence (see LICENSE).

package mpmath

import (
	"math/big"

	"github.com/mpmath/mpmath/internal/libmp"
)

// Arithmetic and comparison of mpf and mpc values with Python's operator
// protocol: the left operand's method, then the right operand's reflected
// method, and a TypeError when both decline.

func typeName(x any) string {
	switch x.(type) {
	case nil:
		return "NoneType"
	case bool:
		return "bool"
	case int, int64, *big.Int:
		return "int"
	case float64:
		return "float"
	case complex128:
		return "complex"
	case *big.Rat:
		return "Fraction"
	case string:
		return "str"
	case Mpf:
		return "mpf"
	case Mpc:
		return "mpc"
	case RawMpf, ManExp:
		return "tuple"
	}
	return "object"
}

func isSeq(x any) bool {
	switch x.(type) {
	case string, RawMpf, ManExp:
		return true
	}
	return false
}

// tryConvert is ctx.convert(x, strings=False) where the caller catches
// TypeError (and ValueError when catchValue is set).
func (c *Context) tryConvert(x any, catchValue bool) (Number, bool) {
	var n Number
	pe := libmp.Try(func() { n = c.convert(x, false) })
	if pe == nil {
		return n, true
	}
	if pe.Class == "TypeError" || (catchValue && libmp.Is(pe, "ValueError")) {
		return nil, false
	}
	panic(pe)
}

func opSymbol(op string) string {
	switch op {
	case "**":
		return "** or pow()"
	case "divmod":
		return "divmod()"
	}
	return op
}

func unsupported(op string, x, y any) {
	libmp.Raise("TypeError", "%s", "unsupported operand type(s) for "+opSymbol(op)+": '"+
		typeName(x)+"' and '"+typeName(y)+"'")
}

// mpfForward is mpf.__op__(self, other); ok is false for NotImplemented.
func (c *Context) mpfForward(op string, s Mpf, other any) (Number, bool) {
	p, r := c.prec, c.rnd
	sv := s.v
	switch o := other.(type) {
	case Mpf:
		switch op {
		case "+":
			return c.mkf(libmp.MpfAdd(sv, o.v, p, r)), true
		case "-":
			return c.mkf(libmp.MpfSub(sv, o.v, p, r)), true
		case "*":
			return c.mkf(libmp.MpfMul(sv, o.v, p, r)), true
		case "/":
			return c.mkf(libmp.MpfDiv(sv, o.v, p, r)), true
		case "%":
			return c.mkf(libmp.MpfMod(sv, o.v, p, r)), true
		case "**":
			var v libmp.MPF
			pe := libmp.Try(func() { v = libmp.MpfPow(sv, o.v, p, r) })
			if pe == nil {
				return c.mkf(v), true
			}
			if pe.Class != "ComplexResult" || c.trapComplex {
				panic(pe)
			}
			return c.mkc(libmp.MpcPowMpf(libmp.MPC{Re: sv, Im: libmp.Fzero}, o.v, p, r)), true
		}
	case Mpc:
		switch op {
		case "+":
			return c.mkc(libmp.MpcAddMpf(o.v, sv, p, r)), true
		case "-":
			return c.mkc(libmp.MpcMpfSub(sv, o.v, p, r)), true
		case "*":
			return c.mkc(libmp.MpcMulMpf(o.v, sv, p, r)), true
		case "/":
			return c.mkc(libmp.MpcMpfDiv(sv, o.v, p, r)), true
		case "%":
			return nil, false
		case "**":
			return c.mkc(libmp.MpcPow(libmp.MPC{Re: sv, Im: libmp.Fzero}, o.v, p, r)), true
		}
	}
	switch op {
	case "//":
		m := c.binop("%", s, other)
		return c.binop("/", c.binop("-", s, m), other), true
	}
	conv, ok := c.tryConvert(other, false)
	if !ok {
		return nil, false
	}
	return c.mpfForward(op, s, conv)
}

// mpfReflected is mpf.__rop__(self, other) for other OP self.
func (c *Context) mpfReflected(op string, s Mpf, other any) (Number, bool) {
	switch op {
	case "+", "*":
		return c.mpfForward(op, s, other)
	case "-", "/", "%", "**":
		o, ok := c.tryConvert(other, true)
		if !ok {
			return nil, false
		}
		if _, isMpc := o.(Mpc); isMpc && op == "%" {
			libmp.Raise("RecursionError", "maximum recursion depth exceeded")
		}
		return c.binop(op, o, s), true
	}
	return nil, false
}

// mpcForward is mpc.__op__(self, other).
func (c *Context) mpcForward(op string, s Mpc, other any) (Number, bool) {
	p, r := c.prec, c.rnd
	sv := s.v
	switch op {
	case "%", "//", "divmod":
		return nil, false
	case "*":
		if b, isInt := toBig(other); isInt {
			return c.mkc(libmp.MpcMulInt(sv, b, p, r)), true
		}
	case "**":
		if b, isInt := toBig(other); isInt {
			return c.mkc(libmp.MpcPowInt(sv, b, p, r)), true
		}
	}
	o, isMpc := other.(Mpc)
	if !isMpc {
		conv, ok := c.tryConvert(other, true)
		if !ok {
			return nil, false
		}
		if f, isMpf := conv.(Mpf); isMpf {
			switch op {
			case "+":
				return c.mkc(libmp.MpcAddMpf(sv, f.v, p, r)), true
			case "-":
				return c.mkc(libmp.MpcSubMpf(sv, f.v, p, r)), true
			case "*":
				return c.mkc(libmp.MpcMulMpf(sv, f.v, p, r)), true
			case "/":
				return c.mkc(libmp.MpcDivMpf(sv, f.v, p, r)), true
			case "**":
				return c.mkc(libmp.MpcPowMpf(sv, f.v, p, r)), true
			}
		}
		o = conv.(Mpc)
	}
	switch op {
	case "+":
		return c.mkc(libmp.MpcAdd(sv, o.v, p, r)), true
	case "-":
		return c.mkc(libmp.MpcSub(sv, o.v, p, r)), true
	case "*":
		return c.mkc(libmp.MpcMul(sv, o.v, p, r)), true
	case "/":
		return c.mkc(libmp.MpcDiv(sv, o.v, p, r)), true
	}
	return c.mkc(libmp.MpcPow(sv, o.v, p, r)), true
}

// mpcReflected is mpc.__rop__(self, other).
func (c *Context) mpcReflected(op string, s Mpc, other any) (Number, bool) {
	switch op {
	case "+":
		return c.mpcForward(op, s, other)
	case "*":
		if b, isInt := toBig(other); isInt {
			return c.mkc(libmp.MpcMulInt(s.v, b, c.prec, c.rnd)), true
		}
		o, ok := c.tryConvert(other, true)
		if !ok {
			return nil, false
		}
		return c.binop(op, o, s), true
	case "-", "/", "**":
		o, ok := c.tryConvert(other, true)
		if !ok {
			return nil, false
		}
		return c.binop(op, o, s), true
	}
	return nil, false
}

func (c *Context) forward(op string, x, y any) (Number, bool) {
	switch t := x.(type) {
	case Mpf:
		return c.mpfForward(op, t, y)
	case Mpc:
		return c.mpcForward(op, t, y)
	}
	return nil, false
}

func (c *Context) reflected(op string, x, y any) (Number, bool) {
	switch t := y.(type) {
	case Mpf:
		return c.mpfReflected(op, t, x)
	case Mpc:
		return c.mpcReflected(op, t, x)
	}
	return nil, false
}

func sameType(x, y any) bool { return typeName(x) == typeName(y) }

// binop evaluates x OP y as Python does.
func (c *Context) binop(op string, x, y any) Number { return c.binop2(op, x, y, false) }

// binop2 also covers the in-place form x OP= y, which reports "OP=" and has no
// sequence fallback.
func (c *Context) binop2(op string, x, y any, inplace bool) Number {
	_, xmp := x.(Number)
	if !xmp {
		if _, ymp := y.(Number); !ymp {
			libmp.Raise("TypeError", "at least one operand must be an mpf or mpc")
		}
	}
	if r, ok := c.forward(op, x, y); ok {
		return r
	}
	if !sameType(x, y) {
		if r, ok := c.reflected(op, x, y); ok {
			return r
		}
	}
	if _, isMpc := x.(Mpc); isMpc {
		if _, isRat := y.(*big.Rat); isRat && (op == "%" || op == "//" || op == "divmod") {
			unsupported(op, complex128(0), complex128(0))
		}
	}
	if inplace {
		unsupported(op+"=", x, y)
	}
	switch op {
	case "+":
		if _, isStr := x.(string); isStr {
			libmp.Raise("TypeError", "%s", "can only concatenate str (not \""+typeName(y)+"\") to str")
		}
		if isSeq(x) {
			libmp.Raise("TypeError", "%s", "can only concatenate tuple (not \""+typeName(y)+"\") to tuple")
		}
	case "*":
		if isSeq(x) {
			libmp.Raise("TypeError", "%s", "can't multiply sequence by non-int of type '"+typeName(y)+"'")
		}
		if isSeq(y) {
			libmp.Raise("TypeError", "%s", "can't multiply sequence by non-int of type '"+typeName(x)+"'")
		}
	}
	unsupported(op, x, y)
	return nil
}

func (c *Context) binopErr(op string, x, y any) (Number, error) {
	return call(func() Number { return c.binop(op, x, y) })
}

// Add is x + y.
func (c *Context) Add(x, y any) (Number, error) { return c.binopErr("+", x, y) }

// Sub is x - y.
func (c *Context) Sub(x, y any) (Number, error) { return c.binopErr("-", x, y) }

// Mul is x * y.
func (c *Context) Mul(x, y any) (Number, error) { return c.binopErr("*", x, y) }

// Div is x / y.
func (c *Context) Div(x, y any) (Number, error) { return c.binopErr("/", x, y) }

// Mod is x % y.
func (c *Context) Mod(x, y any) (Number, error) { return c.binopErr("%", x, y) }

// FloorDiv is x // y.
func (c *Context) FloorDiv(x, y any) (Number, error) { return c.binopErr("//", x, y) }

// Pow is x ** y.
func (c *Context) Pow(x, y any) (Number, error) { return c.binopErr("**", x, y) }

// DivMod is divmod(x, y).
func (c *Context) DivMod(x, y any) (Number, Number, error) {
	var q, m Number
	_, err := call(func() bool {
		if s, ok := x.(Mpf); ok {
			m = c.binop("%", s, y)
			q = c.binop("/", c.binop("-", s, m), y)
			return true
		}
		if _, isMpc := x.(Mpc); isMpc {
			if _, isRat := y.(*big.Rat); isRat {
				unsupported("divmod", complex128(0), complex128(0))
			}
		}
		unsupported("divmod", x, y)
		return false
	})
	return q, m, err
}

// Neg is -x.
func (c *Context) Neg(x Number) (Number, error) {
	return call(func() Number {
		if t, ok := x.(Mpc); ok {
			return c.mkc(libmp.MpcNeg(t.v, c.prec, c.rnd))
		}
		return c.mkf(libmp.MpfNeg(x.(Mpf).v, c.prec, c.rnd))
	})
}

// Pos is +x (x rounded to the working precision).
func (c *Context) Pos(x Number) (Number, error) {
	return call(func() Number {
		if t, ok := x.(Mpc); ok {
			return c.mkc(libmp.MpcPos(t.v, c.prec, c.rnd))
		}
		return c.mkf(libmp.MpfPos(x.(Mpf).v, c.prec, c.rnd))
	})
}

// Abs is abs(x); it is always an mpf.
func (c *Context) Abs(x Number) (Mpf, error) {
	return call(func() Mpf {
		if t, ok := x.(Mpc); ok {
			return c.mkf(libmp.MpcAbs(t.v, c.prec, c.rnd))
		}
		return c.mkf(libmp.MpfAbs(x.(Mpf).v, c.prec, c.rnd))
	})
}

// Conjugate is x.conjugate().
func (c *Context) Conjugate(x Number) (Number, error) {
	return call(func() Number {
		if t, ok := x.(Mpc); ok {
			return c.mkc(libmp.MpcConjugate(t.v, c.prec, c.rnd))
		}
		return x
	})
}

var swapped = map[string]string{"<": ">", ">": "<", "<=": ">=", ">=": "<=", "==": "==", "!=": "!="}

// mpfCmp is mpf._cmp / mpf.__eq__.
func (c *Context) mpfCmp(op string, s Mpf, other any) (bool, bool) {
	if op == "==" || op == "!=" {
		eq, ok := c.mpfEq(s, other)
		if !ok {
			return false, false
		}
		return eq == (op == "=="), true
	}
	var ov libmp.MPF
	var conv Number
	switch o := other.(type) {
	case Mpf:
		ov = o.v
	default:
		n, ok := c.tryConvert(other, true)
		if !ok {
			return false, false
		}
		conv = n
		if f, isMpf := n.(Mpf); isMpf {
			ov, conv = f.v, nil
		}
	}
	if conv != nil {
		if s.v.Eq(libmp.Fnan) {
			return false, true
		}
		libmp.Raise("TypeError", "%s", "cannot unpack non-iterable "+typeName(conv)+" object")
	}
	switch op {
	case "<":
		return libmp.MpfLt(s.v, ov), true
	case "<=":
		return libmp.MpfLe(s.v, ov), true
	case ">":
		return libmp.MpfGt(s.v, ov), true
	}
	return libmp.MpfGe(s.v, ov), true
}

func (c *Context) mpfEq(s Mpf, other any) (bool, bool) {
	switch o := other.(type) {
	case Mpf:
		return libmp.MpfEq(s.v, o.v), true
	case Mpc:
		return o.v.Im.Eq(libmp.Fzero) && libmp.MpfEq(o.v.Re, s.v), true
	}
	n, ok := c.tryConvert(other, false)
	if !ok {
		return false, false
	}
	return c.mpfEq(s, n)
}

func (c *Context) mpcEq(s Mpc, other any) (bool, bool) {
	var o Mpc
	switch t := other.(type) {
	case Mpc:
		o = t
	case string:
		return false, true
	default:
		n, ok := c.tryConvert(other, true)
		if !ok {
			return false, false
		}
		switch m := n.(type) {
		case Mpc:
			o = m
		case Mpf:
			return libmp.MpfEq(s.v.Re, m.v) && libmp.MpfEq(s.v.Im, libmp.Fzero), true
		}
	}
	return libmp.MpfEq(s.v.Re, o.v.Re) && libmp.MpfEq(s.v.Im, o.v.Im), true
}

func (c *Context) cmpForward(op string, x, y any) (bool, bool) {
	switch t := x.(type) {
	case Mpf:
		return c.mpfCmp(op, t, y)
	case Mpc:
		if op == "==" || op == "!=" {
			eq, ok := c.mpcEq(t, y)
			if !ok {
				return false, false
			}
			return eq == (op == "=="), true
		}
	}
	return false, false
}

func (c *Context) cmpop(op string, x, y any) bool {
	_, xmp := x.(Number)
	if !xmp {
		if _, ymp := y.(Number); !ymp {
			libmp.Raise("TypeError", "at least one operand must be an mpf or mpc")
		}
	}
	if r, ok := c.cmpForward(op, x, y); ok {
		return r
	}
	if r, ok := c.cmpForward(swapped[op], y, x); ok {
		return r
	}
	switch op {
	case "==":
		return false
	case "!=":
		return true
	}
	libmp.Raise("TypeError", "%s", "'"+op+"' not supported between instances of '"+typeName(x)+"' and '"+typeName(y)+"'")
	return false
}

func (c *Context) cmpErr(op string, x, y any) (bool, error) {
	return call(func() bool { return c.cmpop(op, x, y) })
}

// Eq is x == y.
func (c *Context) Eq(x, y any) (bool, error) { return c.cmpErr("==", x, y) }

// Ne is x != y.
func (c *Context) Ne(x, y any) (bool, error) { return c.cmpErr("!=", x, y) }

// Lt is x < y.
func (c *Context) Lt(x, y any) (bool, error) { return c.cmpErr("<", x, y) }

// Le is x <= y.
func (c *Context) Le(x, y any) (bool, error) { return c.cmpErr("<=", x, y) }

// Gt is x > y.
func (c *Context) Gt(x, y any) (bool, error) { return c.cmpErr(">", x, y) }

// Ge is x >= y.
func (c *Context) Ge(x, y any) (bool, error) { return c.cmpErr(">=", x, y) }

// Int is int(x).
func (c *Context) Int(x Number) (*big.Int, error) {
	return call(func() *big.Int {
		if t, ok := x.(Mpf); ok {
			return libmp.ToInt(t.v, libmp.RoundFast)
		}
		libmp.Raise("TypeError", "%s", "int() argument must be a string, a bytes-like object or a real number, not 'mpc'")
		return nil
	})
}

// Float is float(x).
func (c *Context) Float(x Number) (float64, error) {
	return call(func() float64 {
		if t, ok := x.(Mpf); ok {
			return libmp.ToFloat(t.v, false, c.rnd)
		}
		libmp.Raise("TypeError", "%s", "float() argument must be a string or a real number, not 'mpc'")
		return 0
	})
}

// Complex is complex(x).
func (c *Context) Complex(x Number) (complex128, error) {
	return call(func() complex128 {
		if t, ok := x.(Mpc); ok {
			return libmp.MpcToComplex(t.v, false, c.rnd)
		}
		return complex(libmp.ToFloat(x.(Mpf).v, false, c.rnd), 0)
	})
}

// Bool is bool(x).
func (c *Context) Bool(x Number) bool {
	if t, ok := x.(Mpc); ok {
		return libmp.MpcIsNonzero(t.v)
	}
	return !x.(Mpf).v.Eq(libmp.Fzero)
}

var hashModulus = new(big.Int).Sub(new(big.Int).Lsh(big.NewInt(1), 61), big.NewInt(1))

// pyIntHash is hash() of a Python int.
func pyIntHash(h *big.Int) int64 {
	if h.IsInt64() && h.Int64() != -1 {
		return h.Int64()
	}
	m := new(big.Int).Rem(new(big.Int).Abs(h), hashModulus)
	if h.Sign() < 0 {
		m.Neg(m)
	}
	v := m.Int64()
	if v == -1 {
		v = -2
	}
	return v
}

// Hash is hash(x); it is not defined here for nan, whose hash is its identity.
func (c *Context) Hash(x Number) (int64, error) {
	return call(func() int64 {
		if t, ok := x.(Mpc); ok {
			if t.v.Re.Eq(libmp.Fnan) || t.v.Im.Eq(libmp.Fnan) {
				libmp.Raise("ValueError", "hash of nan is not defined")
			}
			return pyIntHash(libmp.MpcHash(t.v))
		}
		v := x.(Mpf).v
		if v.Eq(libmp.Fnan) {
			libmp.Raise("ValueError", "hash of nan is not defined")
		}
		return libmp.MpfHash(v)
	})
}

// ManExp is x.man_exp: the unsigned mantissa and the exponent.
func (x Mpf) ManExp() (*big.Int, int, error) {
	type pair struct {
		m *big.Int
		e int
	}
	p, err := call(func() pair { m, e := libmp.ToManExp(x.v, false); return pair{m, e} })
	return p.m, p.e, err
}

// AsIntegerRatio is x.as_integer_ratio().
func (x Mpf) AsIntegerRatio() (*big.Int, *big.Int, error) {
	type pair struct{ p, q *big.Int }
	r, err := call(func() pair { p, q := libmp.ToRational(x.v); return pair{p, q} })
	return r.p, r.q, err
}

// ToFixed is x.to_fixed(prec).
func (x Mpf) ToFixed(prec int) *big.Int { return libmp.ToFixed(x.v, prec) }
