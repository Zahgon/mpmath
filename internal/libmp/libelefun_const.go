// Copyright (c) 2005-2026 Fredrik Johansson and mpmath contributors.
// Go port distributed under the same BSD 3-clause licence (see LICENSE).

package libmp

import (
	"math"
	"math/big"
	"sync"
)

// Constants, powers and roots from mpmath/libmp/libelefun.py.

// memo is constant_memo: the value is kept at the highest precision asked so
// far and lower precisions are truncations of it, as in the original.
type memo struct {
	mu   sync.Mutex
	prec int
	val  *big.Int
	f    func(int) *big.Int
}

func newMemo(f func(int) *big.Int) *memo { return &memo{prec: -1, f: f} }

func (m *memo) get(prec int) *big.Int {
	m.mu.Lock()
	defer m.mu.Unlock()
	if prec <= m.prec {
		return rsh(m.val, m.prec-prec)
	}
	newprec := int(float64(float64(prec)*1.05) + 10)
	m.val = m.f(newprec)
	m.prec = newprec
	return rsh(m.val, newprec-prec)
}

// defMpfConstant is def_mpf_constant.
func defMpfConstant(fixed func(int) *big.Int) func(int, Rnd) MPF {
	return func(prec int, rnd Rnd) MPF {
		wp := prec + 20
		v := fixed(wp)
		if rnd == RoundUp || rnd == RoundCeiling {
			v = addi(v, 1)
		}
		return Normalize(0, v, -wp, v.BitLen(), prec, rnd)
	}
}

func bspAcot(q *big.Int, a, b int, hyperbolic bool) (*big.Int, *big.Int, *big.Int) {
	if b-a == 1 {
		a1 := I(int64(2*a + 3))
		if hyperbolic || a&1 == 1 {
			return I(1), mul(a1, mul(q, q)), a1
		}
		return I(-1), mul(a1, mul(q, q)), a1
	}
	m := floorDivInt(a+b, 2)
	p1, q1, r1 := bspAcot(q, a, m, hyperbolic)
	p2, q2, r2 := bspAcot(q, m, b, hyperbolic)
	return add(mul(q2, p1), mul(r1, p2)), mul(q1, q2), mul(r1, r2)
}

func acotFixed(a int64, prec int, hyperbolic bool) *big.Int {
	N := int(0.35*float64(prec)/math.Log(float64(a)) + 20)
	p, q, _ := bspAcot(I(a), 0, N, hyperbolic)
	return floorDiv(lsh(add(p, q), prec), muli(q, a))
}

func machin(coefs [][2]int64, prec int, hyperbolic bool) *big.Int {
	extraprec := 10
	s := I(0)
	for _, c := range coefs {
		s = add(s, muli(acotFixed(c[1], prec+extraprec, hyperbolic), c[0]))
	}
	return rsh(s, extraprec)
}

var (
	ln2Memo = newMemo(func(prec int) *big.Int {
		return machin([][2]int64{{18, 26}, {-2, 4801}, {8, 8749}}, prec, true)
	})
	ln10Memo = newMemo(func(prec int) *big.Int {
		return machin([][2]int64{{46, 31}, {34, 49}, {20, 161}}, prec, true)
	})
	piMemo  = newMemo(piFixedRaw)
	eMemo   = newMemo(eFixedRaw)
	phiMemo = newMemo(func(prec int) *big.Int {
		prec += 10
		a := add(Isqrt(lsh(bigFive, 2*prec)), lsh(bigOne, prec))
		return rsh(a, 11)
	})
	lnSqrt2piMemo = newMemo(func(prec int) *big.Int {
		wp := prec + 10
		return ToFixed(MpfLn(MpfShift(MpfPi(wp, RoundFast), 1), wp, RoundFast), prec-1)
	})
	sqrtpiMemo = newMemo(func(prec int) *big.Int {
		return SqrtFixed(PiFixed(prec), prec)
	})
)

// Ln2Fixed is floor(ln(2) * 2**prec) through the memo.
func Ln2Fixed(prec int) *big.Int { return ln2Memo.get(prec) }

// Ln10Fixed is ln(10) as a fixed-point number.
func Ln10Fixed(prec int) *big.Int { return ln10Memo.get(prec) }

// PiFixed is floor(pi * 2**prec).
func PiFixed(prec int) *big.Int { return piMemo.get(prec) }

// EFixed is e as a fixed-point number.
func EFixed(prec int) *big.Int { return eMemo.get(prec) }

// PhiFixed is the golden ratio as a fixed-point number.
func PhiFixed(prec int) *big.Int { return phiMemo.get(prec) }

// DegreeFixed is pi/180 as a fixed-point number.
func DegreeFixed(prec int) *big.Int { return floorDiv(PiFixed(prec), I(180)) }

// LnSqrt2piFixed is ln(sqrt(2*pi)) as a fixed-point number.
func LnSqrt2piFixed(prec int) *big.Int { return lnSqrt2piMemo.get(prec) }

// SqrtpiFixed is sqrt(pi) as a fixed-point number.
func SqrtpiFixed(prec int) *big.Int { return sqrtpiMemo.get(prec) }

var (
	chudA = I(13591409)
	chudB = I(545140134)
	chudC = I(640320)
	chudD = I(12)
)

func bsChudnovsky(a, b int) (*big.Int, *big.Int, *big.Int) {
	if b-a == 1 {
		g := I(int64((6*b - 5) * (2*b - 1) * (6*b - 1)))
		bb := I(int64(b))
		p := floorDiv(mul(mul(mul(bb, bb), bb), powi(chudC, 3)), I(24))
		q := mul(g, add(chudA, mul(chudB, bb)))
		if b&1 == 1 {
			q = neg(q)
		}
		return g, p, q
	}
	mid := floorDivInt(a+b, 2)
	g1, p1, q1 := bsChudnovsky(a, mid)
	g2, p2, q2 := bsChudnovsky(mid, b)
	return mul(g1, g2), mul(p1, p2), add(mul(q1, p2), mul(q2, g1))
}

func piFixedRaw(prec int) *big.Int {
	N := int(float64(prec)/3.3219280948/14.181647462 + 2)
	_, p, q := bsChudnovsky(0, N)
	sqrtC := Isqrt(lsh(chudC, 2*prec))
	return floorDiv(mul(mul(p, chudC), sqrtC), mul(add(q, mul(chudA, p)), chudD))
}

func bspe(a, b int) (*big.Int, *big.Int) {
	if b-a == 1 {
		return I(1), I(int64(b))
	}
	m := floorDivInt(a+b, 2)
	p1, q1 := bspe(a, m)
	p2, q2 := bspe(m, b)
	return add(mul(p1, q2), p2), mul(q1, q2)
}

func eFixedRaw(prec int) *big.Int {
	N := int(1.1*float64(prec)/math.Log(float64(prec)) + 20)
	p, q := bspe(0, N)
	return floorDiv(lsh(add(p, q), prec), q)
}

// The mpf constants.
var (
	MpfPhi       = defMpfConstant(PhiFixed)
	MpfPi        = defMpfConstant(PiFixed)
	MpfE         = defMpfConstant(EFixed)
	MpfDegree    = defMpfConstant(DegreeFixed)
	MpfLn2       = defMpfConstant(Ln2Fixed)
	MpfLn10      = defMpfConstant(Ln10Fixed)
	MpfSqrtpi    = defMpfConstant(SqrtpiFixed)
	MpfLnSqrt2pi = defMpfConstant(LnSqrt2piFixed)
)

// MpfPow is s**t; a negative s with a fractional t raises ComplexResult.
func MpfPow(s, t MPF, prec int, rnd Rnd) MPF {
	if s.Sign == 1 && t.Exp < 0 {
		Raise("ComplexResult", "negative number raised to a fractional power")
	}
	if t.Exp >= 0 {
		n := lsh(t.Man, t.Exp)
		if t.Sign == 1 {
			n = neg(n)
		}
		return MpfPowInt(s, n, prec, rnd)
	}
	if t.Exp == -1 {
		if t.Man.Cmp(bigOne) == 0 {
			if t.Sign == 1 {
				return MpfDiv(Fone, MpfSqrt(s, prec+10, ReciprocalRnd(rnd)), prec, rnd)
			}
			return MpfSqrt(s, prec, rnd)
		}
		if t.Sign == 1 {
			return MpfPowInt(MpfSqrt(s, prec+10, ReciprocalRnd(rnd)), neg(t.Man), prec, rnd)
		}
		return MpfPowInt(MpfSqrt(s, prec+10, rnd), t.Man, prec, rnd)
	}
	if s.Eq(Fone) {
		return Fone
	}
	c := MpfLn(s, prec+10, rnd)
	return MpfExp(MpfMul(t, c, 0, RoundFast), prec, rnd)
}

func intPowFixed(y *big.Int, n int, prec int) (*big.Int, int) {
	if n == 2 {
		return mul(y, y), 0
	}
	bc := y.BitLen()
	exp := 0
	workprec := 2 * (prec + 4*intBitLen(n) + 4)
	pm, pe, pbc := I(1), 0, 1
	for {
		if n&1 == 1 {
			pm = mul(pm, y)
			pe = pe + exp
			pbc += bc - 2
			pbc = pbc + rsh(pm, pbc).BitLen()
			if pbc > workprec {
				pm = rsh(pm, pbc-workprec)
				pe += pbc - workprec
				pbc = workprec
			}
			n--
			if n == 0 {
				break
			}
		}
		y = mul(y, y)
		exp = exp + exp
		bc = bc + bc - 2
		bc = bc + rsh(y, bc).BitLen()
		if bc > workprec {
			y = rsh(y, bc-workprec)
			exp += bc - workprec
			bc = workprec
		}
		n = n / 2
	}
	return pm, pe
}

func nthrootFixed(y *big.Int, n int, prec int) *big.Int {
	start := 50
	y1 := rsh(y, prec-n*start)
	var r *big.Int
	f, _ := new(big.Float).SetInt(y1).Float64()
	if !math.IsInf(f, 0) {
		r, _ = new(big.Float).SetFloat64(PyPow(f, 1.0/float64(n))).Int(nil)
	} else {
		y1f := FromInt(y1, start, RoundFast)
		fn := MpfRdivInt(bigOne, FromInt64(int64(n)), start, RoundFast)
		r = ToInt(MpfPow(y1f, fn, start, RoundFast), RoundFast)
	}
	extra := 10
	extra1 := n
	prevp := start
	for _, p := range GiantSteps(start, prec+extra, 2) {
		pm, pe := intPowFixed(r, n-1, prevp)
		r2 := rsh(pm, (n-1)*prevp-p-pe-extra1)
		B := floorDiv(lsh(y, 2*p-prec+extra1), r2)
		r = floorDiv(add(B, muli(lsh(r, p-prevp), int64(n-1))), I(int64(n)))
		prevp = p
	}
	return r
}

// MpfNthroot is the real n-th root of a nonnegative value.
func MpfNthroot(s MPF, n int, prec int, rnd Rnd) MPF {
	if s.Sign == 1 {
		Raise("ComplexResult", "nth root of a negative number")
	}
	if s.Man.Sign() == 0 {
		if s.Eq(Fnan) {
			return Fnan
		}
		if s.Eq(Fzero) {
			if n > 0 {
				return Fzero
			}
			if n == 0 {
				return Fone
			}
			return Finf
		}
		if n == 0 {
			return Fnan
		}
		if n < 0 {
			return Fzero
		}
		return Finf
	}
	flagInverse := false
	extraInverse := 0
	if n < 2 {
		if n == 0 {
			return Fone
		}
		if n == 1 {
			return MpfPos(s, prec, rnd)
		}
		if n == -1 {
			return MpfDiv(Fone, s, prec, rnd)
		}
		rnd = ReciprocalRnd(rnd)
		flagInverse = true
		extraInverse = 5
		prec += extraInverse
		n = -n
	}
	if n > 20 && (n >= 20000 || prec < int(233+28.3*PyPow(float64(n), 0.62))) {
		prec2 := prec + 10
		nth := MpfRdivInt(bigOne, FromInt64(int64(n)), prec2, RoundFast)
		r := MpfPow(s, nth, prec2, rnd)
		s = Normalize(r.Sign, r.Man, r.Exp, r.BC, prec, rnd)
		if flagInverse {
			return MpfDiv(Fone, s, prec-extraInverse, rnd)
		}
		return s
	}
	man, exp, bc := s.Man, s.Exp, s.BC
	prec2 := prec + 2*n - pyModInt(prec, n)
	if n > 10 {
		prec2 += prec2 / 10
		prec2 = prec2 - pyModInt(prec2, n)
	}
	shift := bc - prec2
	sign1 := false
	es := exp + shift
	if es < 0 {
		sign1 = true
		es = -es
	}
	if sign1 {
		shift += es % n
	} else {
		shift -= es % n
	}
	man = rsh(man, shift)
	extra := 10
	exp1 := floorDivInt(exp+shift-(n-1)*prec2, n) - extra
	rndShift := int64(0)
	if flagInverse {
		if rnd == RoundUp || rnd == RoundCeiling {
			rndShift = 1
		}
	} else if rnd == RoundDown || rnd == RoundFloor {
		rndShift = 1
	}
	man = nthrootFixed(addi(man, rndShift), n, prec2)
	s = FromManExp(man, exp1, prec, rnd)
	if flagInverse {
		return MpfDiv(Fone, s, prec-extraInverse, rnd)
	}
	return s
}

// MpfCbrt is the cube root.
func MpfCbrt(s MPF, prec int, rnd Rnd) MPF { return MpfNthroot(s, 3, prec, rnd) }
