// Copyright (c) 2005-2026 Fredrik Johansson and mpmath contributors.
// Go port distributed under the same BSD 3-clause licence (see LICENSE).

package libmp

import (
	"math"
	"math/big"
	"sync"
)

// Logarithms and inverse trigonometric functions from libelefun.py.

const (
	logTaylorPrec      = 2500
	logTaylorShift     = 9
	logAgmMagPrecRatio = 20
	atanTaylorPrec     = 3000
	atanTaylorShift    = 7
	maxLogIntCache     = 2000
)

var cachePrecSteps = func() []int {
	steps := []int{22, 22}
	for k := 1; k <= intBitLen(logTaylorPrec); k++ {
		v := minInt(1<<k, logTaylorPrec) + 20
		for i := 0; i < 1<<(k-1); i++ {
			steps = append(steps, v)
		}
	}
	return steps
}()

type pairKey struct {
	n    string
	prec int
}

type fixedPair struct{ a, b *big.Int }

var (
	cacheMu         sync.Mutex
	logTaylorCache  = map[pairKey]fixedPair{}
	atanTaylorCache = map[pairKey]fixedPair{}
	logIntCache     = map[int]fixedPair{}
)

// LogIntFixed is log(n) as a fixed-point number, cached for small n.
func LogIntFixed(n int, prec int) *big.Int {
	cacheMu.Lock()
	c, ok := logIntCache[n]
	cacheMu.Unlock()
	if ok {
		vprec := int(c.b.Int64())
		if vprec >= prec {
			return rsh(c.a, vprec-prec)
		}
	}
	wp := prec + 10
	v := ToFixed(MpfLn(FromInt64(int64(n)), wp+5, RoundFast), wp)
	if n < maxLogIntCache {
		cacheMu.Lock()
		logIntCache[n] = fixedPair{v, I(int64(wp))}
		cacheMu.Unlock()
	}
	return rsh(v, wp-prec)
}

// AgmFixed is agm(a, b) for fixed-point a, b near unit magnitude.
func AgmFixed(a, b *big.Int, prec int) *big.Int {
	i := 0
	for {
		anew := rsh(add(a, b), 1)
		if i > 4 && absI(sub(a, anew)).Cmp(I(8)) < 0 {
			return a
		}
		b = Isqrt(mul(a, b))
		a = anew
		i++
	}
}

func logAgm(x *big.Int, prec int) *big.Int {
	x2 := rsh(mul(x, x), prec)
	s, a, b := x2, x2, x2
	for a.Sign() != 0 {
		b = rsh(mul(b, x2), prec)
		a = rsh(mul(a, b), prec)
		s = add(s, a)
	}
	s = add(s, lsh(bigOne, prec))
	s = rsh(mul(s, s), prec-2)
	s = rsh(mul(s, Isqrt(lsh(x, prec))), prec)
	t, a, b := x, x, x
	for a.Sign() != 0 {
		b = rsh(mul(b, x2), prec)
		a = rsh(mul(a, b), prec)
		t = add(t, a)
	}
	t = add(lsh(bigOne, prec), lsh(t, 1))
	t = rsh(mul(t, t), prec)
	p := AgmFixed(s, t, prec)
	return floorDiv(lsh(PiFixed(prec), prec), p)
}

func logTaylor(x *big.Int, prec, r int) *big.Int {
	for i := 0; i < r; i++ {
		x = Isqrt(lsh(x, prec))
	}
	one := lsh(bigOne, prec)
	v := floorDiv(lsh(sub(x, one), prec), add(x, one))
	sign := v.Sign() < 0
	if sign {
		v = neg(v)
	}
	v2 := rsh(mul(v, v), prec)
	v4 := rsh(mul(v2, v2), prec)
	s0 := v
	s1 := floorDiv(v, bigThree)
	v = rsh(mul(v, v4), prec)
	k := int64(5)
	for v.Sign() != 0 {
		s0 = add(s0, floorDiv(v, I(k)))
		k += 2
		s1 = add(s1, floorDiv(v, I(k)))
		v = rsh(mul(v, v4), prec)
		k += 2
	}
	s1 = rsh(mul(s1, v2), prec)
	s := lsh(add(s0, s1), 1+r)
	if sign {
		return neg(s)
	}
	return s
}

func logTaylorCached(x *big.Int, prec int) *big.Int {
	n := rsh(x, prec-logTaylorShift)
	cachedPrec := cachePrecSteps[prec]
	dprec := cachedPrec - prec
	key := pairKey{n.String(), cachedPrec}
	cacheMu.Lock()
	c, ok := logTaylorCache[key]
	cacheMu.Unlock()
	if !ok {
		a := lsh(n, cachedPrec-logTaylorShift)
		c = fixedPair{a, logTaylor(a, cachedPrec, 8)}
		cacheMu.Lock()
		logTaylorCache[key] = c
		cacheMu.Unlock()
	}
	a := rsh(c.a, dprec)
	logA := rsh(c.b, dprec)
	u := floorDiv(lsh(sub(x, a), prec), a)
	v := floorDiv(lsh(u, prec), add(lsh(bigTwo, prec), u))
	v2 := rsh(mul(v, v), prec)
	v4 := rsh(mul(v2, v2), prec)
	s0 := v
	s1 := floorDiv(v, bigThree)
	v = rsh(mul(v, v4), prec)
	k := int64(5)
	for v.Sign() != 0 {
		s0 = add(s0, floorDiv(v, I(k)))
		k += 2
		s1 = add(s1, floorDiv(v, I(k)))
		v = rsh(mul(v, v4), prec)
		k += 2
	}
	s1 = rsh(mul(s1, v2), prec)
	s := lsh(add(s0, s1), 1)
	return add(logA, s)
}

// MpfLn is the natural logarithm; a negative x raises ComplexResult.
func MpfLn(x MPF, prec int, rnd Rnd) MPF {
	sign, man, exp, bc := x.Sign, x.Man, x.Exp, x.BC
	if man.Sign() == 0 {
		if x.Eq(Fzero) {
			return Fninf
		}
		if x.Eq(Finf) {
			return Finf
		}
		if x.Eq(Fnan) {
			return Fnan
		}
	}
	if sign == 1 {
		Raise("ComplexResult", "logarithm of a negative number")
	}
	wp := prec + 20
	if man.Cmp(bigOne) == 0 {
		if exp == 0 {
			return Fzero
		}
		return FromManExp(muli(Ln2Fixed(wp), int64(exp)), -wp, prec, rnd)
	}
	mag := exp + bc
	absMag := absInt(mag)
	if absMag <= 1 {
		tsign := 1 - absMag
		var tman *big.Int
		if tsign != 0 {
			tman = sub(lsh(bigOne, bc), man)
		} else {
			tman = sub(man, lsh(bigOne, bc-1))
		}
		tbc := tman.BitLen()
		cancellation := bc - tbc
		if cancellation > wp {
			t := Normalize(tsign, tman, absMag-bc, tbc, tbc, RoundNearest)
			return MpfPerturb(t, tsign, prec, rnd)
		}
		wp += cancellation
		wpb := intBitLen(wp)
		if wpb <= cancellation {
			a := ToFixed(x, wp)
			s := logTaylor(a, wp, 0)
			return FromManExp(s, -wp, prec, rnd)
		}
	}
	if absMag > 10000 {
		if intBitLen(absMag) > wp {
			return FromManExp(muli(Ln2Fixed(wp), int64(exp)), -wp, prec, rnd)
		}
	}
	var m *big.Int
	if wp <= logTaylorPrec {
		m = logTaylorCached(lsh(man, wp-bc), wp)
		if mag != 0 {
			m = add(m, muli(Ln2Fixed(wp), int64(mag)))
		}
	} else {
		optimalMag := floorDivInt(-wp, logAgmMagPrecRatio)
		n := optimalMag - mag
		x = MpfShift(x, n)
		wp += -optimalMag
		m = neg(logAgm(ToFixed(x, wp), wp))
		m = sub(m, muli(Ln2Fixed(wp), int64(n)))
	}
	return FromManExp(m, -wp, prec, rnd)
}

// MpfLog1p is log(1+x).
func MpfLog1p(x MPF, prec int, rnd Rnd) MPF {
	wp := prec + 20
	wp2 := wp * 2
	if x.Exp+x.BC < -wp && (x.Man.Sign() != 0 || x.Exp != 0) {
		x2 := MpfSub(Fone, MpfShift(x, -1), wp2, rnd)
		return MpfMul(x, x2, wp, rnd)
	}
	return MpfLn(MpfAdd(Fone, x, wp2, RoundFast), wp, rnd)
}

// MpfLogHypot is log(sqrt(a**2 + b**2)).
func MpfLogHypot(a, b MPF, prec int, rnd Rnd) MPF {
	if b.Man.Sign() == 0 {
		a, b = b, a
	}
	if a.Man.Sign() == 0 {
		if b.Man.Sign() == 0 {
			if a.Eq(Fzero) && b.Eq(Fzero) {
				return Fninf
			}
			if a.Eq(Fnan) || b.Eq(Fnan) {
				return Fnan
			}
			return Finf
		}
		if a.Eq(Fzero) {
			return MpfLn(MpfAbs(b, 0, RoundFast), prec, rnd)
		}
		if a.Eq(Fnan) {
			return Fnan
		}
		return Finf
	}
	a2 := MpfMul(a, a, 0, RoundFast)
	b2 := MpfMul(b, b, 0, RoundFast)
	extra := 20
	h2 := MpfAdd(a2, b2, prec+extra, RoundFast)
	cancelled := MpfAdd(h2, Fnone, 10, RoundFast)
	magCancelled := cancelled.Exp + cancelled.BC
	if cancelled.Eq(Fzero) || magCancelled < floorDivInt(-extra, 2) {
		h2 = MpfAdd(a2, b2, prec+extra-minInt(a2.Exp, b2.Exp), RoundFast)
	}
	return MpfShift(MpfLn(h2, prec, rnd), -1)
}

func atanNewton(x *big.Int, prec int) *big.Int {
	var r float64
	if prec >= 100 {
		f, _ := new(big.Float).SetInt(rsh(x, prec-53)).Float64()
		r = PyAtan(f / math.Exp2(53))
	} else {
		f, _ := new(big.Float).SetInt(x).Float64()
		r = PyAtan(f / math.Exp2(float64(prec)))
	}
	prevp := 50
	ri, _ := new(big.Float).SetFloat64(r * math.Exp2(53)).Int(nil)
	rr := rsh(ri, 53-prevp)
	extraP := 50
	for _, wp := range GiantSteps(prevp, prec, 2) {
		wp += extraP
		rr = lsh(rr, wp-prevp)
		cos, sin := CosSinFixed(rr, wp)
		tan := floorDiv(lsh(sin, wp), cos)
		a := floorDiv(lsh(sub(tan, rsh(x, prec-wp)), wp), add(lsh(bigOne, wp), rsh(mul(tan, tan), wp)))
		rr = sub(rr, a)
		prevp = wp
	}
	return rsh(rr, prevp-prec)
}

func atanTaylorGetCached(n *big.Int, prec int) (*big.Int, *big.Int) {
	prec2 := (1 << intBitLen(prec-1)) + 20
	dprec := prec2 - prec
	key := pairKey{n.String(), prec2}
	cacheMu.Lock()
	c, ok := atanTaylorCache[key]
	cacheMu.Unlock()
	if !ok {
		a := lsh(n, prec2-atanTaylorShift)
		c = fixedPair{a, atanNewton(a, prec2)}
		cacheMu.Lock()
		atanTaylorCache[key] = c
		cacheMu.Unlock()
	}
	return rsh(c.a, dprec), rsh(c.b, dprec)
}

func atanTaylor(x *big.Int, prec int) *big.Int {
	n := rsh(x, prec-atanTaylorShift)
	a, atanA := atanTaylorGetCached(n, prec)
	d := sub(x, a)
	den := add(add(rsh(mul(a, a), prec), rsh(mul(a, d), prec)), lsh(bigOne, prec))
	v := floorDiv(lsh(d, prec), den)
	s0 := v
	v2 := rsh(mul(v, v), prec)
	v4 := rsh(mul(v2, v2), prec)
	s1 := floorDiv(v, bigThree)
	v = rsh(mul(v, v4), prec)
	k := int64(5)
	for v.Sign() != 0 {
		s0 = add(s0, floorDiv(v, I(k)))
		k += 2
		s1 = add(s1, floorDiv(v, I(k)))
		v = rsh(mul(v, v4), prec)
		k += 2
	}
	s1 = rsh(mul(s1, v2), prec)
	return add(atanA, sub(s0, s1))
}

func atanInf(sign int, prec int, rnd Rnd) MPF {
	if sign == 0 {
		return MpfShift(MpfPi(prec, rnd), -1)
	}
	return MpfNeg(MpfShift(MpfPi(prec, NegativeRnd(rnd)), -1), 0, RoundFast)
}

// MpfAtan is atan(x).
func MpfAtan(x MPF, prec int, rnd Rnd) MPF {
	sign, man, exp, bc := x.Sign, x.Man, x.Exp, x.BC
	if man.Sign() == 0 {
		if x.Eq(Fzero) {
			return Fzero
		}
		if x.Eq(Finf) {
			return atanInf(0, prec, rnd)
		}
		if x.Eq(Fninf) {
			return atanInf(1, prec, rnd)
		}
		return Fnan
	}
	mag := exp + bc
	if mag > prec+20 {
		return atanInf(sign, prec, rnd)
	}
	if -mag > prec+20 {
		return MpfPerturb(x, 1-sign, prec, rnd)
	}
	wp := prec + 30 + absInt(mag)
	reciprocal := false
	if mag >= 2 {
		x = MpfRdivInt(bigOne, x, wp, RoundFast)
		reciprocal = true
	}
	t := ToFixed(x, wp)
	if sign == 1 {
		t = neg(t)
	}
	var a *big.Int
	if wp < atanTaylorPrec {
		a = atanTaylor(t, wp)
	} else {
		a = atanNewton(t, wp)
	}
	if reciprocal {
		a = sub(addi(rsh(PiFixed(wp), 1), 1), a)
	}
	if sign == 1 {
		a = neg(a)
	}
	return FromManExp(a, -wp, prec, rnd)
}

// MpfAtan2 is atan2(y, x).
func MpfAtan2(y, x MPF, prec int, rnd Rnd) MPF {
	if y.Man.Sign() == 0 {
		if y.Eq(Fzero) && !x.Eq(Fnan) {
			if MpfSign(x) >= 0 {
				return Fzero
			}
			return MpfPi(prec, rnd)
		}
		if y.Eq(Finf) || y.Eq(Fninf) {
			if x.Eq(Finf) {
				if y.Eq(Finf) {
					return MpfShift(MpfPi(prec, rnd), -2)
				}
				rnd = NegativeRnd(rnd)
				return MpfNeg(MpfShift(MpfPi(prec, rnd), -2), 0, RoundFast)
			}
			if x.Eq(Fninf) {
				if y.Eq(Finf) {
					return MpfShift(MpfMulInt(MpfPi(prec, rnd), I(3), prec, rnd), -2)
				}
				rnd = NegativeRnd(rnd)
				return MpfNeg(MpfShift(MpfMulInt(MpfPi(prec, rnd), I(3), prec, rnd), -2), 0, RoundFast)
			}
			if y.Eq(Finf) {
				return MpfShift(MpfPi(prec, rnd), -1)
			}
			return MpfNeg(MpfShift(MpfPi(prec, NegativeRnd(rnd)), -1), 0, RoundFast)
		}
		return Fnan
	}
	if y.Sign == 1 {
		return MpfNeg(MpfAtan2(MpfNeg(y, 0, RoundFast), x, prec, NegativeRnd(rnd)), 0, RoundFast)
	}
	if x.Man.Sign() == 0 {
		if x.Eq(Fnan) {
			return Fnan
		}
		if x.Eq(Finf) {
			return Fzero
		}
		if x.Eq(Fninf) {
			return MpfPi(prec, rnd)
		}
		if y.Eq(Fzero) {
			return Fzero
		}
		return MpfShift(MpfPi(prec, rnd), -1)
	}
	tquo := MpfAtan(MpfDiv(y, x, prec+4, RoundFast), prec+4, RoundFast)
	if x.Sign == 1 {
		return MpfAdd(MpfPi(prec+4, RoundFast), tquo, prec, rnd)
	}
	return MpfPos(tquo, prec, rnd)
}

// MpfAsin is asin(x).
func MpfAsin(x MPF, prec int, rnd Rnd) MPF {
	if x.BC+x.Exp > 0 && !x.Eq(Fone) && !x.Eq(Fnone) {
		Raise("ComplexResult", "asin(x) is real only for -1 <= x <= 1")
	}
	wp := prec + 15
	a := MpfMul(x, x, 0, RoundFast)
	b := MpfAdd(Fone, MpfSqrt(MpfSub(Fone, a, wp, RoundFast), wp, RoundFast), wp, RoundFast)
	c := MpfDiv(x, b, wp, RoundFast)
	return MpfShift(MpfAtan(c, prec, rnd), 1)
}

// MpfAcos is acos(x).
func MpfAcos(x MPF, prec int, rnd Rnd) MPF {
	if x.BC+x.Exp > 0 {
		if !x.Eq(Fone) && !x.Eq(Fnone) {
			Raise("ComplexResult", "acos(x) is real only for -1 <= x <= 1")
		}
		if x.Eq(Fnone) {
			return MpfPi(prec, rnd)
		}
	}
	wp := prec + 15
	a := MpfMul(x, x, 0, RoundFast)
	b := MpfSqrt(MpfSub(Fone, a, wp, RoundFast), wp, RoundFast)
	c := MpfDiv(b, MpfAdd(Fone, x, wp, RoundFast), wp, RoundFast)
	return MpfShift(MpfAtan(c, prec, rnd), 1)
}

// MpfAsinh is asinh(x).
func MpfAsinh(x MPF, prec int, rnd Rnd) MPF {
	wp := prec + 20
	mag := x.Exp + x.BC
	if mag < -8 {
		if mag < -wp {
			return MpfPerturb(x, 1-x.Sign, prec, rnd)
		}
		wp += -mag
	}
	q := MpfSqrt(MpfAdd(MpfMul(x, x, 0, RoundFast), Fone, wp, RoundFast), wp, RoundFast)
	q = MpfAdd(MpfAbs(x, 0, RoundFast), q, wp, RoundFast)
	if x.Sign == 1 {
		return MpfNeg(MpfLn(q, prec, NegativeRnd(rnd)), 0, RoundFast)
	}
	return MpfLn(q, prec, rnd)
}

// MpfAcosh is acosh(x).
func MpfAcosh(x MPF, prec int, rnd Rnd) MPF {
	wp := prec + 15
	if MpfCmp(x, Fone) == -1 {
		Raise("ComplexResult", "acosh(x) is real only for x >= 1")
	}
	q := MpfSqrt(MpfAdd(MpfMul(x, x, 0, RoundFast), Fnone, wp, RoundFast), wp, RoundFast)
	return MpfLn(MpfAdd(x, q, wp, RoundFast), prec, rnd)
}

// MpfAtanh is atanh(x).
func MpfAtanh(x MPF, prec int, rnd Rnd) MPF {
	sign, man, exp, bc := x.Sign, x.Man, x.Exp, x.BC
	if man.Sign() == 0 && exp != 0 {
		if x.Eq(Fzero) || x.Eq(Fnan) {
			return x
		}
		Raise("ComplexResult", "atanh(x) is real only for -1 <= x <= 1")
	}
	mag := bc + exp
	if mag > 0 {
		if mag == 1 && man.Cmp(bigOne) == 0 {
			if sign == 1 {
				return Fninf
			}
			return Finf
		}
		Raise("ComplexResult", "atanh(x) is real only for -1 <= x <= 1")
	}
	wp := prec + 15
	if mag < -8 {
		if mag < -wp {
			return MpfPerturb(x, sign, prec, rnd)
		}
		wp += -mag
	}
	a := MpfAdd(x, Fone, wp, RoundFast)
	b := MpfSub(Fone, x, wp, RoundFast)
	return MpfShift(MpfLn(MpfDiv(a, b, wp, RoundFast), prec, rnd), -1)
}

// MpfFibonacci is the Fibonacci function of a real x.
func MpfFibonacci(x MPF, prec int, rnd Rnd) MPF {
	if x.Man.Sign() == 0 {
		if x.Eq(Fninf) {
			return Fnan
		}
		return x
	}
	size := absInt(x.Exp + x.BC)
	if x.Exp >= 0 && x.Sign == 0 {
		if size < 10 || size <= intBitLen(prec) {
			return FromInt(Ifib(int(ToInt(x, RoundFast).Int64())), prec, rnd)
		}
	}
	wp := prec + size + 20
	a := MpfPhi(wp, RoundFast)
	b := MpfAdd(MpfShift(a, 1), Fnone, wp, RoundFast)
	u := MpfPow(a, x, wp, RoundFast)
	v := MpfCosPi(x, wp, RoundFast)
	v = MpfDiv(v, u, wp, RoundFast)
	u = MpfSub(u, v, wp, RoundFast)
	return MpfDiv(u, b, prec, rnd)
}
