// Copyright (c) 2005-2026 Fredrik Johansson and mpmath contributors.
// Go port distributed under the same BSD 3-clause licence (see LICENSE).

package libmp

import (
	"math"
	"math/big"
)

// Exponential, hyperbolic and trigonometric functions from libelefun.py,
// with the pure-Python backend's cutoffs.

const (
	expCoshCutoff    = 600
	expSeriesUCutoff = 1500
	cosSinCachePrec  = 400
	cosSinCacheStep  = 8
)

var cosSinCache = map[string]fixedPair{}

// exponentialSeries is exponential_series: type 0 gives exp(x), 1 gives
// (cosh, sinh) and 2 gives (cos, sin).
func exponentialSeries(x *big.Int, prec int, typ int) (*big.Int, *big.Int) {
	sign := 0
	if x.Sign() < 0 {
		x = neg(x)
		sign = 1
	}
	r := int(0.5 * math.Sqrt(float64(prec)))
	xmag := x.BitLen() - prec
	r = maxInt(0, xmag+r)
	extra := 10 + 2*maxInt(r, -xmag)
	wp := prec + extra
	x = lsh(x, extra-r)
	one := lsh(bigOne, wp)
	alt := typ == 2
	var c *big.Int
	if prec < expSeriesUCutoff {
		x2 := rsh(mul(x, x), wp)
		a := x2
		x4 := rsh(mul(x2, x2), wp)
		s0, s1 := I(0), I(0)
		k := int64(2)
		for a.Sign() != 0 {
			a = floorDiv(a, I((k-1)*k))
			s0 = add(s0, a)
			k += 2
			a = floorDiv(a, I((k-1)*k))
			s1 = add(s1, a)
			k += 2
			a = rsh(mul(a, x4), wp)
		}
		s1 = rsh(mul(x2, s1), wp)
		if alt {
			c = add(sub(s1, s0), one)
		} else {
			c = add(add(s1, s0), one)
		}
	} else {
		u := int(0.3 * PyPow(float64(prec), 0.35))
		x2 := rsh(mul(x, x), wp)
		a := x2
		xpowers := []*big.Int{one, x2}
		for i := 1; i < u; i++ {
			xpowers = append(xpowers, rsh(mul(xpowers[len(xpowers)-1], x2), wp))
		}
		sums := make([]*big.Int, u)
		for i := range sums {
			sums[i] = I(0)
		}
		k := int64(2)
		for a.Sign() != 0 {
			for i := 0; i < u; i++ {
				a = floorDiv(a, I((k-1)*k))
				if alt && k&2 != 0 {
					sums[i] = sub(sums[i], a)
				} else {
					sums[i] = add(sums[i], a)
				}
				k += 2
			}
			a = rsh(mul(a, xpowers[len(xpowers)-1]), wp)
		}
		for i := 1; i < u; i++ {
			sums[i] = rsh(mul(sums[i], xpowers[i]), wp)
		}
		c = I(0)
		for _, s := range sums {
			c = add(c, s)
		}
		c = add(c, one)
	}
	if typ == 0 {
		s := Isqrt(sub(mul(c, c), lsh(one, wp)))
		var v *big.Int
		if sign == 1 {
			v = sub(c, s)
		} else {
			v = add(c, s)
		}
		for i := 0; i < r; i++ {
			v = rsh(mul(v, v), wp)
		}
		return rsh(v, extra), nil
	}
	pshift := wp - 1
	for i := 0; i < r; i++ {
		c = sub(rsh(mul(c, c), pshift), one)
	}
	s := Isqrt(absI(sub(lsh(one, wp), mul(c, c))))
	if sign == 1 {
		s = neg(s)
	}
	return rsh(c, extra), rsh(s, extra)
}

// ExpBasecase is exp(x) for a fixed-point x.
func ExpBasecase(x *big.Int, prec int) *big.Int {
	if prec > expCoshCutoff {
		v, _ := exponentialSeries(x, prec, 0)
		return v
	}
	r := int(math.Sqrt(float64(prec)))
	prec += r
	s0 := lsh(bigOne, prec)
	s1 := s0
	k := int64(2)
	x2 := rsh(mul(x, x), prec)
	a := x2
	for a.Sign() != 0 {
		a = floorDiv(a, I(k))
		s0 = add(s0, a)
		k++
		a = floorDiv(a, I(k))
		s1 = add(s1, a)
		k++
		a = rsh(mul(a, x2), prec)
	}
	s1 = rsh(mul(s1, x), prec)
	s := add(s0, s1)
	u := r
	for r != 0 {
		s = rsh(mul(s, s), prec)
		r--
	}
	return rsh(s, u)
}

func expExpnegBasecase(x *big.Int, prec int) (*big.Int, *big.Int) {
	if prec > expCoshCutoff {
		cosh, sinh := exponentialSeries(x, prec, 1)
		return add(cosh, sinh), sub(cosh, sinh)
	}
	a := ExpBasecase(x, prec)
	b := floorDiv(lsh(bigOne, prec+prec), a)
	return a, b
}

// CosSinBasecase is cos(x), sin(x) for a fixed-point x in [0, pi/2).
func CosSinBasecase(x *big.Int, prec int) (*big.Int, *big.Int) {
	if prec > cosSinCachePrec {
		return exponentialSeries(x, prec, 2)
	}
	precs := prec - cosSinCacheStep
	t := rsh(x, precs)
	key := t.String()
	cacheMu.Lock()
	c, ok := cosSinCache[key]
	cacheMu.Unlock()
	if !ok {
		w := lsh(t, 10+cosSinCachePrec-cosSinCacheStep)
		cosT, sinT := exponentialSeries(w, 10+cosSinCachePrec, 2)
		c = fixedPair{rsh(cosT, 10), rsh(sinT, 10)}
		cacheMu.Lock()
		cosSinCache[key] = c
		cacheMu.Unlock()
	}
	offset := cosSinCachePrec - prec
	cosT := rsh(c.a, offset)
	sinT := rsh(c.b, offset)
	x = sub(x, lsh(t, precs))
	cos := lsh(bigOne, prec)
	sin := x
	k := int64(2)
	a := neg(rsh(mul(x, x), prec))
	for a.Sign() != 0 {
		a = floorDiv(a, I(k))
		cos = add(cos, a)
		k++
		a = rsh(mul(a, x), prec)
		a = floorDiv(a, I(k))
		sin = add(sin, a)
		k++
		a = neg(rsh(mul(a, x), prec))
	}
	return rsh(sub(mul(cos, cosT), mul(sin, sinT)), prec), rsh(add(mul(sin, cosT), mul(cos, sinT)), prec)
}

// MpfExp is exp(x).
func MpfExp(x MPF, prec int, rnd Rnd) MPF {
	sign, man, exp, bc := x.Sign, x.Man, x.Exp, x.BC
	if man.Sign() != 0 {
		mag := bc + exp
		wp := prec + 14
		if sign == 1 {
			man = neg(man)
		}
		if prec > 600 && exp >= 0 {
			e := MpfE(wp+int(1.45*float64(mag)), RoundFast)
			return MpfPowInt(e, lsh(man, exp), prec, rnd)
		}
		if mag < -wp {
			return MpfPerturb(Fone, sign, prec, rnd)
		}
		var t *big.Int
		n := 0
		if mag > 1 {
			wpmod := wp + mag
			t = lsh(man, exp+wpmod)
			lg2 := Ln2Fixed(wpmod)
			var nn *big.Int
			nn, t = divmod(t, lg2)
			n = int(nn.Int64())
			t = rsh(t, mag)
		} else {
			t = lsh(man, exp+wp)
		}
		man = ExpBasecase(t, wp)
		return FromManExp(man, n-wp, prec, rnd)
	}
	if exp == 0 {
		return Fone
	}
	if x.Eq(Fninf) {
		return Fzero
	}
	return x
}

// MpfCoshSinh returns cosh(x) and sinh(x).
func MpfCoshSinh(x MPF, prec int, rnd Rnd) (MPF, MPF) {
	c, s, _ := mpfCoshSinh(x, prec, rnd, false)
	return c, s
}

// MpfTanh is tanh(x).
func MpfTanh(x MPF, prec int, rnd Rnd) MPF {
	_, _, t := mpfCoshSinh(x, prec, rnd, true)
	return t
}

func mpfCoshSinh(x MPF, prec int, rnd Rnd, tanh bool) (MPF, MPF, MPF) {
	sign, man, exp, bc := x.Sign, x.Man, x.Exp, x.BC
	if x.IsSpecial() {
		if tanh {
			if x.Eq(Finf) {
				return MPF{}, MPF{}, Fone
			}
			if x.Eq(Fninf) {
				return MPF{}, MPF{}, Fnone
			}
			return MPF{}, MPF{}, Fnan
		}
		if x.Eq(Finf) {
			return Finf, Finf, MPF{}
		}
		if x.Eq(Fninf) {
			return Finf, Fninf, MPF{}
		}
		return Fnan, Fnan, MPF{}
	}
	mag := exp + bc
	wp := prec + 14
	if mag < -4 {
		if mag < -wp {
			if tanh {
				return MPF{}, MPF{}, MpfPerturb(x, 1-sign, prec, rnd)
			}
			return MpfPerturb(Fone, 0, prec, rnd), MpfPerturb(x, sign, prec, rnd), MPF{}
		}
		wp += -mag
	}
	if mag > 10 {
		if 3*(1<<uint(minInt(mag-1, 62))) > wp || mag-1 >= 62 {
			if tanh {
				one := Fone
				if sign == 1 {
					one = Fnone
				}
				return MPF{}, MPF{}, MpfPerturb(one, 1-sign, prec, rnd)
			}
			c := MpfShift(MpfExp(MpfAbs(x, 0, RoundFast), prec, rnd), -1)
			s := c
			if sign == 1 {
				s = MpfNeg(s, 0, RoundFast)
			}
			return c, s, MPF{}
		}
	}
	var t *big.Int
	n := 0
	if mag > 1 {
		wpmod := wp + mag
		t = lsh(man, exp+wpmod)
		lg2 := Ln2Fixed(wpmod)
		var nn *big.Int
		nn, t = divmod(t, lg2)
		n = int(nn.Int64())
		t = rsh(t, mag)
	} else {
		t = lsh(man, exp+wp)
	}
	a, b := expExpnegBasecase(t, wp)
	cosh := add(a, rsh(b, 2*n))
	sinh := sub(a, rsh(b, 2*n))
	if sign == 1 {
		sinh = neg(sinh)
	}
	if tanh {
		m := floorDiv(lsh(sinh, wp), cosh)
		return MPF{}, MPF{}, FromManExp(m, -wp, prec, rnd)
	}
	return FromManExp(cosh, n-wp-1, prec, rnd), FromManExp(sinh, n-wp-1, prec, rnd), MPF{}
}

// MpfCosh is cosh(x).
func MpfCosh(x MPF, prec int, rnd Rnd) MPF { c, _ := MpfCoshSinh(x, prec, rnd); return c }

// MpfSinh is sinh(x).
func MpfSinh(x MPF, prec int, rnd Rnd) MPF { _, s := MpfCoshSinh(x, prec, rnd); return s }

// low2 is n & 3 with Python's two's-complement semantics.
func low2(n *big.Int) int { return int(new(big.Int).And(n, bigThree).Int64()) }

func modPi2(man *big.Int, exp, mag, wp int) (*big.Int, int, int) {
	if mag > 0 {
		i := 0
		for {
			cancellationPrec := 20 << i
			wpmod := wp + mag + cancellationPrec
			pi2 := PiFixed(wpmod - 1)
			pi4 := rsh(pi2, 1)
			t := lsh(man, wpmod+exp)
			n, y := divmod(t, pi2)
			var small *big.Int
			if y.Cmp(pi4) > 0 {
				small = sub(pi2, y)
			} else {
				small = y
			}
			if rsh(small, wp+mag-10).Sign() != 0 {
				return rsh(y, mag), low2(n), wpmod - mag
			}
			i++
		}
	}
	wp += -mag
	return lsh(man, exp+wp), 0, wp
}

// mpfCosSin computes cos and/or sin (of pi*x when pi is set); which is 0 for
// both, 1 cos, 2 sin and 3 tan (returned as the first value).
func mpfCosSin(x MPF, prec int, rnd Rnd, which int, pi bool) (MPF, MPF) {
	sign, man, exp, bc := x.Sign, x.Man, x.Exp, x.BC
	pick := func(c, s MPF) (MPF, MPF) {
		switch which {
		case 1:
			return c, MPF{}
		case 2, 3:
			return s, MPF{}
		}
		return c, s
	}
	if man.Sign() == 0 {
		if exp != 0 {
			return pick(Fnan, Fnan)
		}
		return pick(Fone, Fzero)
	}
	mag := bc + exp
	wp := prec + 10
	if mag < 0 && mag < -wp {
		if pi {
			x = MpfMul(x, MpfPi(wp, RoundFast), 0, RoundFast)
		}
		c := MpfPerturb(Fone, 1, prec, rnd)
		s := MpfPerturb(x, 1-sign, prec, rnd)
		if which == 3 {
			return MpfPerturb(x, sign, prec, rnd), MPF{}
		}
		return pick(c, s)
	}
	var t *big.Int
	var n int
	if pi {
		if exp >= -1 {
			var c, s MPF
			if exp == -1 {
				c = Fzero
				if (man.Bit(1) == 1) != (sign == 1) {
					s = Fnone
				} else {
					s = Fone
				}
			} else if exp == 0 {
				c, s = Fnone, Fzero
			} else {
				c, s = Fone, Fzero
			}
			if which == 3 {
				return MpfDiv(s, c, prec, rnd), MPF{}
			}
			return pick(c, s)
		}
		nn := rsh(addi(rsh(man, -exp-2), 1), 1)
		man = sub(man, lsh(nn, -exp-1))
		mag2 := man.BitLen() + exp
		wp = prec + 10 - mag2
		t = lsh(man, exp+wp)
		t = rsh(mul(t, PiFixed(wp)), wp)
		n = low2(nn)
	} else {
		t, n, wp = modPi2(man, exp, mag, wp)
	}
	c, s := CosSinBasecase(t, wp)
	switch n & 3 {
	case 1:
		c, s = neg(s), c
	case 2:
		c, s = neg(c), neg(s)
	case 3:
		c, s = s, neg(c)
	}
	if sign == 1 {
		s = neg(s)
	}
	switch which {
	case 0:
		return FromManExp(c, -wp, prec, rnd), FromManExp(s, -wp, prec, rnd)
	case 1:
		return FromManExp(c, -wp, prec, rnd), MPF{}
	case 2:
		return FromManExp(s, -wp, prec, rnd), MPF{}
	}
	return FromRational(s, c, prec, rnd), MPF{}
}

// MpfCosSin returns cos(x) and sin(x).
func MpfCosSin(x MPF, prec int, rnd Rnd) (MPF, MPF) { return mpfCosSin(x, prec, rnd, 0, false) }

// MpfCos is cos(x).
func MpfCos(x MPF, prec int, rnd Rnd) MPF { v, _ := mpfCosSin(x, prec, rnd, 1, false); return v }

// MpfSin is sin(x).
func MpfSin(x MPF, prec int, rnd Rnd) MPF { v, _ := mpfCosSin(x, prec, rnd, 2, false); return v }

// MpfTan is tan(x).
func MpfTan(x MPF, prec int, rnd Rnd) MPF { v, _ := mpfCosSin(x, prec, rnd, 3, false); return v }

// MpfCosSinPi returns cos(pi*x) and sin(pi*x).
func MpfCosSinPi(x MPF, prec int, rnd Rnd) (MPF, MPF) { return mpfCosSin(x, prec, rnd, 0, true) }

// MpfCosPi is cos(pi*x).
func MpfCosPi(x MPF, prec int, rnd Rnd) MPF { v, _ := mpfCosSin(x, prec, rnd, 1, true); return v }

// MpfSinPi is sin(pi*x).
func MpfSinPi(x MPF, prec int, rnd Rnd) MPF { v, _ := mpfCosSin(x, prec, rnd, 2, true); return v }

// CosSinFixed is cos and sin of a fixed-point x.
func CosSinFixed(x *big.Int, prec int) (*big.Int, *big.Int) {
	pi2 := PiFixed(prec - 1)
	n, t := divmod(x, pi2)
	c, s := CosSinBasecase(t, prec)
	switch low2(n) {
	case 1:
		return neg(s), c
	case 2:
		return neg(c), neg(s)
	case 3:
		return s, neg(c)
	}
	return c, s
}

// ExpFixed is exp of a fixed-point x.
func ExpFixed(x *big.Int, prec int) *big.Int {
	ln2 := Ln2Fixed(prec)
	n, t := divmod(x, ln2)
	v := ExpBasecase(t, prec)
	return lsh(v, int(n.Int64()))
}
