// Copyright (c) 2005-2026 Fredrik Johansson and mpmath contributors.
// Go port distributed under the same BSD 3-clause licence (see LICENSE).

package libmp

import "math/big"

// Taylor and Stirling series behind the gamma function (gammazeta.py).

const (
	maxGammaTaylorPrec      = 5000
	gammaStirlingBeta       = 0.2
	smallFactorialCacheSize = 150
)

// gammaTaylorCache keeps the original dict's insertion order, which decides
// which higher-precision entry gets truncated for reuse.
var (
	gammaTaylorCache = map[int][]*big.Int{}
	gammaTaylorOrder []int
)

type stirlingCoef struct {
	p, q   *big.Int
	pb, qb int
}

var gammaStirlingCache = map[int]stirlingCoef{}

var smallFactorialCache = func() []MPF {
	v := make([]MPF, smallFactorialCacheSize+1)
	for n := range v {
		v[n] = FromInt(Ifac(n), 0, RoundFast)
	}
	return v
}()

func zetaArray(N int, prec int) []*big.Int {
	extra := 30
	wp := prec + extra
	f := RoundFast
	zv := make([]*big.Int, N+2)
	for i := range zv {
		zv[i] = I(0)
	}
	pi := PiFixed(wp)
	one := lsh(bigOne, wp)
	zv[0] = floorDiv(neg(one), bigTwo)
	f2pi := MpfShift(MpfPi(wp, f), 1)
	exp2pi := MpfExp(f2pi, wp, f)
	exp2piK := exp2pi
	var exps3 [][2]*big.Int
	for k := 1; ; k++ {
		tp := wp - 9*k
		if tp < 1 {
			break
		}
		q1 := MpfDiv(Fone, MpfSub(exp2piK, Fone, tp, f), tp, f)
		q2 := MpfMul(exp2piK, MpfMul(q1, q1, tp, f), tp, f)
		q1i := ToFixed(q1, wp)
		q2i := ToFixed(q2, wp)
		q2i = rsh(mul(muli(q2i, int64(k)), pi), wp)
		exps3 = append(exps3, [2]*big.Int{q1i, q2i})
		exp2piK = MpfMul(exp2piK, exp2pi, wp, f)
	}
	for n := 3; n <= N; n += 2 {
		s := I(0)
		k := 1
		for _, e := range exps3 {
			kn := ipow(int64(k), n)
			var t *big.Int
			if n%4 == 3 {
				t = floorDiv(e[0], kn)
			} else {
				U := I(int64((n - 1) / 4))
				t = floorDiv(add(e[0], floorDiv(e[1], U)), kn)
			}
			if t.Sign() == 0 {
				break
			}
			s = add(s, t)
			k++
		}
		zv[n] = muli(s, -2)
	}
	B := make([]MPF, N+2)
	for k := range B {
		B[k] = MpfAbs(MpfBernoulli(k, wp, f, false), 0, f)
	}
	fpi := MpfPowInt(MpfShift(MpfPi(wp, f), 1), I(2), wp, f)
	piPow := MpfDiv(fpi, FromInt64(4), wp, f)
	for n := 2; n < N+2; n += 2 {
		z := MpfMul(B[n], piPow, wp, f)
		zv[n] = ToFixed(z, wp)
		piPow = MpfMul(piPow, fpi, wp, f)
		piPow = MpfDiv(piPow, FromInt64(int64((n+1)*(n+2))), wp, f)
	}
	reciprocalPi := floorDiv(lsh(one, wp), pi)
	for n := 3; n <= N; n += 4 {
		U := (n - 3) / 4
		s := floorDiv(muli(zv[4*U+4], int64(4*U+7)), I(4))
		for k := 1; k <= U; k++ {
			s = sub(s, rsh(mul(zv[4*k], zv[4*U+4-4*k]), wp))
		}
		zv[n] = add(zv[n], rsh(muli(mul(s, reciprocalPi), 2), wp))
	}
	for n := 5; n <= N; n += 4 {
		U := (n - 1) / 4
		s := muli(zv[4*U+2], int64(2*U+1))
		for k := 1; k <= 2*U; k++ {
			t := muli(mul(zv[2*k], zv[4*U+2-2*k]), int64(2*k))
			if k&1 == 1 {
				t = neg(t)
			}
			s = add(s, rsh(t, wp))
		}
		zv[n] = add(zv[n], floorDiv(rsh(mul(s, reciprocalPi), wp), I(int64(2*U))))
	}
	for i := range zv {
		zv[i] = rsh(zv[i], extra)
	}
	return zv
}

func gammaTaylorCoefficients(inprec int) ([]*big.Int, int) {
	var prec int
	switch {
	case inprec < 400:
		prec = inprec + (10 - pyModInt(inprec, 10))
	case inprec < 1000:
		prec = inprec + (30 - pyModInt(inprec, 30))
	default:
		prec = inprec
	}
	if c, ok := gammaTaylorCache[prec]; ok {
		return c, prec
	}
	var N int
	if prec < 1000 {
		N = int(PyPow(float64(prec), 0.76) + 2)
	} else {
		N = int(PyPow(float64(prec), 0.787) + 2)
	}
	for _, cprec := range gammaTaylorOrder {
		if cprec > prec {
			src := gammaTaylorCache[cprec]
			start := len(src) - N
			if start < 0 {
				start = 0
			}
			coeffs := make([]*big.Int, 0, len(src)-start)
			for _, x := range src[start:] {
				coeffs = append(coeffs, rsh(x, cprec-prec))
			}
			if inprec < 1000 {
				gammaTaylorStore(prec, coeffs)
			}
			return coeffs, prec
		}
	}
	if prec > 1000 {
		prec = int(float64(prec) * 1.2)
	}
	wp := prec + 20
	A := make([]*big.Int, N)
	for i := range A {
		A[i] = I(0)
	}
	A[1] = lsh(bigOne, wp)
	A[2] = EulerFixed(wp)
	zv := zetaArray(N, wp)
	for k := 3; k < N; k++ {
		a := rsh(mul(neg(A[2]), A[k-1]), wp)
		for j := 2; j < k; j++ {
			t := mul(zv[j], A[k-j])
			if j&1 == 1 {
				t = neg(t)
			}
			a = add(a, rsh(t, wp))
		}
		a = floorDiv(a, I(int64(1-k)))
		A[k] = a
	}
	out := make([]*big.Int, 0, N)
	for i := len(A) - 1; i >= 1; i-- {
		out = append(out, rsh(A[i], 20))
	}
	gammaTaylorStore(prec, out)
	return gammaTaylorCoefficients(inprec)
}

func gammaTaylorStore(prec int, v []*big.Int) {
	if _, ok := gammaTaylorCache[prec]; !ok {
		gammaTaylorOrder = append(gammaTaylorOrder, prec)
	}
	gammaTaylorCache[prec] = v
}

func gammaFixedTaylor(xmpf MPF, x *big.Int, wp, prec int, rnd Rnd, typ int) MPF {
	nearestInt := rsh(addi(rsh(x, wp-1), 1), 1)
	one := lsh(bigOne, wp)
	coeffs, cwp := gammaTaylorCoefficients(wp)
	f := RoundFast
	if nearestInt.Sign() > 0 {
		r := one
		for i := int64(0); i < nearestInt.Int64()-1; i++ {
			x = sub(x, one)
			r = rsh(mul(r, x), wp)
		}
		x = sub(x, one)
		p := I(0)
		for _, c := range coeffs {
			p = add(c, rsh(mul(x, p), wp))
		}
		p = rsh(p, cwp-wp)
		switch typ {
		case 0:
			return FromManExp(floorDiv(lsh(r, wp), p), -wp, prec, rnd)
		case 2:
			return MpfShift(FromRational(p, lsh(r, wp), prec, rnd), wp)
		case 3:
			return MpfLn(MpfAbs(FromManExp(floorDiv(lsh(r, wp), p), -wp, 0, f), 0, f), prec, rnd)
		}
		return MPF{}
	}
	r := one
	for i := int64(0); i < -nearestInt.Int64(); i++ {
		r = rsh(mul(r, x), wp)
		x = add(x, one)
	}
	p := I(0)
	for _, c := range coeffs {
		p = add(c, rsh(mul(x, p), wp))
	}
	p = rsh(p, cwp-wp)
	if wp-x.BitLen() > 10 {
		g := MpfAdd(xmpf, FromInt(neg(nearestInt), 0, f), 0, f)
		rr := FromManExp(mul(p, r), -wp-wp, 0, f)
		rr = MpfMul(rr, g, wp, f)
		switch typ {
		case 0:
			return MpfDiv(Fone, rr, prec, rnd)
		case 2:
			return MpfPos(rr, prec, rnd)
		case 3:
			return MpfLn(MpfAbs(MpfDiv(Fone, rr, wp, f), 0, f), prec, rnd)
		}
		return MPF{}
	}
	rr := FromManExp(mul(mul(x, p), r), -3*wp, 0, f)
	switch typ {
	case 0:
		return MpfDiv(Fone, rr, prec, rnd)
	case 2:
		return MpfPos(rr, prec, rnd)
	case 3:
		return MpfNeg(MpfLn(MpfAbs(rr, 0, f), prec, rnd), 0, f)
	}
	return MPF{}
}

func stirlingCoefficient(n int) stirlingCoef {
	if c, ok := gammaStirlingCache[n]; ok {
		return c
	}
	p, q := Bernfrac(n, false)
	q = mul(q, I(int64(n)*int64(n-1)))
	c := stirlingCoef{p, q, p.BitLen(), q.BitLen()}
	gammaStirlingCache[n] = c
	return c
}

func stirlingTerm(t *big.Int, mulBy, div int64) *big.Int {
	return floorDiv(muli(t, mulBy), I(div))
}

func realStirlingSeries(x *big.Int, prec int) *big.Int {
	t := floorDiv(lsh(bigOne, prec+prec), x)
	u := rsh(mul(t, t), prec)
	s := sub(LnSqrt2piFixed(prec), x)
	step := func() { t = rsh(mul(t, u), prec) }
	s = add(s, stirlingTerm(t, 1, 12))
	step()
	s = sub(s, stirlingTerm(t, 1, 360))
	step()
	s = add(s, stirlingTerm(t, 1, 1260))
	step()
	s = sub(s, stirlingTerm(t, 1, 1680))
	step()
	if t.Sign() == 0 {
		return s
	}
	s = add(s, stirlingTerm(t, 1, 1188))
	step()
	s = sub(s, stirlingTerm(t, 691, 360360))
	step()
	s = add(s, stirlingTerm(t, 1, 156))
	step()
	if t.Sign() == 0 {
		return s
	}
	s = sub(s, stirlingTerm(t, 3617, 122400))
	step()
	s = add(s, stirlingTerm(t, 43867, 244188))
	step()
	s = sub(s, stirlingTerm(t, 174611, 125400))
	step()
	if t.Sign() == 0 {
		return s
	}
	k := 22
	usize := u.BitLen()
	texp := 0
	for {
		c := stirlingCoefficient(k)
		shift := -texp
		term := rsh(floorDiv(mul(t, c.p), c.q), shift)
		if term.Sign() == 0 {
			break
		}
		s = add(s, term)
		t = rsh(mul(t, u), usize)
		texp -= prec - usize
		k += 2
	}
	return s
}

func complexStirlingSeries(x, y *big.Int, prec int) (*big.Int, *big.Int) {
	m := rsh(add(mul(x, x), mul(y, y)), prec)
	tre := floorDiv(lsh(x, prec), m)
	tim := floorDiv(lsh(neg(y), prec), m)
	ure := rsh(sub(mul(tre, tre), mul(tim, tim)), prec)
	uim := rsh(mul(tim, tre), prec-1)
	sre := sub(LnSqrt2piFixed(prec), x)
	sim := neg(y)
	step := func() {
		tre, tim = rsh(sub(mul(tre, ure), mul(tim, uim)), prec), rsh(add(mul(tre, uim), mul(tim, ure)), prec)
	}
	addT := func(sign int, mulBy, div int64) {
		a, b := stirlingTerm(tre, mulBy, div), stirlingTerm(tim, mulBy, div)
		if sign > 0 {
			sre, sim = add(sre, a), add(sim, b)
		} else {
			sre, sim = sub(sre, a), sub(sim, b)
		}
	}
	small := func() bool { return add(absI(tre), absI(tim)).Cmp(I(5)) < 0 }
	addT(1, 1, 12)
	step()
	addT(-1, 1, 360)
	step()
	addT(1, 1, 1260)
	step()
	addT(-1, 1, 1680)
	step()
	if small() {
		return sre, sim
	}
	addT(1, 1, 1188)
	step()
	addT(-1, 691, 360360)
	step()
	addT(1, 1, 156)
	step()
	if small() {
		return sre, sim
	}
	addT(-1, 3617, 122400)
	step()
	addT(1, 43867, 244188)
	step()
	addT(-1, 174611, 125400)
	step()
	if small() {
		return sre, sim
	}
	k := 22
	mx := absI(ure)
	if absI(uim).Cmp(mx) > 0 {
		mx = absI(uim)
	}
	usize := mx.BitLen()
	texp := 0
	for {
		c := stirlingCoefficient(k)
		shift := -texp
		termre := rsh(floorDiv(mul(tre, c.p), c.q), shift)
		termim := rsh(floorDiv(mul(tim, c.p), c.q), shift)
		if add(absI(termre), absI(termim)).Cmp(I(5)) < 0 {
			break
		}
		sre = add(sre, termre)
		sim = add(sim, termim)
		tre, tim = rsh(sub(mul(tre, ure), mul(tim, uim)), usize), rsh(add(mul(tre, uim), mul(tim, ure)), usize)
		texp -= prec - usize
		k += 2
	}
	return sre, sim
}
