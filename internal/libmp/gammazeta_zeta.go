// Copyright (c) 2005-2026 Fredrik Johansson and mpmath contributors.
// Go port distributed under the same BSD 3-clause licence (see LICENSE).

package libmp

import (
	"math/big"
)

// The Riemann zeta function from gammazeta.py.

var borweinCache = map[int][]*big.Int{}

func borweinCoefficients(n int) []*big.Int {
	if ds, ok := borweinCache[n]; ok {
		return ds
	}
	ds := make([]*big.Int, n+1)
	d := I(1)
	s := I(1)
	ds[0] = s
	for i := 1; i <= n; i++ {
		d = mul(d, I(4*int64(n+i-1)*int64(n-i+1)))
		d = floorDiv(d, I(int64(2*i)*int64(2*i-1)))
		s = add(s, d)
		ds[i] = s
	}
	borweinCache[n] = ds
	return ds
}

type zetaIntEntry struct {
	wp  int
	val MPF
}

var zetaIntCache = map[int]zetaIntEntry{}

// MpfZetaInt is zeta(s) for an integer s.
func MpfZetaInt(s int, prec int, rnd Rnd) MPF {
	wp := prec + 20
	if c, ok := zetaIntCache[s]; ok && c.wp >= wp {
		return MpfPos(c.val, prec, rnd)
	}
	if s < 2 {
		if s == 1 {
			Raise("ValueError", "zeta(1) pole")
		}
		if s == 0 {
			return MpfNeg(Fhalf, 0, RoundFast)
		}
		return MpfDiv(MpfBernoulli(-s+1, wp, RoundFast, false), FromInt64(int64(s-1)), prec, rnd)
	}
	if s >= wp {
		return MpfPerturb(Fone, 0, prec, rnd)
	} else if float64(s) >= float64(wp)*0.431 {
		one := lsh(bigOne, wp)
		t := add(one, lsh(bigOne, wp-s))
		t = add(t, floorDiv(one, ipow(3, s)))
		t = add(t, lsh(bigOne, maxInt(0, wp-s*2)))
		return FromManExp(t, -wp, prec, rnd)
	} else {
		m := float64(wp)/float64(s-1) + 1
		if m < 30 {
			neededTerms := int(PyPow(2.0, m) + 1)
			if float64(neededTerms) < float64(int(float64(wp)/2.54+5))/10 {
				t := Fone
				for _, k := range ListPrimes(neededTerms) {
					powprec := int(float64(wp) - float64(float64(s)*pyLog2Int(k)))
					if powprec < 2 {
						break
					}
					a := MpfSub(Fone, MpfPowInt(FromInt64(int64(k)), I(int64(-s)), powprec, RoundFast), wp, RoundFast)
					t = MpfMul(t, a, wp, RoundFast)
				}
				return MpfDiv(Fone, t, wp, RoundFast)
			}
		}
	}
	n := int(float64(wp)/2.54 + 5)
	d := borweinCoefficients(n)
	t := I(0)
	for k := 0; k < n; k++ {
		v := sub(d[k], d[n])
		if k&1 == 1 {
			v = neg(v)
		}
		t = add(t, floorDiv(lsh(v, wp), ipow(int64(k+1), s)))
	}
	t = floorDiv(lsh(t, wp), neg(d[n]))
	t = floorDiv(lsh(t, wp), sub(lsh(bigOne, wp), lsh(bigOne, wp+1-s)))
	if c, ok := zetaIntCache[s]; !ok || c.wp < wp {
		zetaIntCache[s] = zetaIntEntry{wp, FromManExp(t, -wp-wp, 0, RoundFast)}
	}
	return FromManExp(t, -wp-wp, prec, rnd)
}

// MpfZeta is zeta(s) (or the Dirichlet eta function when alt is set).
func MpfZeta(s MPF, prec int, rnd Rnd, alt bool) MPF {
	sign, man, exp, bc := s.Sign, s.Man, s.Exp, s.BC
	f := RoundFast
	if man.Sign() == 0 {
		if s.Eq(Fzero) {
			if alt {
				return Fhalf
			}
			return MpfNeg(Fhalf, 0, f)
		}
		if s.Eq(Finf) {
			return Fone
		}
		return Fnan
	}
	wp := prec + 20
	altInt := 0
	if alt {
		altInt = 1
	}
	if sign == 0 && float64(exp+bc) > PyLog(float64(wp))/PyLog(2)+2 {
		return MpfPerturb(Fone, altInt, prec, rnd)
	} else if exp >= 0 {
		if alt {
			if s.Eq(Fone) {
				return MpfLn2(prec, rnd)
			}
			z := MpfZetaInt(int(ToInt(s, f).Int64()), wp, NegativeRnd(rnd))
			q := MpfSub(Fone, MpfPow(Ftwo, MpfSub(Fone, s, wp, f), wp, f), wp, f)
			return MpfMul(z, q, prec, rnd)
		}
		return MpfZetaInt(int(ToInt(s, f).Int64()), prec, rnd)
	}
	if sign == 1 {
		if alt {
			q := MpfSub(Fone, MpfPow(Ftwo, MpfSub(Fone, s, wp, f), wp, f), wp, f)
			return MpfMul(MpfZeta(s, wp, f, false), q, prec, rnd)
		}
		y := MpfSub(Fone, s, 10*wp, f)
		a := MpfGamma(y, wp, f)
		b := MpfZeta(y, wp, f, false)
		c := MpfSinPi(MpfShift(s, -1), wp, f)
		wp2 := wp + maxInt(0, exp+bc)
		pi := MpfPi(wp+wp2, f)
		d := MpfDiv(MpfPow(MpfShift(pi, 1), s, wp2, f), pi, wp2, f)
		return MpfMul(a, MpfMul(b, MpfMul(c, d, wp, f), wp, f), prec, rnd)
	}
	r := MpfSub(Fone, s, wp, f)
	ar := MpfAbs(r, 0, f)
	poleDist := -2 * (ar.Exp + ar.BC)
	if poleDist > wp {
		if alt {
			return MpfLn2(prec, rnd)
		}
		q := MpfNeg(MpfDiv(Fone, r, wp, f), 0, f)
		return MpfAdd(q, MpfEuler(wp, f), prec, rnd)
	}
	wp += maxInt(0, poleDist)
	n := int(float64(wp)/2.54 + 5)
	d := borweinCoefficients(n)
	t := I(0)
	sf := ToFixed(s, wp)
	for k := 0; k < n; k++ {
		u := rsh(mul(neg(sf), LogIntFixed(k+1, wp)), wp)
		eman := ExpFixed(u, wp)
		w := mul(sub(d[k], d[n]), eman)
		if k&1 == 1 {
			t = sub(t, w)
		} else {
			t = add(t, w)
		}
	}
	t = floorDiv(t, neg(d[n]))
	tf := FromManExp(t, -wp, wp, f)
	if alt {
		return MpfPos(tf, prec, rnd)
	}
	q := MpfSub(Fone, MpfPow(Ftwo, MpfSub(Fone, s, wp, f), wp, f), wp, f)
	return MpfDiv(tf, q, prec, rnd)
}

// MpcZeta is zeta(s) for complex s; force skips the size guard.
func MpcZeta(s MPC, prec int, rnd Rnd, alt bool, force bool) MPC {
	re, im := s.Re, s.Im
	f := RoundFast
	if im.Eq(Fzero) {
		return MPC{MpfZeta(re, prec, rnd, alt), Fzero}
	}
	if !force && MpfGt(MpcAbs(s, 10, f), FromInt64(int64(prec))) {
		Raise("NotImplementedError", "")
	}
	wp := prec + 20
	r := MpcSub(MpcOne, s, wp, f)
	ar := MpcAbs(r, 10, f)
	poleDist := -2 * (ar.Exp + ar.BC)
	if poleDist > wp {
		if alt {
			q := MpfLn2(wp, f)
			y := MpfMul(q, MpfEuler(wp, f), wp, f)
			g := MpfShift(MpfMul(q, q, wp, f), -1)
			g = MpfSub(y, g, 0, f)
			z := MpcMulMpf(r, MpfNeg(g, 0, f), wp, f)
			z = MpcAddMpf(z, q, wp, f)
			return MpcPos(z, prec, rnd)
		}
		q := MpcNeg(MpcDiv(MpcOne, r, wp, f), 0, f)
		q = MpcAddMpf(q, MpfEuler(wp, f), wp, f)
		return MpcPos(q, prec, rnd)
	}
	wp += maxInt(0, poleDist)
	if MpfLt(re, Fzero) {
		if alt {
			q := MpcSub(MpcOne, MpcPow(MpcTwo, MpcSub(MpcOne, s, wp, f), wp, f), wp, f)
			return MpcMul(MpcZeta(s, wp, f, false, false), q, prec, rnd)
		}
		y := MpcSub(MpcOne, s, 10*wp, f)
		a := MpcGamma(y, wp, f)
		b := MpcZeta(y, wp, f, false, false)
		c := MpcSinPi(MpcShift(s, -1), wp, f)
		mag := maxInt(re.Exp+re.BC, im.Exp+im.BC)
		wp2 := wp + maxInt(0, mag)
		pi := MpfPi(wp+wp2, f)
		pi2 := MPC{MpfShift(pi, 1), Fzero}
		d := MpcDivMpf(MpcPow(pi2, s, wp2, f), pi, wp2, f)
		return MpcMul(a, MpcMul(b, MpcMul(c, d, wp, f), wp, f), prec, rnd)
	}
	n := int(float64(wp)/2.54 + 5)
	imAbs, _ := new(big.Float).SetInt(absI(ToInt(im, f))).Float64()
	n += int(0.9 * imAbs)
	d := borweinCoefficients(n)
	ref := ToFixed(re, wp)
	imf := ToFixed(im, wp)
	tre, tim := I(0), I(0)
	one2wp := lsh(bigOne, 2*wp)
	criticalLine := re.Eq(Fhalf)
	wp2 := wp + wp
	for k := 0; k < n; k++ {
		log := LogIntFixed(k+1, wp)
		var w *big.Int
		if criticalLine {
			w = floorDiv(one2wp, Isqrt(lsh(I(int64(k+1)), wp2)))
		} else {
			w = ExpFixed(rsh(mul(neg(ref), log), wp), wp)
		}
		if k&1 == 1 {
			w = mul(w, sub(d[n], d[k]))
		} else {
			w = mul(w, sub(d[k], d[n]))
		}
		wre, wim := CosSinFixed(rsh(mul(neg(imf), log), wp), wp)
		tre = add(tre, rsh(mul(w, wre), wp))
		tim = add(tim, rsh(mul(w, wim), wp))
	}
	tre = floorDiv(tre, neg(d[n]))
	tim = floorDiv(tim, neg(d[n]))
	t := MPC{FromManExp(tre, -wp, wp, f), FromManExp(tim, -wp, wp, f)}
	if alt {
		return MpcPos(t, prec, rnd)
	}
	q := MpcSub(MpcOne, MpcPow(MpcTwo, r, wp, f), wp, f)
	return MpcDiv(t, q, prec, rnd)
}

func MpfAltzeta(s MPF, prec int, rnd Rnd) MPF { return MpfZeta(s, prec, rnd, true) }
func MpcAltzeta(s MPC, prec int, rnd Rnd) MPC { return MpcZeta(s, prec, rnd, true, false) }

func powFixed(x *big.Int, n int, wp int) *big.Int {
	if n == 1 {
		return x
	}
	y := lsh(bigOne, wp)
	for n != 0 {
		if n&1 == 1 {
			y = rsh(mul(y, x), wp)
			n--
		}
		x = rsh(mul(x, x), wp)
		n /= 2
	}
	return y
}

var sieveCache, primesCache, multCache []int

func primesieve(n int) ([]int, []int, []int) {
	if n < len(sieveCache) {
		return sieveCache, primesCache, multCache
	}
	sieve := make([]int, n+1)
	mult := make([]int, n+1)
	primes := ListPrimes(n)
	for _, p := range primes {
		for k := p; k <= n; k += p {
			sieve[k] = p
		}
	}
	for i, p := range sieve {
		if i >= 2 {
			m := 1
			q := i / p
			for q%p == 0 {
				q /= p
				m++
			}
			mult[i] = m
		}
	}
	sieveCache, primesCache, multCache = sieve, primes, mult
	return sieve, primes, mult
}

func zetasumSieved(criticalLine bool, sre, sim *big.Int, a, n, wp int) (*big.Int, *big.Int) {
	if a < 1 {
		Raise("ValueError", "a cannot be less than 1")
	}
	sieve, primes, mult := primesieve(a + n)
	basicPowers := map[int][][2]*big.Int{}
	one := lsh(bigOne, wp)
	one2wp := lsh(bigOne, 2*wp)
	wp2 := wp + wp
	for _, p := range primes {
		if p*2 > a+n {
			break
		}
		log := LogIntFixed(p, wp)
		cos, sin := CosSinFixed(rsh(mul(neg(sim), log), wp), wp)
		var u *big.Int
		if criticalLine {
			u = floorDiv(one2wp, Isqrt(lsh(I(int64(p)), wp2)))
		} else {
			u = ExpFixed(rsh(mul(neg(sre), log), wp), wp)
		}
		pre := rsh(mul(u, cos), wp)
		pim := rsh(mul(u, sin), wp)
		basicPowers[p] = [][2]*big.Int{{pre, pim}}
		tre, tim := pre, pim
		top := int(PyLog(float64(a+n))/PyLog(float64(p)) + 0.01)
		for m := 1; m <= top; m++ {
			tre, tim = rsh(sub(mul(pre, tre), mul(pim, tim)), wp), rsh(add(mul(pim, tre), mul(pre, tim)), wp)
			basicPowers[p] = append(basicPowers[p], [2]*big.Int{tre, tim})
		}
	}
	xre, xim := I(0), I(0)
	if a == 1 {
		xre = add(xre, one)
	}
	aa := maxInt(a, 2)
	for k := aa; k <= a+n; k++ {
		p := sieve[k]
		var tre, tim *big.Int
		if bp, ok := basicPowers[p]; ok {
			m := mult[k]
			tre, tim = bp[m-1][0], bp[m-1][1]
			kk := k
			for {
				pm := 1
				for i := 0; i < m; i++ {
					pm *= p
				}
				kk /= pm
				if kk == 1 {
					break
				}
				p = sieve[kk]
				m = mult[kk]
				pre, pim := basicPowers[p][m-1][0], basicPowers[p][m-1][1]
				tre, tim = rsh(sub(mul(pre, tre), mul(pim, tim)), wp), rsh(add(mul(pim, tre), mul(pre, tim)), wp)
			}
		} else {
			log := LogIntFixed(k, wp)
			cos, sin := CosSinFixed(rsh(mul(neg(sim), log), wp), wp)
			var u *big.Int
			if criticalLine {
				u = floorDiv(one2wp, Isqrt(lsh(I(int64(k)), wp2)))
			} else {
				u = ExpFixed(rsh(mul(neg(sre), log), wp), wp)
			}
			tre = rsh(mul(u, cos), wp)
			tim = rsh(mul(u, sin), wp)
		}
		xre = add(xre, tre)
		xim = add(xim, tim)
	}
	return xre, xim
}

const zetasumSieveCutoff = 10

// MpcZetasum is mpc_zetasum: the sums of (a+k)^-s (times log powers for the
// requested derivatives) for k = 0..n, and the reflected sums if asked.
func MpcZetasum(s MPC, a, n int, derivatives []int, reflect bool, prec int) ([]MPC, []MPC) {
	wp := prec + 10
	haveDerivatives := !(len(derivatives) == 1 && derivatives[0] == 0)
	haveOneDerivative := len(derivatives) == 1
	criticalLine := s.Re.Eq(Fhalf)
	sre := ToFixed(s.Re, wp)
	sim := ToFixed(s.Im, wp)
	if a > 0 && n > zetasumSieveCutoff && !haveDerivatives && !reflect {
		re, im := zetasumSieved(criticalLine, sre, sim, a, n, wp)
		return []MPC{{FromManExp(re, -wp, prec, RoundNearest), FromManExp(im, -wp, prec, RoundNearest)}}, nil
	}
	maxd := derivatives[0]
	for _, d := range derivatives {
		maxd = maxInt(maxd, d)
	}
	if !haveOneDerivative {
		derivatives = make([]int, maxd+1)
		for i := range derivatives {
			derivatives[i] = i
		}
	}
	nd := len(derivatives)
	zeros := func() []*big.Int {
		v := make([]*big.Int, nd)
		for i := range v {
			v[i] = I(0)
		}
		return v
	}
	xre, xim := zeros(), zeros()
	var yre, yim []*big.Int
	if reflect {
		yre, yim = zeros(), zeros()
	}
	one2wp := lsh(bigOne, 2*wp)
	wp2 := wp + wp
	for w := a; w <= a+n; w++ {
		log := LogIntFixed(w, wp)
		cos, sin := CosSinFixed(rsh(mul(neg(sim), log), wp), wp)
		var u *big.Int
		if criticalLine {
			u = floorDiv(one2wp, Isqrt(lsh(I(int64(w)), wp2)))
		} else {
			u = ExpFixed(rsh(mul(neg(sre), log), wp), wp)
		}
		xtre := rsh(mul(u, cos), wp)
		xtim := rsh(mul(u, sin), wp)
		var ytre, ytim *big.Int
		if reflect {
			reciprocal := floorDiv(one2wp, muli(u, int64(w)))
			ytre = rsh(mul(reciprocal, cos), wp)
			ytim = rsh(mul(reciprocal, sin), wp)
		}
		if haveDerivatives {
			if haveOneDerivative {
				log = powFixed(log, maxd, wp)
				xre[0] = add(xre[0], rsh(mul(xtre, log), wp))
				xim[0] = add(xim[0], rsh(mul(xtim, log), wp))
				if reflect {
					yre[0] = add(yre[0], rsh(mul(ytre, log), wp))
					yim[0] = add(yim[0], rsh(mul(ytim, log), wp))
				}
			} else {
				t := lsh(bigOne, wp)
				for _, d := range derivatives {
					xre[d] = add(xre[d], rsh(mul(xtre, t), wp))
					xim[d] = add(xim[d], rsh(mul(xtim, t), wp))
					if reflect {
						yre[d] = add(yre[d], rsh(mul(ytre, t), wp))
						yim[d] = add(yim[d], rsh(mul(ytim, t), wp))
					}
					t = rsh(mul(t, log), wp)
				}
			}
		} else {
			xre[0] = add(xre[0], xtre)
			xim[0] = add(xim[0], xtim)
			if reflect {
				yre[0] = add(yre[0], ytre)
				yim[0] = add(yim[0], ytim)
			}
		}
	}
	if haveDerivatives {
		if haveOneDerivative {
			if maxd%2 == 1 {
				xre[0], xim[0] = neg(xre[0]), neg(xim[0])
				if reflect {
					yre[0], yim[0] = neg(yre[0]), neg(yim[0])
				}
			}
		} else {
			for _, d := range derivatives {
				if d&1 == 1 {
					xre[d], xim[d] = neg(xre[d]), neg(xim[d])
					if reflect {
						yre[d], yim[d] = neg(yre[d]), neg(yim[d])
					}
				}
			}
		}
	}
	var xs, ys []MPC
	for i := range xre {
		xs = append(xs, MPC{FromManExp(xre[i], -wp, prec, RoundNearest), FromManExp(xim[i], -wp, prec, RoundNearest)})
	}
	for i := range yre {
		ys = append(ys, MPC{FromManExp(yre[i], -wp, prec, RoundNearest), FromManExp(yim[i], -wp, prec, RoundNearest)})
	}
	return xs, ys
}
