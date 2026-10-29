// Copyright (c) 2005-2026 Fredrik Johansson and mpmath contributors.
// Go port distributed under the same BSD 3-clause licence (see LICENSE).

package libmp

import (
	"math"
	"math/big"
)

// Constants and Bernoulli numbers from mpmath/libmp/gammazeta.py.

// pyLog2Int is math.log(n, 2) for a positive int n.
func pyLog2Int(n int) float64 { return PyLog(float64(n)) / PyLog(2) }

func ipow(n int64, k int) *big.Int { return new(big.Int).Exp(I(n), I(int64(k)), nil) }

var (
	catalanMemo = newMemo(func(prec int) *big.Int {
		prec += 20
		one := lsh(bigOne, prec)
		a := one
		s := I(0)
		t := I(1)
		for n := int64(1); t.Sign() != 0; n++ {
			bn := I(n)
			n3 := mul(bn, mul(bn, bn))
			a = mul(a, muli(mul(n3, I(2*n-1)), 32))
			d := add(I(3-16*n), muli(mul(bn, bn), 16))
			a = floorDiv(a, mul(d, d))
			t = mul(a, add(muli(mul(bn, bn), 40), I(-24*n+3)))
			if (n-1)&1 == 1 {
				t = neg(t)
			}
			t = floorDiv(t, mul(n3, I(2*n-1)))
			s = add(s, t)
		}
		return rsh(s, 20+6)
	})

	khinchinMemo = newMemo(func(prec int) *big.Int {
		wp := int(float64(prec) + math.Sqrt(float64(prec)) + 15)
		s := I(0)
		fac := FromInt64(4)
		ONE := lsh(bigOne, wp)
		t := ONE
		pi := MpfPi(wp, RoundFast)
		twopi2 := MpfShift(MpfMul(pi, pi, wp, RoundFast), 2)
		pipow := twopi2
		n := int64(1)
		for {
			zeta2n := MpfAbs(MpfBernoulli(int(2*n), wp, RoundFast, false), 0, RoundFast)
			zeta2n = MpfMul(zeta2n, pipow, wp, RoundFast)
			zeta2n = MpfDiv(zeta2n, fac, wp, RoundFast)
			z := ToFixed(zeta2n, wp)
			term := rsh(floorDiv(mul(sub(z, ONE), t), I(n)), wp)
			if term.Cmp(I(100)) < 0 {
				break
			}
			s = add(s, term)
			t = add(t, sub(floorDiv(ONE, I(2*n+1)), floorDiv(ONE, I(2*n))))
			n++
			fac = MpfMulInt(fac, I((2*n)*(2*n-1)), wp, RoundFast)
			pipow = MpfMul(pipow, twopi2, wp, RoundFast)
		}
		s = floorDiv(lsh(s, wp), Ln2Fixed(wp))
		K := MpfExp(FromManExp(s, -wp, 0, RoundFast), wp, RoundFast)
		return ToFixed(K, prec)
	})

	glaisherMemo = newMemo(func(prec int) *big.Int {
		wp := prec + 30
		N := int(float64(0.33*float64(prec)) + 5)
		ONE := lsh(bigOne, wp)
		s := I(0)
		for k := 2; k < N; k++ {
			s = add(s, floorDiv(LogIntFixed(k, wp), I(int64(k)*int64(k))))
		}
		logN := LogIntFixed(N, wp)
		bN := I(int64(N))
		s = add(s, floorDiv(add(ONE, logN), bN))
		s = add(s, floorDiv(logN, muli(mul(bN, bN), 2)))
		pN := mul(bN, mul(bN, bN))
		a, b := I(1), I(-2)
		j := int64(3)
		fac := FromInt64(2)
		k := int64(1)
		for {
			D := floorDiv(add(lsh(a, wp), mul(b, logN)), pN)
			Df := FromManExp(D, -wp, 0, RoundFast)
			B := MpfBernoulli(int(2*k), wp, RoundFast, false)
			term := MpfMul(B, Df, wp, RoundFast)
			term = MpfDiv(term, fac, wp, RoundFast)
			tf := ToFixed(term, wp)
			if absI(tf).Cmp(I(100)) < 0 {
				break
			}
			s = sub(s, tf)
			for r := 0; r < 2; r++ {
				a, b, pN, j = sub(b, muli(a, j)), muli(b, -j), mul(pN, bN), j+1
			}
			k++
			fac = MpfMulInt(fac, I((2*k)*(2*k-1)), wp, RoundFast)
		}
		pi := PiFixed(wp)
		s = muli(s, 6)
		s = floorDiv(lsh(s, wp), rsh(mul(pi, pi), wp))
		s = add(s, EulerFixed(wp))
		s = add(s, ToFixed(MpfLn(FromManExp(muli(pi, 2), -wp, 0, RoundFast), wp, RoundFast), wp))
		s = floorDiv(s, I(12))
		A := MpfExp(FromManExp(s, -wp, 0, RoundFast), wp, RoundFast)
		return ToFixed(A, prec)
	})

	aperyMemo = newMemo(func(prec int) *big.Int {
		prec += 20
		d := lsh(bigOne, prec)
		term := lsh(I(77), prec)
		s := I(0)
		for n := int64(1); term.Sign() != 0; n++ {
			s = add(s, term)
			d = mul(d, ipow(n, 10))
			d = floorDiv(d, mul(ipow(2*n+1, 5), ipow(2*n, 5)))
			term = mul(I(205*n*n+250*n+77), d)
			if n&1 == 1 {
				term = neg(term)
			}
		}
		return rsh(s, 20+6)
	})

	eulerMemo = newMemo(func(prec int) *big.Int {
		extra := 30
		prec += extra
		p := int(PyLog(float64(prec)/4*PyLog(2))/PyLog(2)) + 1
		n := lsh(bigOne, p)
		n2 := mul(n, n)
		A := muli(Ln2Fixed(prec), int64(-p))
		U := A
		B := lsh(bigOne, prec)
		V := B
		for k := int64(1); ; k++ {
			bk := I(k)
			B = floorDiv(mul(B, n2), mul(bk, bk))
			A = floorDiv(add(floorDiv(mul(A, n2), bk), B), bk)
			U = add(U, A)
			V = add(V, B)
			if absI(A).Cmp(I(100)) < 0 && absI(B).Cmp(I(100)) < 0 {
				break
			}
		}
		return floorDiv(lsh(U, prec-extra), V)
	})

	mertensMemo = newMemo(func(prec int) *big.Int {
		wp := prec + 20
		s := MpfEuler(wp, RoundFast)
		for m := 2; ; m++ {
			t := MpfZetaInt(m, wp, RoundFast)
			if t.Eq(Fone) {
				break
			}
			t = MpfLn(t, wp, RoundFast)
			t = MpfMulInt(t, I(int64(Moebius(m))), wp, RoundFast)
			t = MpfDiv(t, FromInt64(int64(m)), wp, RoundFast)
			s = MpfAdd(s, t, 0, RoundFast)
		}
		return ToFixed(s, prec)
	})

	twinprimeMemo = newMemo(func(prec int) *big.Int {
		In := func(n int) *big.Int {
			s := I(0)
			for d := 1; d <= n; d++ {
				if n%d == 0 {
					s = add(s, muli(lsh(bigOne, n/d), int64(Moebius(d))))
				}
			}
			return floorDiv(s, I(int64(n)))
		}
		wp := 2*prec + 30
		res := Fone
		var primes, ppowers []MPF
		for _, p := range []int64{2, 3, 5, 7} {
			primes = append(primes, FromRational(I(1), I(p), wp, RoundFast))
		}
		for _, p := range primes {
			ppowers = append(ppowers, MpfMul(p, p, wp, RoundFast))
		}
		for n := 2; ; n++ {
			a := MpfZetaInt(n, wp, RoundFast)
			for i := 0; i < 4; i++ {
				a = MpfMul(a, MpfSub(Fone, ppowers[i], 0, RoundFast), wp, RoundFast)
				ppowers[i] = MpfMul(ppowers[i], primes[i], wp, RoundFast)
			}
			a = MpfPowInt(a, neg(In(n)), wp, RoundFast)
			if MpfPos(a, prec+10, RoundNearest).Eq(Fone) {
				break
			}
			res = MpfMul(res, a, wp, RoundFast)
		}
		res = MpfMul(res, FromInt64(3*15*35), wp, RoundFast)
		res = MpfDiv(res, FromInt64(4*16*36), wp, RoundFast)
		return ToFixed(res, prec)
	})
)

func CatalanFixed(prec int) *big.Int   { return catalanMemo.get(prec) }
func KhinchinFixed(prec int) *big.Int  { return khinchinMemo.get(prec) }
func GlaisherFixed(prec int) *big.Int  { return glaisherMemo.get(prec) }
func AperyFixed(prec int) *big.Int     { return aperyMemo.get(prec) }
func EulerFixed(prec int) *big.Int     { return eulerMemo.get(prec) }
func MertensFixed(prec int) *big.Int   { return mertensMemo.get(prec) }
func TwinprimeFixed(prec int) *big.Int { return twinprimeMemo.get(prec) }

var (
	MpfEuler     = defMpfConstant(EulerFixed)
	MpfApery     = defMpfConstant(AperyFixed)
	MpfKhinchin  = defMpfConstant(KhinchinFixed)
	MpfGlaisher  = defMpfConstant(GlaisherFixed)
	MpfCatalan   = defMpfConstant(CatalanFixed)
	MpfMertens   = defMpfConstant(MertensFixed)
	MpfTwinprime = defMpfConstant(TwinprimeFixed)
)

const maxBernoulliCache = 3000

type bernState struct {
	numbers   map[int]MPF
	m         int
	bin, bin1 *big.Int
}

var bernoulliCache = map[int]*bernState{}

// bernoulliSize estimates the size of B_n (even n > 2 only).
func bernoulliSize(n int) int {
	lgn := pyLog2Int(n)
	return int(2.326 + float64(0.5*lgn) + float64(float64(n)*(lgn-4.094)))
}

var bernoulliPrecCutoff = bernoulliSize(maxBernoulliCache)

// MpfBernoulli is the n-th Bernoulli number.
func MpfBernoulli(n int, prec int, rnd Rnd, plus bool) MPF {
	if n < 2 {
		if n < 0 {
			Raise("ValueError", "Bernoulli numbers only defined for n >= 0")
		}
		if n == 0 {
			return Fone
		}
		if plus {
			return Fhalf
		}
		return MpfNeg(Fhalf, 0, RoundFast)
	}
	if n&1 == 1 {
		return Fzero
	}
	if prec > bernoulliPrecCutoff && float64(prec) > float64(float64(bernoulliSize(n))*1.1)+1000 {
		p, q := Bernfrac(n, false)
		return FromRational(p, q, prec, rnd)
	}
	if n > maxBernoulliCache {
		return mpfBernoulliHuge(n, prec, rnd)
	}
	wp := prec + 30
	wp += 32 - (prec & 31)
	st, ok := bernoulliCache[wp]
	if ok {
		if v, ok := st.numbers[n]; ok {
			return MpfPos(v, prec, rnd)
		}
		if n-st.m > 10 {
			return mpfBernoulliHuge(n, prec, rnd)
		}
	} else {
		if n > 10 {
			return mpfBernoulliHuge(n, prec, rnd)
		}
		st = &bernState{numbers: map[int]MPF{0: Fone}, m: 2, bin: I(10), bin1: I(1)}
		bernoulliCache[wp] = st
	}
	f3, f6 := FromInt64(3), FromInt64(6)
	for st.m <= n {
		m := st.m
		szbm := bernoulliSize(m)
		s := I(0)
		sexp := maxInt(0, szbm) - wp
		var a *big.Int
		if m < 6 {
			a = I(0)
		} else {
			a = st.bin1
		}
		for j := 1; j <= m/6; j++ {
			u := st.numbers[m-6*j]
			uman := u.Man
			if u.Sign == 1 {
				uman = neg(uman)
			}
			s = add(s, lsh(mul(a, uman), u.Exp-sexp))
			j6 := int64(6 * j)
			mm := int64(m)
			num := I(1)
			for _, f := range []int64{mm - 5 - j6, mm - 4 - j6, mm - 3 - j6, mm - 2 - j6, mm - 1 - j6, mm - j6} {
				num = muli(num, f)
			}
			den := I(1)
			for _, f := range []int64{4 + j6, 5 + j6, 6 + j6, 7 + j6, 8 + j6, 9 + j6} {
				den = muli(den, f)
			}
			a = floorDiv(mul(a, num), den)
		}
		var b MPF
		switch m % 6 {
		case 0, 2:
			b = MpfRdivInt(I(int64(m+3)), f3, wp, RoundFast)
		case 4:
			b = MpfRdivInt(I(int64(-m-3)), f6, wp, RoundFast)
		}
		sf := FromManExp(s, sexp, wp, RoundFast)
		b = MpfDiv(MpfSub(b, sf, wp, RoundFast), FromInt(st.bin, 0, RoundFast), wp, RoundFast)
		st.numbers[m] = b
		m += 2
		mm := int64(m)
		st.bin = floorDiv(mul(st.bin, I((mm+2)*(mm+3))), I(mm*(mm-1)))
		if m > 6 {
			st.bin1 = floorDiv(mul(st.bin1, I((2+mm)*(3+mm))), I((mm-7)*(mm-6)))
		}
		st.m = m
	}
	return MpfPos(st.numbers[n], prec, rnd)
}

func mpfBernoulliHuge(n int, prec int, rnd Rnd) MPF {
	wp := prec + 10
	piprec := wp + int(pyLog2Int(n))
	v := MpfGammaInt(I(int64(n+1)), wp, RoundFast)
	v = MpfMul(v, MpfZetaInt(n, wp, RoundFast), wp, RoundFast)
	v = MpfMul(v, MpfPowInt(MpfPi(piprec, RoundFast), I(int64(-n)), wp, RoundFast), 0, RoundFast)
	v = MpfShift(v, 1-n)
	if n&3 == 0 {
		v = MpfNeg(v, 0, RoundFast)
	}
	return MpfPos(v, prec, rnd)
}
