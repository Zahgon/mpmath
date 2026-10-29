// Copyright (c) 2005-2026 Fredrik Johansson and mpmath contributors.
// Go port distributed under the same BSD 3-clause licence (see LICENSE).

package libmp

import "math/big"

// bernfrac, harmonic numbers and polygamma functions from gammazeta.py.

// pyLog2Big is math.log(q, 2) for a positive int q, following CPython's
// handling of ints too large for a float.
func pyLog2Big(q *big.Int) float64 {
	if q.BitLen() <= 1000 {
		f, _ := new(big.Float).SetInt(q).Float64()
		return PyLog(f) / PyLog(2)
	}
	e := q.BitLen()
	f, _ := new(big.Float).SetMantExp(new(big.Float).SetInt(q), -e).Float64()
	return (PyLog(f) + float64(PyLog(2)*float64(e))) / PyLog(2)
}

// Bernfrac is B_n as an exact reduced fraction p/q.
func Bernfrac(n int, plus bool) (*big.Int, *big.Int) {
	if n < 3 {
		switch n {
		case 0:
			return I(1), I(1)
		case 1:
			if plus {
				return I(1), I(2)
			}
			return I(-1), I(2)
		default:
			return I(1), I(6)
		}
	}
	if n&1 == 1 {
		return I(0), I(1)
	}
	q := I(1)
	for _, k := range ListPrimes(n + 1) {
		if n%(k-1) == 0 {
			q = muli(q, int64(k))
		}
	}
	prec := bernoulliSize(n) + int(pyLog2Big(q)) + 20
	b := MpfBernoulli(n, prec, RoundFast, false)
	p := MpfMul(b, FromInt(q, 0, RoundFast), 0, RoundFast)
	return ToInt(p, RoundNearest), q
}

func MpfHarmonic(x MPF, prec int, rnd Rnd) MPF {
	if x.Eq(Fzero) || x.Eq(Fnan) || x.Eq(Finf) {
		return x
	}
	a := MpfPsi0(MpfAdd(Fone, x, prec+5, RoundFast), prec, RoundFast)
	return MpfAdd(a, MpfEuler(prec+5, rnd), prec, rnd)
}

func MpcHarmonic(z MPC, prec int, rnd Rnd) MPC {
	if z.Im.Eq(Fzero) {
		return MPC{MpfHarmonic(z.Re, prec, rnd), Fzero}
	}
	a := MpcPsi0(MpcAddMpf(z, Fone, prec+5, RoundFast), prec, RoundFast)
	return MpcAddMpf(a, MpfEuler(prec+5, rnd), prec, rnd)
}

// MpfPsi0 is the digamma function of a real argument.
func MpfPsi0(x MPF, prec int, rnd Rnd) MPF {
	sign, man, exp, bc := x.Sign, x.Man, x.Exp, x.BC
	wp := prec + 10
	if man.Sign() == 0 {
		if x.Eq(Finf) {
			return x
		}
		if x.Eq(Fninf) || x.Eq(Fnan) {
			return Fnan
		}
	}
	if x.Eq(Fzero) || (exp >= 0 && sign == 1) {
		Raise("ValueError", "polygamma pole")
	}
	if exp+bc < -5 {
		v := MpfPsi0(MpfAdd(x, Fone, prec, rnd), prec, rnd)
		return MpfSub(v, MpfDiv(Fone, x, wp, rnd), prec, rnd)
	}
	if sign == 1 && exp+bc > 3 {
		c, s := MpfCosSinPi(x, wp, RoundFast)
		q := MpfMul(MpfDiv(c, s, wp, RoundFast), MpfPi(wp, RoundFast), wp, RoundFast)
		p := MpfPsi0(MpfSub(Fone, x, wp, RoundFast), wp, RoundFast)
		return MpfSub(p, q, prec, rnd)
	}
	if sign == 0 && bc+exp > wp {
		return MpfLn(MpfSub(x, Fone, wp, RoundFast), prec, rnd)
	}
	m := ToInt(x, RoundFast)
	n := int(0.11*float64(wp)) + 2
	s := I(0)
	xf := ToFixed(x, wp)
	one := lsh(bigOne, wp)
	if m.Cmp(I(int64(n))) < 0 {
		for k := m.Int64(); k < int64(n); k++ {
			s = sub(s, floorDiv(lsh(one, wp), xf))
			xf = add(xf, one)
		}
	}
	xf = sub(xf, one)
	s = add(s, ToFixed(MpfLn(FromManExp(xf, -wp, wp, RoundFast), wp, RoundFast), wp))
	s = add(s, floorDiv(lsh(one, wp), muli(xf, 2)))
	x2 := rsh(mul(xf, xf), wp)
	t := one
	prev := I(0)
	for k := 1; ; k++ {
		t = rsh(mul(t, x2), wp)
		b := MpfBernoulli(2*k, wp, RoundFast, false)
		offset := b.Exp + 2*wp
		var term *big.Int
		den := muli(t, int64(2*k))
		if offset >= 0 {
			term = floorDiv(lsh(b.Man, offset), den)
		} else {
			term = floorDiv(rsh(b.Man, -offset), den)
		}
		if k&1 == 1 {
			s = sub(s, term)
		} else {
			s = add(s, term)
		}
		if k > 2 && term.Cmp(prev) >= 0 {
			break
		}
		prev = term
	}
	return FromManExp(s, -wp, wp, rnd)
}

// MpcPsi0 is the digamma function of a complex argument.
func MpcPsi0(z MPC, prec int, rnd Rnd) MPC {
	re, im := z.Re, z.Im
	if im.Eq(Fzero) {
		return MPC{MpfPsi0(re, prec, rnd), Fzero}
	}
	wp := prec + 20
	sign, exp, bc := re.Sign, re.Exp, re.BC
	f := RoundFast
	if sign == 1 && exp+bc > 3 {
		c := MpcCosPi(z, wp, f)
		s := MpcSinPi(z, wp, f)
		q := MpcMulMpf(MpcDiv(c, s, wp, f), MpfPi(wp, f), wp, f)
		p := MpcPsi0(MpcSub(MpcOne, z, wp, f), wp, f)
		return MpcSub(p, q, prec, rnd)
	}
	if sign == 0 && bc+exp > wp {
		return MpcLn(MpcSub(z, MpcOne, wp, f), prec, rnd)
	}
	w := ToInt(re, RoundFast)
	n := int(0.11*float64(wp)) + 2
	s := MpcZero
	if w.Cmp(I(int64(n))) < 0 {
		for k := w.Int64(); k < int64(n); k++ {
			s = MpcSub(s, MpcReciprocal(z, wp, f), wp, f)
			z = MpcAddMpf(z, Fone, wp, f)
		}
	}
	z = MpcSub(z, MpcOne, wp, f)
	s = MpcAdd(s, MpcLn(z, wp, f), wp, f)
	s = MpcAdd(s, MpcDiv(MpcHalf, z, wp, f), wp, f)
	z2 := MpcSquare(z, wp, f)
	t := MpcOne
	szprev := Fzero
	eps := MpfShift(Fone, -wp+2)
	for k := 1; ; k++ {
		t = MpcMul(t, z2, wp, f)
		bern := MpfBernoulli(2*k, wp, f, false)
		term := MpcMpfDiv(bern, MpcMulInt(t, I(int64(2*k)), wp, f), wp, f)
		s = MpcSub(s, term, wp, f)
		szterm := MpcAbs(term, 10, f)
		if k > 2 && (MpfLe(szterm, eps) || MpfLe(szprev, szterm)) {
			break
		}
		szprev = szterm
	}
	return s
}

// MpfPsi is the polygamma function of integer order m >= 0.
func MpfPsi(m int, x MPF, prec int, rnd Rnd) MPF {
	if m == 0 {
		return MpfPsi0(x, prec, RoundFast)
	}
	return MpcPsi(m, MPC{x, Fzero}, prec, rnd).Re
}

// MpcPsi is the polygamma function of integer order m >= 0.
func MpcPsi(m int, z MPC, prec int, rnd Rnd) MPC {
	if m == 0 {
		return MpcPsi0(z, prec, rnd)
	}
	re, im := z.Re, z.Im
	wp := prec + 20
	f := RoundFast
	if im.Man.Sign() == 0 {
		if isInfNanMPF(im) {
			return MPC{Fnan, Fnan}
		}
	}
	if re.Man.Sign() == 0 {
		if re.Eq(Finf) && im.Eq(Fzero) {
			return MpcZero
		}
		if re.Eq(Fnan) {
			return MPC{Fnan, Fnan}
		}
	}
	w := ToInt(re, RoundFast)
	n := int(float64(0.4*float64(wp)) + float64(4*float64(m)))
	s := MpcZero
	if w.Cmp(I(int64(n))) < 0 {
		for k := w.Int64(); k < int64(n); k++ {
			t := MpcPowInt(z, I(int64(-m-1)), wp, f)
			s = MpcAdd(s, t, wp, f)
			z = MpcAddMpf(z, Fone, wp, f)
		}
	}
	zm := MpcPowInt(z, I(int64(-m)), wp, f)
	z2 := MpcPowInt(z, I(-2), wp, f)
	s = MpcAdd(s, MpcDivMpf(zm, FromInt64(int64(m)), wp, f), wp, f)
	s = MpcAdd(s, MpcMulMpf(MpcDiv(zm, z, wp, f), Fhalf, wp, f), wp, f)
	a := I(int64(m + 1))
	b := I(2)
	magn := MpcAbs(s, 10, f)
	eps := MpfShift(Fone, magn.Exp+magn.BC-wp+2)
	for k := int64(1); ; k++ {
		zm = MpcMul(zm, z2, wp, f)
		bern := MpfBernoulli(int(2*k), wp, f, false)
		scal := MpfMulInt(bern, a, wp, f)
		scal = MpfDiv(scal, FromInt(b, 0, f), wp, f)
		term := MpcMulMpf(zm, scal, wp, f)
		s = MpcAdd(s, term, wp, f)
		szterm := MpcAbs(term, 10, f)
		if k > 2 && MpfLe(szterm, eps) {
			break
		}
		mm := int64(m)
		a = mul(a, I((mm+2*k)*(mm+2*k+1)))
		b = mul(b, I((2*k+1)*(2*k+2)))
	}
	v := MpcMulMpf(s, MpfGamma(FromInt64(int64(m+1)), wp, f), prec, rnd)
	if m&1 == 0 {
		v = MPC{MpfNeg(v.Re, 0, f), MpfNeg(v.Im, 0, f)}
	}
	return v
}
