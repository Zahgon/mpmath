// Copyright (c) 2005-2026 Fredrik Johansson and mpmath contributors.
// Go port distributed under the same BSD 3-clause licence (see LICENSE).

package libmp

import (
	"math"
	"math/big"
)

// Complex arithmetic from libmpc.py. An MPC is a (real, imag) pair.

type MPC struct{ Re, Im MPF }

var (
	MpcOne  = MPC{Fone, Fzero}
	MpcZero = MPC{Fzero, Fzero}
	MpcTwo  = MPC{Ftwo, Fzero}
	MpcHalf = MPC{Fhalf, Fzero}
)

func isInfMPF(x MPF) bool    { return x.Eq(Finf) || x.Eq(Fninf) }
func isInfNanMPF(x MPF) bool { return isInfMPF(x) || x.Eq(Fnan) }

// Eq is structural equality of the raw pairs.
func (z MPC) Eq(w MPC) bool { return z.Re.Eq(w.Re) && z.Im.Eq(w.Im) }

func MpcIsInf(z MPC) bool    { return isInfMPF(z.Re) || isInfMPF(z.Im) }
func MpcIsInfNan(z MPC) bool { return isInfNanMPF(z.Re) || isInfNanMPF(z.Im) }

// MpcToStr is mpc_to_str.
func MpcToStr(z MPC, dps int, o StrOpts) string {
	rs := ToStr(z.Re, dps, o)
	if z.Im.Sign != 0 {
		return rs + " - " + ToStr(MpfNeg(z.Im, 0, RoundFast), dps, o) + "j"
	}
	return rs + " + " + ToStr(z.Im, dps, o) + "j"
}

// MpcToComplex is mpc_to_complex.
func MpcToComplex(z MPC, strict bool, rnd Rnd) complex128 {
	return complex(ToFloat(z.Re, strict, rnd), ToFloat(z.Im, strict, rnd))
}

// MpcHash is mpc_hash before CPython reduces the returned int.
func MpcHash(z MPC) *big.Int {
	h := add(I(MpfHash(z.Re)), muli(I(MpfHash(z.Im)), 1000003))
	if h.Cmp(I(-1)) == 0 {
		return I(-2)
	}
	return h
}

func MpcConjugate(z MPC, prec int, rnd Rnd) MPC { return MPC{z.Re, MpfNeg(z.Im, prec, rnd)} }
func MpcIsNonzero(z MPC) bool                   { return !z.Eq(MpcZero) }

func MpcAdd(z, w MPC, prec int, rnd Rnd) MPC {
	return MPC{MpfAdd(z.Re, w.Re, prec, rnd), MpfAdd(z.Im, w.Im, prec, rnd)}
}
func MpcAddMpf(z MPC, x MPF, prec int, rnd Rnd) MPC { return MPC{MpfAdd(z.Re, x, prec, rnd), z.Im} }
func MpcSub(z, w MPC, prec int, rnd Rnd) MPC {
	return MPC{MpfSub(z.Re, w.Re, prec, rnd), MpfSub(z.Im, w.Im, prec, rnd)}
}
func MpcSubMpf(z MPC, p MPF, prec int, rnd Rnd) MPC { return MPC{MpfSub(z.Re, p, prec, rnd), z.Im} }
func MpcMpfSub(p MPF, z MPC, prec int, rnd Rnd) MPC {
	return MPC{MpfSub(p, z.Re, prec, rnd), MpfNeg(z.Im, prec, rnd)}
}
func MpcPos(z MPC, prec int, rnd Rnd) MPC {
	return MPC{MpfPos(z.Re, prec, rnd), MpfPos(z.Im, prec, rnd)}
}
func MpcNeg(z MPC, prec int, rnd Rnd) MPC {
	return MPC{MpfNeg(z.Re, prec, rnd), MpfNeg(z.Im, prec, rnd)}
}
func MpcShift(z MPC, n int) MPC           { return MPC{MpfShift(z.Re, n), MpfShift(z.Im, n)} }
func MpcAbs(z MPC, prec int, rnd Rnd) MPF { return MpfHypot(z.Re, z.Im, prec, rnd) }
func MpcArg(z MPC, prec int, rnd Rnd) MPF { return MpfAtan2(z.Im, z.Re, prec, rnd) }
func MpcFloor(z MPC, prec int, rnd Rnd) MPC {
	return MPC{MpfFloor(z.Re, prec, rnd), MpfFloor(z.Im, prec, rnd)}
}
func MpcCeil(z MPC, prec int, rnd Rnd) MPC {
	return MPC{MpfCeil(z.Re, prec, rnd), MpfCeil(z.Im, prec, rnd)}
}
func MpcNint(z MPC, prec int, rnd Rnd) MPC {
	return MPC{MpfNint(z.Re, prec, rnd), MpfNint(z.Im, prec, rnd)}
}
func MpcFrac(z MPC, prec int, rnd Rnd) MPC {
	return MPC{MpfFrac(z.Re, prec, rnd), MpfFrac(z.Im, prec, rnd)}
}

func exactMul(a, b MPF) MPF { return MpfMul(a, b, 0, RoundFast) }

// MpcMul is complex multiplication with each part rounded separately.
func MpcMul(z, w MPC, prec int, rnd Rnd) MPC {
	a, b, c, d := z.Re, z.Im, w.Re, w.Im
	p, q, r, s := exactMul(a, c), exactMul(b, d), exactMul(a, d), exactMul(b, c)
	return MPC{MpfSub(p, q, prec, rnd), MpfAdd(r, s, prec, rnd)}
}

func MpcSquare(z MPC, prec int, rnd Rnd) MPC {
	a, b := z.Re, z.Im
	p, q := exactMul(a, a), exactMul(b, b)
	r := MpfMul(a, b, prec, rnd)
	return MPC{MpfSub(p, q, prec, rnd), MpfShift(r, 1)}
}

func MpcMulMpf(z MPC, p MPF, prec int, rnd Rnd) MPC {
	return MPC{MpfMul(z.Re, p, prec, rnd), MpfMul(z.Im, p, prec, rnd)}
}

func MpcMulInt(z MPC, n *big.Int, prec int, rnd Rnd) MPC {
	return MPC{MpfMulInt(z.Re, n, prec, rnd), MpfMulInt(z.Im, n, prec, rnd)}
}

func MpcDiv(z, w MPC, prec int, rnd Rnd) MPC {
	if MpcIsInf(w) && !MpcIsInfNan(z) {
		return MpcZero
	}
	a, b, c, d := z.Re, z.Im, w.Re, w.Im
	wp := prec + 10
	mag := MpfAdd(exactMul(c, c), exactMul(d, d), wp, RoundFast)
	t := MpfAdd(exactMul(a, c), exactMul(b, d), wp, RoundFast)
	u := MpfSub(exactMul(b, c), exactMul(a, d), wp, RoundFast)
	return MPC{MpfDiv(t, mag, prec, rnd), MpfDiv(u, mag, prec, rnd)}
}

func MpcDivMpf(z MPC, p MPF, prec int, rnd Rnd) MPC {
	return MPC{MpfDiv(z.Re, p, prec, rnd), MpfDiv(z.Im, p, prec, rnd)}
}

func MpcReciprocal(z MPC, prec int, rnd Rnd) MPC {
	if MpcIsInf(z) {
		return MpcZero
	}
	a, b := z.Re, z.Im
	m := MpfAdd(exactMul(a, a), exactMul(b, b), prec+10, RoundFast)
	return MPC{MpfDiv(a, m, prec, rnd), MpfNeg(MpfDiv(b, m, prec, rnd), 0, RoundFast)}
}

func MpcMpfDiv(p MPF, z MPC, prec int, rnd Rnd) MPC {
	if MpcIsInf(z) && !isInfNanMPF(p) {
		return MpcZero
	}
	a, b := z.Re, z.Im
	m := MpfAdd(exactMul(a, a), exactMul(b, b), prec+10, RoundFast)
	re := MpfDiv(exactMul(a, p), m, prec, rnd)
	im := MpfDiv(MpfNeg(exactMul(b, p), 0, RoundFast), m, prec, rnd)
	return MPC{re, im}
}

// ComplexIntPow is (a+b*I)**n exactly for n >= 0.
func ComplexIntPow(a, b *big.Int, n int64) (*big.Int, *big.Int) {
	wre, wim := I(1), I(0)
	for n != 0 {
		if n&1 != 0 {
			wre, wim = sub(mul(wre, a), mul(wim, b)), add(mul(wim, a), mul(wre, b))
			n--
		}
		a, b = sub(mul(a, a), mul(b, b)), muli(mul(a, b), 2)
		n /= 2
	}
	return wre, wim
}

func MpcPow(z, w MPC, prec int, rnd Rnd) MPC {
	if w.Im.Eq(Fzero) {
		return MpcPowMpf(z, w.Re, prec, rnd)
	}
	return MpcExp(MpcMul(MpcLn(z, prec+10, RoundFast), w, prec+10, RoundFast), prec, rnd)
}

func MpcPowMpf(z MPC, p MPF, prec int, rnd Rnd) MPC {
	if p.Exp >= 0 {
		n := lsh(p.Man, p.Exp)
		if p.Sign == 1 {
			n = neg(n)
		}
		return MpcPowInt(z, n, prec, rnd)
	}
	if p.Exp == -1 {
		sqrtz := MpcSqrt(z, prec+10, RoundFast)
		n := p.Man
		if p.Sign == 1 {
			n = neg(n)
		}
		return MpcPowInt(sqrtz, n, prec, rnd)
	}
	return MpcExp(MpcMulMpf(MpcLn(z, prec+10, RoundFast), p, prec+10, RoundFast), prec, rnd)
}

func MpcPowInt(z MPC, n *big.Int, prec int, rnd Rnd) MPC {
	a, b := z.Re, z.Im
	if b.Eq(Fzero) {
		return MPC{MpfPowInt(a, n, prec, rnd), Fzero}
	}
	if a.Eq(Fzero) {
		v := MpfPowInt(b, n, prec, rnd)
		switch pyMod(n, I(4)).Int64() {
		case 0:
			return MPC{v, Fzero}
		case 1:
			return MPC{Fzero, v}
		case 2:
			return MPC{MpfNeg(v, 0, RoundFast), Fzero}
		default:
			return MPC{Fzero, MpfNeg(v, 0, RoundFast)}
		}
	}
	if n.IsInt64() {
		switch n.Int64() {
		case 0:
			return MpcOne
		case 1:
			return MpcPos(z, prec, rnd)
		case 2:
			return MpcSquare(z, prec, rnd)
		case -1:
			return MpcReciprocal(z, prec, rnd)
		}
	}
	if n.Sign() < 0 {
		return MpcReciprocal(MpcPowInt(z, neg(n), prec+4, RoundFast), prec, rnd)
	}
	aman, aexp, abc := a.Man, a.Exp, a.BC
	bman, bexp, bbc := b.Man, b.Exp, b.BC
	if a.Sign == 1 {
		aman = neg(aman)
	}
	if b.Sign == 1 {
		bman = neg(bman)
	}
	de := aexp - bexp
	exactSize := mul(n, I(int64(absInt(de)+maxInt(abc, bbc))))
	if exactSize.Cmp(I(10000)) < 0 && minInt(abc, bbc) >= 0 {
		if de > 0 {
			aman = lsh(aman, de)
			aexp = bexp
		} else {
			bman = lsh(bman, -de)
			bexp = aexp
		}
		nn := n.Int64()
		re, im := ComplexIntPow(aman, bman, nn)
		return MPC{FromManExp(re, int(nn)*aexp, prec, rnd), FromManExp(im, int(nn)*bexp, prec, rnd)}
	}
	return MpcExp(MpcMulInt(MpcLn(z, prec+10, RoundFast), n, prec+10, RoundFast), prec, rnd)
}

// MpcSqrt is the principal complex square root.
func MpcSqrt(z MPC, prec int, rnd Rnd) MPC {
	a, b := z.Re, z.Im
	if b.Eq(Fzero) {
		if a.Eq(Fzero) {
			return z
		}
		if a.Sign == 1 {
			return MPC{Fzero, MpfSqrt(MpfNeg(a, 0, RoundFast), prec, rnd)}
		}
		return MPC{MpfSqrt(a, prec, rnd), Fzero}
	}
	if isInfMPF(b) {
		return MPC{Finf, b}
	}
	wp := prec + 20
	var re, im MPF
	if a.Sign == 0 {
		t := MpfAdd(MpcAbs(z, wp, RoundFast), a, wp, RoundFast)
		u := MpfShift(t, -1)
		re = MpfSqrt(u, prec, rnd)
		v := MpfShift(t, 1)
		w := MpfSqrt(v, wp, RoundFast)
		im = MpfDiv(b, w, prec, rnd)
	} else {
		t := MpfSub(MpcAbs(z, wp, RoundFast), a, wp, RoundFast)
		u := MpfShift(t, -1)
		im = MpfSqrt(u, prec, rnd)
		v := MpfShift(t, 1)
		w := MpfSqrt(v, wp, RoundFast)
		re = MpfDiv(b, w, prec, rnd)
		if b.Sign == 1 {
			re = MpfNeg(re, 0, RoundFast)
			im = MpfNeg(im, 0, RoundFast)
		}
	}
	return MPC{re, im}
}

// pyComplexPowReal emulates CPython's (re + im*1j) ** e for a non-integral
// real e; ok is false where CPython raises OverflowError.
func pyComplexPowReal(re, im, e float64) (float64, float64, bool) {
	if re == 0 && im == 0 {
		return 0, 0, true
	}
	vabs := PyHypot(re, im)
	l := PyPow(vabs, e)
	phase := PyAtan2(im, re) * e
	r, i := l*PyCos(phase), l*PySin(phase)
	if math.IsInf(r, 0) || math.IsInf(i, 0) {
		return 0, 0, false
	}
	return r, i, true
}

func bigToFloat(x *big.Int) (float64, bool) {
	f, _ := new(big.Float).SetInt(x).Float64()
	return f, !math.IsInf(f, 0)
}

func mpcNthrootFixed(a, b *big.Int, n int, prec int) (*big.Int, *big.Int) {
	start := 50
	a1 := rsh(a, prec-n*start)
	b1 := rsh(b, prec-n*start)
	var re, im *big.Int
	fa, oka := bigToFloat(a1)
	fb, okb := bigToFloat(b1)
	var rr, ri float64
	ok := oka && okb
	if ok {
		rr, ri, ok = pyComplexPowReal(fa, fb, 1.0/float64(n))
	}
	if ok {
		re, _ = new(big.Float).SetFloat64(rr).Int(nil)
		im, _ = new(big.Float).SetFloat64(ri).Int(nil)
	} else {
		a1f := FromInt(a1, start, RoundFast)
		b1f := FromInt(b1, start, RoundFast)
		nth := MpfRdivInt(bigOne, FromInt64(int64(n)), start, RoundFast)
		r := MpcPow(MPC{a1f, b1f}, MPC{nth, Fzero}, start, RoundFast)
		re = ToInt(r.Re, RoundFast)
		im = ToInt(r.Im, RoundFast)
	}
	extra := 10
	prevp := start
	extra1 := n
	nb := I(int64(n))
	n1 := I(int64(n - 1))
	for _, p := range GiantSteps(start, prec+extra, 2) {
		re2, im2 := ComplexIntPow(re, im, int64(n-1))
		re2 = rsh(re2, (n-1)*prevp-p-extra1)
		im2 = rsh(im2, (n-1)*prevp-p-extra1)
		r4 := rsh(add(mul(re2, re2), mul(im2, im2)), p+extra1)
		ap := rsh(a, prec-p)
		bp := rsh(b, prec-p)
		rec := rsh(add(mul(ap, re2), mul(bp, im2)), p)
		imc := rsh(add(mul(neg(ap), im2), mul(bp, re2)), p)
		reb := floorDiv(lsh(rec, p), r4)
		imb := floorDiv(lsh(imc, p), r4)
		re = floorDiv(add(reb, mul(n1, lsh(re, p-prevp))), nb)
		im = floorDiv(add(imb, mul(n1, lsh(im, p-prevp))), nb)
		prevp = p
	}
	return re, im
}

// MpcNthroot is the complex n-th root.
func MpcNthroot(z MPC, n int, prec int, rnd Rnd) MPC {
	a, b := z.Re, z.Im
	if a.Sign == 0 && b.Eq(Fzero) {
		return MPC{MpfNthroot(a, n, prec, rnd), Fzero}
	}
	if n < 2 {
		switch n {
		case 0:
			return MpcOne
		case 1:
			return MpcPos(z, prec, rnd)
		case -1:
			return MpcDiv(MpcOne, z, prec, rnd)
		}
		inverse := MpcNthroot(z, -n, prec+5, ReciprocalRnd(rnd))
		return MpcDiv(MpcOne, inverse, prec, rnd)
	}
	if n <= 20 {
		prec2 := int(1.2 * float64(prec+10))
		pf := MpcAbs(z, prec, RoundFast)
		if pf.Exp+pf.BC > -10 && pf.Exp+pf.BC < prec {
			af := ToFixed(a, prec2)
			bf := ToFixed(b, prec2)
			re, im := mpcNthrootFixed(af, bf, n, prec2)
			extra := 10
			return MPC{FromManExp(re, -prec2-extra, prec2, rnd), FromManExp(im, -prec2-extra, prec2, rnd)}
		}
	}
	fn := FromInt64(int64(n))
	prec2 := prec + 10 + 10
	nth := MpfRdivInt(bigOne, fn, prec2, RoundFast)
	r := MpcPow(z, MPC{nth, Fzero}, prec2, rnd)
	re := Normalize(r.Re.Sign, r.Re.Man, r.Re.Exp, r.Re.BC, prec, rnd)
	im := Normalize(r.Im.Sign, r.Im.Man, r.Im.Exp, r.Im.BC, prec, rnd)
	return MPC{re, im}
}

func MpcCbrt(z MPC, prec int, rnd Rnd) MPC { return MpcNthroot(z, 3, prec, rnd) }

func MpcExp(z MPC, prec int, rnd Rnd) MPC {
	a, b := z.Re, z.Im
	if a.Eq(Fzero) {
		c, s := MpfCosSin(b, prec, rnd)
		return MPC{c, s}
	}
	if b.Eq(Fzero) {
		return MPC{MpfExp(a, prec, rnd), Fzero}
	}
	mag := MpfExp(a, prec+4, rnd)
	c, s := MpfCosSin(b, prec+4, rnd)
	return MPC{MpfMul(mag, c, prec, rnd), MpfMul(mag, s, prec, rnd)}
}

func MpcLn(z MPC, prec int, rnd Rnd) MPC {
	return MPC{MpfLogHypot(z.Re, z.Im, prec, rnd), MpcArg(z, prec, rnd)}
}

func MpcCos(z MPC, prec int, rnd Rnd) MPC {
	a, b := z.Re, z.Im
	if b.Eq(Fzero) {
		return MPC{MpfCos(a, prec, rnd), Fzero}
	}
	if a.Eq(Fzero) {
		return MPC{MpfCosh(b, prec, rnd), Fzero}
	}
	wp := prec + 6
	c, s := MpfCosSin(a, wp, RoundFast)
	ch, sh := MpfCoshSinh(b, wp, RoundFast)
	return MPC{MpfMul(c, ch, prec, rnd), MpfNeg(MpfMul(s, sh, prec, rnd), 0, RoundFast)}
}

func MpcSin(z MPC, prec int, rnd Rnd) MPC {
	a, b := z.Re, z.Im
	if b.Eq(Fzero) {
		return MPC{MpfSin(a, prec, rnd), Fzero}
	}
	if a.Eq(Fzero) {
		return MPC{Fzero, MpfSinh(b, prec, rnd)}
	}
	wp := prec + 6
	c, s := MpfCosSin(a, wp, RoundFast)
	ch, sh := MpfCoshSinh(b, wp, RoundFast)
	return MPC{MpfMul(s, ch, prec, rnd), MpfMul(c, sh, prec, rnd)}
}

func MpcTan(z MPC, prec int, rnd Rnd) MPC {
	a, b := z.Re, z.Im
	if b.Eq(Finf) {
		return MPC{Fzero, Fone}
	}
	if b.Eq(Fninf) {
		return MPC{Fzero, Fnone}
	}
	if isInfMPF(a) {
		if b.Eq(Fzero) {
			return MPC{Fnan, Fzero}
		}
		return MPC{Fnan, Fnan}
	}
	if b.Eq(Fzero) {
		return MPC{MpfTan(a, prec, rnd), Fzero}
	}
	if a.Eq(Fzero) {
		return MPC{Fzero, MpfTanh(b, prec, rnd)}
	}
	wp := prec + 15
	a = MpfShift(a, 1)
	b = MpfShift(b, 1)
	c, s := MpfCosSin(a, wp, RoundFast)
	ch, sh := MpfCoshSinh(b, wp, RoundFast)
	mag := MpfAdd(c, ch, wp, RoundFast)
	return MPC{MpfDiv(s, mag, prec, rnd), MpfDiv(sh, mag, prec, rnd)}
}

func MpcCosPi(z MPC, prec int, rnd Rnd) MPC {
	a, b := z.Re, z.Im
	if b.Eq(Fzero) {
		return MPC{MpfCosPi(a, prec, rnd), Fzero}
	}
	b = MpfMul(b, MpfPi(prec+5, RoundFast), prec+5, RoundFast)
	if a.Eq(Fzero) {
		return MPC{MpfCosh(b, prec, rnd), Fzero}
	}
	wp := prec + 6
	c, s := MpfCosSinPi(a, wp, RoundFast)
	ch, sh := MpfCoshSinh(b, wp, RoundFast)
	return MPC{MpfMul(c, ch, prec, rnd), MpfNeg(MpfMul(s, sh, prec, rnd), 0, RoundFast)}
}

func MpcSinPi(z MPC, prec int, rnd Rnd) MPC {
	a, b := z.Re, z.Im
	if b.Eq(Fzero) {
		return MPC{MpfSinPi(a, prec, rnd), Fzero}
	}
	b = MpfMul(b, MpfPi(prec+5, RoundFast), prec+5, RoundFast)
	if a.Eq(Fzero) {
		return MPC{Fzero, MpfSinh(b, prec, rnd)}
	}
	wp := prec + 6
	c, s := MpfCosSinPi(a, wp, RoundFast)
	ch, sh := MpfCoshSinh(b, wp, RoundFast)
	return MPC{MpfMul(s, ch, prec, rnd), MpfMul(c, sh, prec, rnd)}
}

func MpcCosSin(z MPC, prec int, rnd Rnd) (MPC, MPC) {
	a, b := z.Re, z.Im
	if a.Eq(Fzero) {
		ch, sh := MpfCoshSinh(b, prec, rnd)
		return MPC{ch, Fzero}, MPC{Fzero, sh}
	}
	if b.Eq(Fzero) {
		c, s := MpfCosSin(a, prec, rnd)
		return MPC{c, Fzero}, MPC{s, Fzero}
	}
	wp := prec + 6
	c, s := MpfCosSin(a, wp, RoundFast)
	ch, sh := MpfCoshSinh(b, wp, RoundFast)
	cre, cim := MpfMul(c, ch, prec, rnd), MpfMul(s, sh, prec, rnd)
	sre, sim := MpfMul(s, ch, prec, rnd), MpfMul(c, sh, prec, rnd)
	return MPC{cre, MpfNeg(cim, 0, RoundFast)}, MPC{sre, sim}
}

func MpcCosSinPi(z MPC, prec int, rnd Rnd) (MPC, MPC) {
	a, b := z.Re, z.Im
	if b.Eq(Fzero) {
		c, s := MpfCosSinPi(a, prec, rnd)
		return MPC{c, Fzero}, MPC{s, Fzero}
	}
	b = MpfMul(b, MpfPi(prec+5, RoundFast), prec+5, RoundFast)
	if a.Eq(Fzero) {
		ch, sh := MpfCoshSinh(b, prec, rnd)
		return MPC{ch, Fzero}, MPC{Fzero, sh}
	}
	wp := prec + 6
	c, s := MpfCosSinPi(a, wp, RoundFast)
	ch, sh := MpfCoshSinh(b, wp, RoundFast)
	cre, cim := MpfMul(c, ch, prec, rnd), MpfMul(s, sh, prec, rnd)
	sre, sim := MpfMul(s, ch, prec, rnd), MpfMul(c, sh, prec, rnd)
	return MPC{cre, MpfNeg(cim, 0, RoundFast)}, MPC{sre, sim}
}

func MpcCosh(z MPC, prec int, rnd Rnd) MPC {
	return MpcCos(MPC{z.Im, MpfNeg(z.Re, 0, RoundFast)}, prec, rnd)
}

func MpcSinh(z MPC, prec int, rnd Rnd) MPC {
	r := MpcSin(MPC{z.Im, z.Re}, prec, rnd)
	return MPC{r.Im, r.Re}
}

func MpcTanh(z MPC, prec int, rnd Rnd) MPC {
	r := MpcTan(MPC{z.Im, z.Re}, prec, rnd)
	return MPC{r.Im, r.Re}
}

func MpcAtan(z MPC, prec int, rnd Rnd) MPC {
	a, b := z.Re, z.Im
	wp := prec + 15
	x := MPC{MpfAdd(Fone, b, wp, RoundFast), MpfNeg(a, 0, RoundFast)}
	y := MPC{MpfSub(Fone, b, wp, RoundFast), a}
	l1 := MpcLn(x, wp, RoundFast)
	l2 := MpcLn(y, wp, RoundFast)
	d := MpcSub(l1, l2, prec, rnd)
	v := MPC{MpfNeg(MpfShift(d.Im, -1), 0, RoundFast), MpfShift(d.Re, -1)}
	if v.Im.Eq(Fnan) && MpcIsInf(z) {
		v.Im = Fzero
	}
	return v
}

var (
	betaCrossover  = FromFloat(0.6417, 53, RoundFast)
	alphaCrossover = FromFloat(1.5, 53, RoundFast)
)

func acosAsin(z MPC, prec int, rnd Rnd, n int) MPC {
	a, b := z.Re, z.Im
	wp := prec + 10
	if b.Eq(Fzero) {
		am := MpfSub(Fone, MpfAbs(a, 0, RoundFast), wp, RoundFast)
		if am.Sign == 0 {
			if n == 0 {
				return MPC{MpfAcos(a, prec, rnd), Fzero}
			}
			return MPC{MpfAsin(a, prec, rnd), Fzero}
		}
		if a.Sign == 1 {
			pi := MpfPi(prec, rnd)
			c := MpfAcosh(MpfNeg(a, 0, RoundFast), prec, rnd)
			if n == 0 {
				return MPC{pi, MpfNeg(c, 0, RoundFast)}
			}
			return MPC{MpfNeg(MpfShift(pi, -1), 0, RoundFast), c}
		}
		c := MpfAcosh(a, prec, rnd)
		if n == 0 {
			return MPC{Fzero, c}
		}
		pi := MpfPi(prec, rnd)
		return MPC{MpfShift(pi, -1), MpfNeg(c, 0, RoundFast)}
	}
	asign, bsign := 0, 0
	if a.Sign == 1 {
		a = MpfNeg(a, 0, RoundFast)
		asign = 1
	}
	if b.Sign == 1 {
		b = MpfNeg(b, 0, RoundFast)
		bsign = 1
	}
	f := RoundFast
	am := MpfSub(Fone, a, wp, f)
	ap := MpfAdd(Fone, a, wp, f)
	r := MpfHypot(ap, b, wp, f)
	s := MpfHypot(am, b, wp, f)
	alpha := MpfShift(MpfAdd(r, s, wp, f), -1)
	beta := MpfDiv(a, alpha, wp, f)
	b2 := MpfMul(b, b, wp, f)
	var re, im MPF
	if MpfSub(betaCrossover, beta, wp, f).Sign == 0 {
		if n == 0 {
			re = MpfAcos(beta, wp, f)
		} else {
			re = MpfAsin(beta, wp, f)
		}
	} else {
		Ax := MpfAdd(alpha, a, wp, f)
		if am.Sign == 0 {
			c := MpfDiv(b2, MpfAdd(r, ap, wp, f), wp, f)
			d := MpfAdd(s, am, wp, f)
			re = MpfShift(MpfMul(Ax, MpfAdd(c, d, wp, f), wp, f), -1)
			if n == 0 {
				re = MpfAtan(MpfDiv(MpfSqrt(re, wp, f), a, wp, f), wp, f)
			} else {
				re = MpfAtan(MpfDiv(a, MpfSqrt(re, wp, f), wp, f), wp, f)
			}
		} else {
			c := MpfDiv(Ax, MpfAdd(r, ap, wp, f), wp, f)
			d := MpfDiv(Ax, MpfSub(s, am, wp, f), wp, f)
			re = MpfShift(MpfAdd(c, d, wp, f), -1)
			re = MpfMul(b, MpfSqrt(re, wp, f), wp, f)
			if n == 0 {
				re = MpfAtan(MpfDiv(re, a, wp, f), wp, f)
			} else {
				re = MpfAtan(MpfDiv(a, re, wp, f), wp, f)
			}
		}
	}
	if MpfSub(alphaCrossover, alpha, wp, f).Sign == 0 {
		c1 := MpfDiv(b2, MpfAdd(r, ap, wp, f), wp, f)
		var Am1 MPF
		if MpfNeg(am, 0, f).Sign == 1 {
			c2 := MpfAdd(s, am, wp, f)
			c2 = MpfDiv(b2, c2, wp, f)
			Am1 = MpfShift(MpfAdd(c1, c2, wp, f), -1)
		} else {
			c2 := MpfSub(s, am, wp, f)
			Am1 = MpfShift(MpfAdd(c1, c2, wp, f), -1)
		}
		im = MpfMul(Am1, MpfAdd(alpha, Fone, wp, f), wp, f)
		im = MpfLog1p(MpfAdd(Am1, MpfSqrt(im, wp, f), wp, f), wp, f)
	} else {
		im = MpfSqrt(MpfSub(MpfMul(alpha, alpha, wp, f), Fone, wp, f), wp, f)
		im = MpfLn(MpfAdd(alpha, im, wp, f), wp, f)
	}
	if asign == 1 {
		if n == 0 {
			re = MpfSub(MpfPi(wp, f), re, wp, f)
		} else {
			re = MpfNeg(re, 0, f)
		}
	}
	if bsign == 0 && n == 0 {
		im = MpfNeg(im, 0, f)
	}
	if bsign == 1 && n == 1 {
		im = MpfNeg(im, 0, f)
	}
	if re.BC >= 0 {
		re = Normalize(re.Sign, re.Man, re.Exp, re.BC, prec, rnd)
	}
	if im.BC >= 0 {
		im = Normalize(im.Sign, im.Man, im.Exp, im.BC, prec, rnd)
	}
	if re.Eq(Fnan) && MpcIsInf(z) {
		a, b := z.Re, z.Im
		if isInfMPF(a) {
			if isInfMPF(b) {
				re = MpfShift(MpfPi(prec, rnd), -2)
				if a.Eq(Fninf) {
					if n == 0 {
						re = MpfMulInt(re, I(3), prec, rnd)
					} else {
						re = MpfNeg(re, 0, f)
					}
				}
			} else if b.Eq(Fnan) {
				if n == 0 {
					im = Finf
				} else {
					im = Fninf
				}
			} else {
				if n == 0 {
					if a.Eq(Finf) {
						re = Fzero
					} else {
						re = MpfPi(prec, rnd)
					}
				} else {
					re = MpfShift(MpfPi(prec, rnd), -1)
					if a.Eq(Fninf) {
						re = MpfNeg(re, 0, f)
					}
				}
			}
		} else {
			if n == 0 {
				return MPC{Fnan, MpfNeg(b, 0, f)}
			}
			return MPC{Fnan, b}
		}
	}
	return MPC{re, im}
}

func MpcAcos(z MPC, prec int, rnd Rnd) MPC { return acosAsin(z, prec, rnd, 0) }
func MpcAsin(z MPC, prec int, rnd Rnd) MPC { return acosAsin(z, prec, rnd, 1) }

func MpcAsinh(z MPC, prec int, rnd Rnd) MPC {
	r := MpcAsin(MPC{z.Im, MpfNeg(z.Re, 0, RoundFast)}, prec, rnd)
	return MPC{MpfNeg(r.Im, 0, RoundFast), r.Re}
}

func MpcAcosh(z MPC, prec int, rnd Rnd) MPC {
	r := MpcAcos(z, prec, rnd)
	if r.Im.Sign == 1 || r.Im.Eq(Fzero) {
		return MPC{MpfNeg(r.Im, 0, RoundFast), r.Re}
	}
	return MPC{r.Im, MpfNeg(r.Re, 0, RoundFast)}
}

func MpcAtanh(z MPC, prec int, rnd Rnd) MPC {
	wp := prec + 15
	a := MpcAdd(z, MpcOne, wp, RoundFast)
	b := MpcSub(MpcOne, z, wp, RoundFast)
	a = MpcLn(a, wp, RoundFast)
	b = MpcLn(b, wp, RoundFast)
	v := MpcShift(MpcSub(a, b, wp, RoundFast), -1)
	if v.Re.Eq(Fnan) && MpcIsInf(z) {
		v.Re = Fzero
	}
	return v
}

func MpcFibonacci(z MPC, prec int, rnd Rnd) MPC {
	if z.Im.Eq(Fzero) {
		return MPC{MpfFibonacci(z.Re, prec, rnd), Fzero}
	}
	size := absInt(z.Re.Exp + z.Re.BC)
	wp := prec + size + 20
	a := MpfPhi(wp, RoundFast)
	b := MpfAdd(MpfShift(a, 1), Fnone, wp, RoundFast)
	u := MpcPow(MPC{a, Fzero}, z, wp, RoundFast)
	v := MpcCosPi(z, wp, RoundFast)
	v = MpcDiv(v, u, wp, RoundFast)
	u = MpcSub(u, v, wp, RoundFast)
	return MpcDivMpf(u, b, prec, rnd)
}

func MpfExpj(x MPF, prec int, rnd Rnd) MPF {
	Raise("ComplexResult", "")
	return MPF{}
}

func MpcExpj(z MPC, prec int, rnd Rnd) MPC {
	re, im := z.Re, z.Im
	if im.Eq(Fzero) {
		c, s := MpfCosSin(re, prec, rnd)
		return MPC{c, s}
	}
	if re.Eq(Fzero) {
		return MPC{MpfExp(MpfNeg(im, 0, RoundFast), prec, rnd), Fzero}
	}
	ey := MpfExp(MpfNeg(im, 0, RoundFast), prec+10, RoundFast)
	c, s := MpfCosSin(re, prec+10, RoundFast)
	return MPC{MpfMul(ey, c, prec, rnd), MpfMul(ey, s, prec, rnd)}
}

func MpfExpjpi(x MPF, prec int, rnd Rnd) MPF {
	Raise("ComplexResult", "")
	return MPF{}
}

func MpcExpjpi(z MPC, prec int, rnd Rnd) MPC {
	re, im := z.Re, z.Im
	if im.Eq(Fzero) {
		c, s := MpfCosSinPi(re, prec, rnd)
		return MPC{c, s}
	}
	wp := prec + 10
	if im.Man.Sign() != 0 {
		wp += maxInt(0, im.Exp+im.BC)
	}
	im = MpfNeg(MpfMul(MpfPi(wp, RoundFast), im, wp, RoundFast), 0, RoundFast)
	if re.Eq(Fzero) {
		return MPC{MpfExp(im, prec, rnd), Fzero}
	}
	ey := MpfExp(im, prec+10, RoundFast)
	c, s := MpfCosSinPi(re, prec+10, RoundFast)
	return MPC{MpfMul(ey, c, prec, rnd), MpfMul(ey, s, prec, rnd)}
}
