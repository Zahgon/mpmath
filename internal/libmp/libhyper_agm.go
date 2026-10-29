// Copyright (c) 2005-2026 Fredrik Johansson and mpmath contributors.
// Go port distributed under the same BSD 3-clause licence (see LICENSE).

package libmp

func mpcTupleEq(z, w MPC) bool { return z.Re.Eq(w.Re) && z.Im.Eq(w.Im) }

func MpfAgm(a, b MPF, prec int, rnd Rnd) MPF {
	if a.Sign != 0 || b.Sign != 0 {
		Raise("ComplexResult", "agm of a negative number")
	}
	if a.Man.Sign() == 0 || b.Man.Sign() == 0 {
		if a.Eq(Fnan) || b.Eq(Fnan) {
			return Fnan
		}
		if a.Eq(Finf) {
			if b.Eq(Fzero) {
				return Fnan
			}
			return Finf
		}
		if b.Eq(Finf) {
			if a.Eq(Fzero) {
				return Fnan
			}
			return Finf
		}
		return Fzero
	}
	wp := prec + 20
	amag := a.Exp + a.BC
	bmag := b.Exp + b.BC
	absMagDelta := absInt(amag - bmag)
	if absMagDelta > 10 {
		for absMagDelta > 10 {
			a, b = MpfShift(MpfAdd(a, b, wp, RoundFast), -1), MpfSqrt(MpfMul(a, b, wp, RoundFast), wp, RoundFast)
			absMagDelta /= 2
		}
		amag = a.Exp + a.BC
		bmag = b.Exp + b.BC
	}
	minMag := minInt(amag, bmag)
	maxMag := maxInt(amag, bmag)
	n := 0
	if minMag < -8 {
		n = -minMag
	} else if maxMag > 20 {
		n = -maxMag
	}
	if n != 0 {
		a = MpfShift(a, n)
		b = MpfShift(b, n)
	}
	g := AgmFixed(ToFixed(a, wp), ToFixed(b, wp), wp)
	return FromManExp(g, -wp-n, prec, rnd)
}

func MpfAgm1(a MPF, prec int, rnd Rnd) MPF { return MpfAgm(Fone, a, prec, rnd) }

func MpcAgm(a, b MPC, prec int, rnd Rnd) MPC {
	if MpcIsInfNan(a) || MpcIsInfNan(b) {
		return MPC{Fnan, Fnan}
	}
	if mpcTupleEq(a, MpcZero) || mpcTupleEq(b, MpcZero) {
		return MPC{Fzero, Fzero}
	}
	if mpcTupleEq(MpcNeg(a, 0, RoundFast), b) {
		return MPC{Fzero, Fzero}
	}
	wp := prec + 20
	eps := MpfShift(Fone, -wp+10)
	for {
		a1 := MpcShift(MpcAdd(a, b, wp, RoundFast), -1)
		b1 := MpcSqrt(MpcMul(a, b, wp, RoundFast), wp, RoundFast)
		a, b = a1, b1
		size := MpcAbs(a, 10, RoundFast)
		if s := MpcAbs(b, 10, RoundFast); MpfGt(s, size) {
			size = s
		}
		err := MpcAbs(MpcSub(a, b, 10, RoundFast), 10, RoundFast)
		if size.Eq(Fzero) || MpfLt(err, MpfMul(eps, size, 0, RoundFast)) {
			return a
		}
	}
}

func MpcAgm1(a MPC, prec int, rnd Rnd) MPC { return MpcAgm(MpcOne, a, prec, rnd) }

func MpfEllipk(x MPF, prec int, rnd Rnd) MPF {
	if x.Man.Sign() == 0 {
		if x.Eq(Fzero) {
			return MpfShift(MpfPi(prec, rnd), -1)
		}
		if x.Eq(Fninf) {
			return Fzero
		}
		if x.Eq(Fnan) {
			return x
		}
	}
	if x.Eq(Fone) {
		return Finf
	}
	wp := prec + 15
	a := MpfSqrt(MpfSub(Fone, x, wp, RoundFast), wp, RoundFast)
	v := MpfAgm1(a, wp, RoundFast)
	r := MpfDiv(MpfPi(wp, RoundFast), v, prec, rnd)
	return MpfShift(r, -1)
}

func MpcEllipk(z MPC, prec int, rnd Rnd) MPC {
	if z.Im.Eq(Fzero) {
		if z.Re.Eq(Finf) {
			return MpcZero
		}
		if MpfLe(z.Re, Fone) {
			return MPC{MpfEllipk(z.Re, prec, rnd), Fzero}
		}
	}
	wp := prec + 15
	a := MpcSqrt(MpcSub(MpcOne, z, wp, RoundFast), wp, RoundFast)
	v := MpcAgm1(a, wp, RoundFast)
	r := MpcMpfDiv(MpfPi(wp, RoundFast), v, prec, rnd)
	return MpcShift(r, -1)
}

func MpfEllipe(x MPF, prec int, rnd Rnd) MPF {
	if x.Man.Sign() == 0 {
		if x.Eq(Fzero) {
			return MpfShift(MpfPi(prec, rnd), -1)
		}
		if x.Eq(Fninf) {
			return Finf
		}
		if x.Eq(Fnan) {
			return x
		}
		if x.Eq(Finf) {
			Raise("ComplexResult", "")
		}
	}
	if x.Eq(Fone) {
		return Fone
	}
	wp := prec + 20
	mag := x.Exp + x.BC
	if mag < -wp {
		return MpfShift(MpfPi(prec, rnd), -1)
	}
	p := maxInt(mag, 0) - wp
	h := MpfShift(Fone, p)
	K := MpfEllipk(x, 2*wp, RoundFast)
	Kh := MpfEllipk(MpfSub(x, h, 0, RoundFast), 2*wp, RoundFast)
	Kdiff := MpfShift(MpfSub(K, Kh, 0, RoundFast), -p)
	t := MpfSub(Fone, x, 0, RoundFast)
	b := MpfMul(Kdiff, MpfShift(x, 1), wp, RoundFast)
	return MpfMul(t, MpfAdd(K, b, 0, RoundFast), prec, rnd)
}

func MpcEllipe(z MPC, prec int, rnd Rnd) MPC {
	if z.Im.Eq(Fzero) {
		if z.Re.Eq(Finf) {
			return MPC{Fzero, Finf}
		}
		if MpfLe(z.Re, Fone) {
			return MPC{MpfEllipe(z.Re, prec, rnd), Fzero}
		}
	}
	wp := prec + 15
	mag := MpcAbs(z, 1, RoundFast)
	p := maxInt(mag.Exp+mag.BC, 0) - wp
	h := MpfShift(Fone, p)
	K := MpcEllipk(z, 2*wp, RoundFast)
	Kh := MpcEllipk(MpcAddMpf(z, h, 2*wp, RoundFast), 2*wp, RoundFast)
	Kdiff := MpcShift(MpcSub(Kh, K, wp, RoundFast), -p)
	t := MpcSub(MpcOne, z, wp, RoundFast)
	b := MpcMul(Kdiff, MpcShift(z, 1), wp, RoundFast)
	return MpcMul(t, MpcAdd(K, b, wp, RoundFast), prec, rnd)
}
