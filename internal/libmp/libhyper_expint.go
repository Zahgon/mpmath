// Copyright (c) 2005-2026 Fredrik Johansson and mpmath contributors.
// Go port distributed under the same BSD 3-clause licence (see LICENSE).

package libmp

import "math/big"

// MpfExpint computes E_n(x), or the upper incomplete gamma function Gamma(n,x)
// when gamma is set. hasIm reports whether the optional branch-cut imaginary
// part im is present.
func MpfExpint(n int, x MPF, prec int, rnd Rnd, gamma bool) (re, im MPF, hasIm bool) {
	if x.Man.Sign() == 0 {
		if gamma {
			if x.Eq(Fzero) {
				if n <= 0 {
					return Finf, MPF{}, false
				}
				return MpfGammaInt(I(int64(n)), prec, rnd), MPF{}, false
			}
			if x.Eq(Finf) {
				return Fzero, MPF{}, false
			}
			return Fnan, Fnan, true
		}
		if x.Eq(Fzero) {
			if n > 1 {
				return FromRational(I(1), I(int64(n-1)), prec, rnd), MPF{}, false
			}
			return Finf, MPF{}, false
		}
		if x.Eq(Finf) {
			return Fzero, MPF{}, false
		}
		return Fnan, Fnan, true
	}
	nOrig := n
	if gamma {
		n = 1 - n
	}
	wp := prec + 20
	xmag := x.Exp + x.BC
	if xmag < -10 {
		Raise("NotImplementedError", "")
	}
	nmag := intBitLen(n)
	haveImag := n > 0 && x.Sign != 0
	negx := MpfNeg(x, 0, RoundFast)
	if n == 0 || 2*nmag-xmag < -wp {
		v := MpfExp(negx, wp, RoundFast)
		if gamma {
			re = MpfMul(v, MpfPowInt64(x, int64(nOrig-1), wp, RoundFast), prec, rnd)
		} else {
			re = MpfDiv(v, x, prec, rnd)
		}
	} else {
		canUseAsymp := -3*wp < n && n <= 0
		if !canUseAsymp {
			xi := absI(ToInt(x, RoundFast))
			d := sub(xi, I(int64(n)))
			m := 2 * wp
			if d.Cmp(I(int64(m))) < 0 {
				m = int(d.Int64())
			}
			m = minInt(maxInt(1, m), 2*wp)
			siz := -n*nmag + (m+n)*intBitLen(m+n) - m*xmag - floorDivInt(144*m, 100)
			tol := -wp - 10
			canUseAsymp = siz < tol
		}
		if canUseAsymp {
			r := floorDiv(lsh(I(-1), wp+wp), ToFixed(x, wp))
			m := int64(n)
			t := muli(r, m)
			s := lsh(bigOne, wp)
			for m != 0 && t.Sign() != 0 {
				s = add(s, t)
				m++
				t = rsh(mul(muli(r, m), t), wp)
			}
			v := MpfExp(negx, wp, RoundFast)
			if gamma {
				v = MpfMul(v, MpfPowInt64(x, int64(nOrig-1), wp, RoundFast), wp, RoundFast)
			} else {
				v = MpfDiv(v, x, wp, RoundFast)
			}
			re = MpfMul(v, FromManExp(s, -wp, 0, RoundFast), prec, rnd)
		} else if n == 1 {
			re = MpfNeg(MpfEi(negx, prec, rnd), 0, RoundFast)
		} else if n > 0 && n < 3*wp {
			T1 := MpfNeg(MpfEi(negx, wp, RoundFast), 0, RoundFast)
			if gamma {
				if nOrig&1 != 0 {
					T1 = MpfNeg(T1, 0, RoundFast)
				}
			} else {
				T1 = MpfMul(T1, MpfPowInt64(negx, int64(n-1), wp, RoundFast), wp, RoundFast)
			}
			r := ToFixed(x, wp)
			t := r
			facs := make([]*big.Int, n-1)
			facs[0] = I(1)
			for k := 1; k < n-1; k++ {
				facs[k] = muli(facs[k-1], int64(k))
			}
			for i, j := 0, len(facs)-1; i < j; i, j = i+1, j-1 {
				facs[i], facs[j] = facs[j], facs[i]
			}
			s := lsh(facs[0], wp)
			for k := 1; k < n-1; k++ {
				if k&1 == 1 {
					s = sub(s, mul(facs[k], t))
				} else {
					s = add(s, mul(facs[k], t))
				}
				t = rsh(mul(t, r), wp)
			}
			T2 := FromManExp(s, -wp, wp, RoundFast)
			T2 = MpfMul(T2, MpfExp(negx, wp, RoundFast), 0, RoundFast)
			if gamma {
				T2 = MpfMul(T2, MpfPowInt64(x, int64(nOrig), wp, RoundFast), wp, RoundFast)
			}
			R := MpfAdd(T1, T2, 0, RoundFast)
			re = MpfDiv(R, FromInt(Ifac(n-1), 0, RoundFast), prec, rnd)
		} else {
			Raise("NotImplementedError", "")
		}
	}
	if haveImag {
		M := FromInt(neg(Ifac(n-1)), 0, RoundFast)
		if gamma {
			im = MpfDiv(MpfPi(wp, RoundFast), M, prec, rnd)
			if nOrig&1 != 0 {
				im = MpfNeg(im, 0, RoundFast)
			}
		} else {
			im = MpfDiv(MpfMul(MpfPi(wp, RoundFast), MpfPowInt64(negx, int64(nOrig-1), wp, RoundFast), wp, RoundFast), M, prec, rnd)
		}
		return re, im, true
	}
	return re, MPF{}, false
}

func mpfCiSiTaylor(x MPF, wp int, which int) MPF {
	xf := ToFixed(x, wp)
	x2 := rsh(neg(mul(xf, xf)), wp)
	var s, t *big.Int
	var k int64
	if which == 0 {
		s, t, k = I(0), lsh(bigOne, wp), 2
	} else {
		s, t, k = xf, xf, 3
	}
	for t.Sign() != 0 {
		t = rsh(floorDiv(mul(t, x2), I(k*(k-1))), wp)
		s = add(s, floorDiv(t, I(k)))
		k += 2
	}
	return FromManExp(s, -wp, 0, RoundFast)
}

func mpcCiSiTaylor(re, im MPF, wp int, which int) MPC {
	var mag int
	haveMag := false
	if re.Man.Sign() != 0 {
		mag, haveMag = re.Exp+re.BC, true
	} else if im.Man.Sign() != 0 {
		mag, haveMag = im.Exp+im.BC, true
	}
	if im.Man.Sign() != 0 {
		mag = maxInt(mag, im.Exp+im.BC)
	}
	if !haveMag {
		Raise("UnboundLocalError", "cannot access local variable 'mag' where it is not associated with a value")
	}
	if mag > 2 || mag < -wp {
		Raise("NotImplementedError", "")
	}
	wp += 2 - mag
	zre := ToFixed(re, wp)
	zim := ToFixed(im, wp)
	z2re := rsh(sub(mul(zim, zim), mul(zre, zre)), wp)
	z2im := rsh(mul(muli(zre, -2), zim), wp)
	var sre, sim, tre, tim *big.Int
	var k int64
	if which == 0 {
		sre, sim, tre, tim, k = I(0), I(0), lsh(bigOne, wp), I(0), 2
	} else {
		sre, sim, tre, tim, k = zre, zim, zre, zim, 3
	}
	two := I(2)
	for {
		m := absI(tre)
		if a := absI(tim); a.Cmp(m) > 0 {
			m = a
		}
		if m.Cmp(two) <= 0 {
			break
		}
		f := I(k * (k - 1))
		tre, tim = rsh(floorDiv(sub(mul(tre, z2re), mul(tim, z2im)), f), wp),
			rsh(floorDiv(add(mul(tre, z2im), mul(tim, z2re)), f), wp)
		K := I(k)
		sre = add(sre, floorDiv(tre, K))
		sim = add(sim, floorDiv(tim, K))
		k += 2
	}
	return MPC{FromManExp(sre, -wp, 0, RoundFast), FromManExp(sim, -wp, 0, RoundFast)}
}

// MpfCiSi returns (Ci(x), Si(x)); which 0 computes only Ci, 1 only Si and 2
// both. A part that was not computed is reported absent through its flag.
func MpfCiSi(x MPF, prec int, rnd Rnd, which int) (ci MPF, hasCi bool, si MPF, hasSi bool) {
	wp := prec + 20
	if x.Man.Sign() == 0 {
		if x.Eq(Fzero) {
			return Fninf, true, Fzero, true
		}
		if x.Eq(Fnan) {
			return x, true, x, true
		}
		ci, hasCi = Fzero, true
		if which != 0 {
			if x.Eq(Finf) {
				si, hasSi = MpfShift(MpfPi(prec, rnd), -1), true
			}
			if x.Eq(Fninf) {
				si, hasSi = MpfNeg(MpfShift(MpfPi(prec, NegativeRnd(rnd)), -1), 0, RoundFast), true
			}
		}
		return
	}
	mag := x.Exp + x.BC
	if mag < -wp {
		if which != 0 {
			si, hasSi = MpfPerturb(x, 1-x.Sign, prec, rnd), true
		}
		if which != 1 {
			y := MpfEuler(wp, RoundFast)
			xabs := MpfAbs(x, 0, RoundFast)
			ci, hasCi = MpfAdd(y, MpfLn(xabs, wp, RoundFast), prec, rnd), true
		}
		return
	} else if mag > wp {
		if which != 0 {
			if x.Sign != 0 {
				si = MpfNeg(MpfPi(prec, NegativeRnd(rnd)), 0, RoundFast)
			} else {
				si = MpfPi(prec, rnd)
			}
			si, hasSi = MpfShift(si, -1), true
		}
		if which != 1 {
			ci, hasCi = MpfDiv(MpfSin(x, wp, RoundFast), x, prec, rnd), true
		}
		return
	} else {
		wp += absInt(mag)
	}
	asymptotic := float64(mag-1) > pyLog2Int(wp)
	if !asymptotic {
		if which != 0 {
			si, hasSi = MpfPos(mpfCiSiTaylor(x, wp, 1), prec, rnd), true
		}
		if which != 1 {
			c := mpfCiSiTaylor(x, wp, 0)
			c = MpfAdd(c, MpfEuler(wp, RoundFast), wp, RoundFast)
			ci, hasCi = MpfAdd(c, MpfLn(MpfAbs(x, 0, RoundFast), wp, RoundFast), prec, rnd), true
		}
		return
	}
	sign := x.Sign
	x = MpfAbs(x, 0, RoundFast)
	xf := ToFixed(x, wp)
	xr := floorDiv(lsh(bigOne, 2*wp), xf)
	s1 := lsh(bigOne, wp)
	s2 := xr
	t := xr
	k := int64(2)
	for t.Sign() != 0 {
		t = neg(t)
		t = rsh(muli(mul(t, xr), k), wp)
		k++
		s1 = add(s1, t)
		t = rsh(muli(mul(t, xr), k), wp)
		k++
		s2 = add(s2, t)
	}
	S1 := MpfDiv(FromManExp(s1, -wp, 0, RoundFast), x, wp, RoundFast)
	S2 := MpfDiv(FromManExp(s2, -wp, 0, RoundFast), x, wp, RoundFast)
	cos, sin := MpfCosSin(x, wp, RoundFast)
	if which != 0 {
		v := MpfAdd(MpfMul(cos, S1, 0, RoundFast), MpfMul(sin, S2, 0, RoundFast), wp, RoundFast)
		v = MpfSub(MpfShift(MpfPi(wp, RoundFast), -1), v, wp, RoundFast)
		if sign != 0 {
			v = MpfNeg(v, 0, RoundFast)
		}
		si, hasSi = MpfPos(v, prec, rnd), true
	}
	if which != 1 {
		ci, hasCi = MpfSub(MpfMul(sin, S1, 0, RoundFast), MpfMul(cos, S2, 0, RoundFast), prec, rnd), true
	}
	return
}

func MpfCi(x MPF, prec int, rnd Rnd) MPF {
	if MpfSign(x) < 0 {
		Raise("ComplexResult", "")
	}
	ci, _, _, _ := MpfCiSi(x, prec, rnd, 0)
	return ci
}

func MpfSi(x MPF, prec int, rnd Rnd) MPF {
	_, _, si, _ := MpfCiSi(x, prec, rnd, 1)
	return si
}

func MpcCi(z MPC, prec int, rnd Rnd) MPC {
	if z.Im.Eq(Fzero) {
		ci, _, _, _ := MpfCiSi(z.Re, prec, rnd, 0)
		if MpfSign(z.Re) < 0 {
			return MPC{ci, MpfPi(prec, rnd)}
		}
		return MPC{ci, Fzero}
	}
	wp := prec + 20
	c := mpcCiSiTaylor(z.Re, z.Im, wp, 0)
	c.Re = MpfAdd(c.Re, MpfEuler(wp, RoundFast), wp, RoundFast)
	return MpcAdd(c, MpcLn(z, wp, RoundFast), prec, rnd)
}

func MpcSi(z MPC, prec int, rnd Rnd) MPC {
	if z.Im.Eq(Fzero) {
		_, _, si, _ := MpfCiSi(z.Re, prec, rnd, 1)
		return MPC{si, Fzero}
	}
	wp := prec + 20
	return MpcPos(mpcCiSiTaylor(z.Re, z.Im, wp, 1), prec, rnd)
}

func MpfBesseljn(n int, x MPF, prec int, rnd Rnd) MPF {
	prec += 50
	negate := n < 0 && n&1 != 0
	mag := x.Exp + x.BC
	n = absInt(n)
	wp := prec + 20 + n*intBitLen(n)
	if mag < 0 {
		wp -= n * mag
	}
	xf := ToFixed(x, wp)
	x2 := rsh(mul(xf, xf), wp)
	var s, t *big.Int
	if n == 0 {
		s = lsh(bigOne, wp)
	} else {
		s = rsh(floorDiv(powi(xf, n), Ifac(n)), (n-1)*wp+n)
	}
	t = s
	for k := int64(1); t.Sign() != 0; k++ {
		t = rsh(floorDiv(mul(t, x2), I(-4*k*(k+int64(n)))), wp)
		s = add(s, t)
	}
	if negate {
		s = neg(s)
	}
	return FromManExp(s, -wp, prec, rnd)
}

func MpcBesseljn(n int, z MPC, prec int, rnd Rnd) MPC {
	negate := n < 0 && n&1 != 0
	n = absInt(n)
	origprec := prec
	mag := maxInt(z.Re.Exp+z.Re.BC, z.Im.Exp+z.Im.BC)
	prec += 20 + n*intBitLen(n) + absInt(mag)
	if mag < 0 {
		prec -= n * mag
	}
	zre := ToFixed(z.Re, prec)
	zim := ToFixed(z.Im, prec)
	z2re := rsh(sub(mul(zre, zre), mul(zim, zim)), prec)
	z2im := rsh(mul(zre, zim), prec-1)
	var sre, sim, tre, tim *big.Int
	if n == 0 {
		sre, tre = lsh(bigOne, prec), lsh(bigOne, prec)
		sim, tim = I(0), I(0)
	} else {
		re, im := ComplexIntPow(zre, zim, int64(n))
		sh := (n-1)*prec + n
		sre = rsh(floorDiv(re, Ifac(n)), sh)
		sim = rsh(floorDiv(im, Ifac(n)), sh)
		tre, tim = sre, sim
	}
	three := I(3)
	for k := int64(1); add(absI(tre), absI(tim)).Cmp(three) > 0; k++ {
		p := I(-4 * k * (k + int64(n)))
		tre, tim = sub(mul(tre, z2re), mul(tim, z2im)), add(mul(tim, z2re), mul(tre, z2im))
		tre = rsh(floorDiv(tre, p), prec)
		tim = rsh(floorDiv(tim, p), prec)
		sre = add(sre, tre)
		sim = add(sim, tim)
	}
	if negate {
		sre, sim = neg(sre), neg(sim)
	}
	return MPC{FromManExp(sre, -prec, origprec, rnd), FromManExp(sim, -prec, origprec, rnd)}
}
