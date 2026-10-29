// Copyright (c) 2005-2026 Fredrik Johansson and mpmath contributors.
// Go port distributed under the same BSD 3-clause licence (see LICENSE).

package libmp

import (
	"math"
	"math/big"
)

// pyIntToFloat is Python's float(int): round half to even, OverflowError past
// the largest double.
func pyIntToFloat(n *big.Int) float64 {
	f, _ := new(big.Float).SetInt(n).Float64()
	if math.IsInf(f, 0) {
		Raise("OverflowError", "int too large to convert to float")
	}
	return f
}

func MpfErf(x MPF, prec int, rnd Rnd) MPF {
	if x.Man.Sign() == 0 {
		switch {
		case x.Eq(Fzero):
			return Fzero
		case x.Eq(Finf):
			return Fone
		case x.Eq(Fninf):
			return Fnone
		}
		return Fnan
	}
	size := x.Exp + x.BC
	if size > 3 && float64(2*(size-1))+0.528766 > pyLog2Int(prec) {
		if x.Sign != 0 {
			return MpfPerturb(Fnone, 0, prec, rnd)
		}
		return MpfPerturb(Fone, 1, prec, rnd)
	}
	if size < -prec {
		x = MpfShift(x, 1)
		c := MpfSqrt(MpfPi(prec+20, RoundFast), prec+20, RoundFast)
		return MpfDiv(x, c, prec, rnd)
	}
	wp := prec + absInt(size) + 25
	t := absI(ToFixed(x, wp))
	t2 := rsh(mul(t, t), wp)
	s, term, k := t, I(12345), 1
	for term.Sign() != 0 {
		t = floorDiv(rsh(mul(t, t2), wp), I(int64(k)))
		term = floorDiv(t, I(int64(2*k+1)))
		if k&1 == 1 {
			s = sub(s, term)
		} else {
			s = add(s, term)
		}
		k++
	}
	s = floorDiv(lsh(s, wp+1), SqrtFixed(PiFixed(wp), wp))
	if x.Sign != 0 {
		s = neg(s)
	}
	return FromManExp(s, -wp, prec, rnd)
}

func erfcCheckSeries(x MPF, prec int) bool {
	n := ToInt(x, RoundFast)
	return pyIntToFloat(mul(n, n))*1.44 > float64(prec)
}

func MpfErfc(x MPF, prec int, rnd Rnd) MPF {
	if x.Man.Sign() == 0 {
		switch {
		case x.Eq(Fzero):
			return Fone
		case x.Eq(Finf):
			return Fzero
		case x.Eq(Fninf):
			return Ftwo
		}
		return Fnan
	}
	wp := prec + 20
	mag := x.BC + x.Exp
	wp += maxInt(0, 2*mag)
	regularErf := x.Sign != 0 || mag < 2
	if regularErf || !erfcCheckSeries(x, wp) {
		if regularErf {
			return MpfSub(Fone, MpfErf(x, prec+10, NegativeRnd(rnd)), prec, rnd)
		}
		n := addi(ToInt(x, RoundFast), 1)
		extra := int(pyIntToFloat(mul(n, n)) * 1.44)
		return MpfSub(Fone, MpfErf(x, prec+extra+10, RoundFast), prec, rnd)
	}
	s := lsh(bigOne, wp)
	term := s
	termPrev := I(0)
	fx := ToFixed(x, wp)
	t := rsh(muli(mul(fx, fx), 2), wp)
	for k := 1; ; k++ {
		term = floorDiv(lsh(muli(term, int64(2*k-1)), wp), t)
		if (k > 4 && term.Cmp(termPrev) > 0) || term.Sign() == 0 {
			break
		}
		if k&1 == 1 {
			s = sub(s, term)
		} else {
			s = add(s, term)
		}
		termPrev = term
	}
	s = floorDiv(lsh(s, wp), SqrtFixed(PiFixed(wp), wp))
	sf := FromManExp(s, -wp, wp, RoundFast)
	z := MpfExp(MpfNeg(MpfMul(x, x, wp, RoundFast), wp, RoundFast), wp, RoundFast)
	return MpfDiv(MpfMul(z, sf, wp, RoundFast), x, prec, rnd)
}

func eiTaylor(x *big.Int, prec int) *big.Int {
	s, t := x, x
	for k := int64(2); t.Sign() != 0; k++ {
		t = floorDiv(rsh(mul(t, x), prec), I(k))
		s = add(s, floorDiv(t, I(k)))
	}
	return s
}

func complexEiTaylor(zre, zim *big.Int, prec int) (*big.Int, *big.Int) {
	sre, tre := zre, zre
	sim, tim := zim, zim
	five := I(5)
	for k := int64(2); add(absI(tre), absI(tim)).Cmp(five) > 0; k++ {
		K := I(k)
		tre, tim = rsh(floorDiv(sub(mul(tre, zre), mul(tim, zim)), K), prec),
			rsh(floorDiv(add(mul(tre, zim), mul(tim, zre)), K), prec)
		sre = add(sre, floorDiv(tre, K))
		sim = add(sim, floorDiv(tim, K))
	}
	return sre, sim
}

func eiAsymptotic(x *big.Int, prec int) *big.Int {
	one := lsh(bigOne, prec)
	x = floorDiv(lsh(one, prec), x)
	t := x
	s := add(one, x)
	for k := int64(2); t.Sign() != 0; k++ {
		t = rsh(mul(muli(t, k), x), prec)
		s = add(s, t)
	}
	return s
}

func complexEiAsymptotic(zre, zim *big.Int, prec int) (*big.Int, *big.Int) {
	one := lsh(bigOne, prec)
	M := rsh(add(mul(zim, zim), mul(zre, zre)), prec)
	xre := floorDiv(lsh(zre, prec), M)
	xim := floorDiv(lsh(neg(zim), prec), M)
	tre, tim := xre, xim
	sre := add(one, xre)
	sim := xim
	thousand := I(1000)
	for k := int64(2); add(absI(tre), absI(tim)).Cmp(thousand) > 0; {
		tre, tim = rsh(muli(sub(mul(tre, xre), mul(tim, xim)), k), prec),
			rsh(muli(add(mul(tre, xim), mul(tim, xre)), k), prec)
		sre = add(sre, tre)
		sim = add(sim, tim)
		k++
		if k > int64(prec) {
			Raise("NoConvergence", "")
		}
	}
	return sre, sim
}

func mpfEi(x MPF, prec int, rnd Rnd, e1 bool) MPF {
	if e1 {
		x = MpfNeg(x, 0, RoundFast)
	}
	if e1 && x.Sign == 0 {
		if x.Eq(Fzero) {
			return Finf
		}
		Raise("ComplexResult", "E1(x) for x < 0")
	}
	var v MPF
	if x.Man.Sign() != 0 {
		xabs := MPF{0, x.Man, x.Exp, x.BC}
		xmag := x.Exp + x.BC
		wp := prec + 20
		canUseAsymp := xmag > wp
		if !canUseAsymp {
			var xabsint *big.Int
			if x.Exp >= 0 {
				xabsint = lsh(x.Man, x.Exp)
			} else {
				xabsint = rsh(x.Man, -x.Exp)
			}
			canUseAsymp = xabsint.Cmp(I(int64(int(float64(wp)*0.693)+10))) > 0
		}
		if canUseAsymp {
			if xmag > wp {
				v = Fone
			} else {
				v = FromManExp(eiAsymptotic(ToFixed(x, wp), wp), -wp, 0, RoundFast)
			}
			v = MpfMul(v, MpfExp(x, wp, RoundFast), wp, RoundFast)
			v = MpfDiv(v, x, prec, rnd)
		} else {
			wp += 2 * int(ToInt(xabs, RoundFast).Int64())
			u := ToFixed(x, wp)
			vi := add(eiTaylor(u, wp), EulerFixed(wp))
			t1 := FromManExp(vi, -wp, 0, RoundFast)
			t2 := MpfLn(xabs, wp, RoundFast)
			v = MpfAdd(t1, t2, prec, rnd)
		}
	} else {
		switch {
		case x.Eq(Fzero):
			v = Fninf
		case x.Eq(Finf):
			v = Finf
		case x.Eq(Fninf):
			v = Fzero
		default:
			v = Fnan
		}
	}
	if e1 {
		v = MpfNeg(v, 0, RoundFast)
	}
	return v
}

func mpcEi(z MPC, prec int, rnd Rnd, e1 bool) MPC {
	if e1 {
		z = MpcNeg(z, 0, RoundFast)
	}
	a, b := z.Re, z.Im
	if b.Eq(Fzero) {
		if e1 {
			x := MpfNeg(mpfEi(a, prec, rnd, false), 0, RoundFast)
			y := Fzero
			if a.Sign == 0 {
				y = MpfNeg(MpfPi(prec, rnd), 0, RoundFast)
			}
			return MPC{x, y}
		}
		return MPC{mpfEi(a, prec, rnd, false), Fzero}
	}
	if !a.Eq(Fzero) {
		if a.Man.Sign() == 0 || b.Man.Sign() == 0 {
			return MPC{Fnan, Fnan}
		}
	}
	wp := prec + 40
	amag := a.Exp + a.BC
	bmag := b.Exp + b.BC
	zmag := maxInt(amag, bmag)
	canUseAsymp := zmag > wp
	if !canUseAsymp {
		zabsint := add(absI(ToInt(a, RoundFast)), absI(ToInt(b, RoundFast)))
		canUseAsymp = zabsint.Cmp(I(int64(int(float64(wp)*0.693)+20))) > 0
	}
	if canUseAsymp {
		var v MPC
		var done bool
		err := Try(func() {
			if zmag > wp {
				v = MPC{Fone, Fzero}
			} else {
				zre := ToFixed(a, wp)
				zim := ToFixed(b, wp)
				vre, vim := complexEiAsymptotic(zre, zim, wp)
				v = MPC{FromManExp(vre, -wp, 0, RoundFast), FromManExp(vim, -wp, 0, RoundFast)}
			}
			v = MpcMul(v, MpcExp(z, wp, RoundFast), wp, RoundFast)
			v = MpcDiv(v, z, wp, RoundFast)
			if e1 {
				v = MpcNeg(v, prec, rnd)
			} else {
				x, y := v.Re, v.Im
				if b.Sign != 0 {
					v = MPC{MpfPos(x, prec, rnd), MpfSub(y, MpfPi(wp, RoundFast), prec, rnd)}
				} else {
					v = MPC{MpfPos(x, prec, rnd), MpfAdd(y, MpfPi(wp, RoundFast), prec, rnd)}
				}
			}
			done = true
		})
		if err != nil && err.Class != "NoConvergence" {
			panic(err)
		}
		if done {
			return v
		}
	}
	wp += 2 * int(ToInt(MpcAbs(z, 5, RoundFast), RoundFast).Int64())
	zre := ToFixed(a, wp)
	zim := ToFixed(b, wp)
	vre, vim := complexEiTaylor(zre, zim, wp)
	vre = add(vre, EulerFixed(wp))
	v := MPC{FromManExp(vre, -wp, 0, RoundFast), FromManExp(vim, -wp, 0, RoundFast)}
	var u MPC
	if e1 {
		u = MpcLn(MpcNeg(z, 0, RoundFast), wp, RoundFast)
	} else {
		u = MpcLn(z, wp, RoundFast)
	}
	v = MpcAdd(v, u, prec, rnd)
	if e1 {
		v = MpcNeg(v, 0, RoundFast)
	}
	return v
}

func MpfEi(x MPF, prec int, rnd Rnd) MPF { return mpfEi(x, prec, rnd, false) }
func MpcEi(z MPC, prec int, rnd Rnd) MPC { return mpcEi(z, prec, rnd, false) }
func MpfE1(x MPF, prec int, rnd Rnd) MPF { return mpfEi(x, prec, rnd, true) }
func MpcE1(z MPC, prec int, rnd Rnd) MPC { return mpcEi(z, prec, rnd, true) }
