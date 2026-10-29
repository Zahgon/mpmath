// Copyright (c) 2005-2026 Fredrik Johansson and mpmath contributors.
// Go port distributed under the same BSD 3-clause licence (see LICENSE).

package libmp

import (
	"math"
	"math/big"
)

// The gamma function and its variants from gammazeta.py. typ is 0 for
// gamma(x), 1 for x!, 2 for 1/gamma(x) and 3 for log|gamma(x)|.

// MpfGamma is mpf_gamma.
func MpfGamma(x MPF, prec int, rnd Rnd) MPF { return mpfGamma(x, prec, rnd, 0) }

func mpfGamma(x MPF, prec int, rnd Rnd, typ int) MPF {
	sign, man, exp, bc := x.Sign, x.Man, x.Exp, x.BC
	f := RoundFast
	if man.Sign() == 0 {
		if x.Eq(Fzero) {
			if typ == 1 {
				return Fone
			}
			if typ == 2 {
				return Fzero
			}
			Raise("ValueError", "gamma function pole")
		}
		if x.Eq(Finf) {
			if typ == 2 {
				return Fzero
			}
			return Finf
		}
		return Fnan
	}
	if typ == 3 {
		wp := prec + 20
		if exp+bc > wp && sign == 0 {
			return MpfSub(MpfMul(x, MpfLn(x, wp, f), wp, f), x, prec, rnd)
		}
	}
	isInteger := exp >= 0
	var n *big.Int
	if isInteger {
		if sign == 1 {
			if typ == 2 {
				return Fzero
			}
			Raise("ValueError", "gamma function pole")
		}
		n = lsh(man, exp)
		if n.Cmp(I(smallFactorialCacheSize)) < 0 {
			k := int(n.Int64())
			switch typ {
			case 0:
				return MpfPos(smallFactorialCache[k-1], prec, rnd)
			case 1:
				return MpfPos(smallFactorialCache[k], prec, rnd)
			case 2:
				return MpfDiv(Fone, smallFactorialCache[k-1], prec, rnd)
			case 3:
				return MpfLn(smallFactorialCache[k-1], prec, rnd)
			}
		}
	} else {
		n = rsh(man, -exp)
	}
	mag := exp + bc
	gammaSize := muli(n, int64(mag))
	var wp int
	if typ == 3 {
		wp = prec + 20
	} else {
		wp = prec + gammaSize.BitLen() + 20
	}
	if mag < -wp {
		switch typ {
		case 0:
			return MpfSub(MpfDiv(Fone, x, wp, f), MpfShift(Fone, -wp), prec, rnd)
		case 1:
			return MpfSub(Fone, x, prec, rnd)
		case 2:
			return MpfAdd(x, MpfShift(Fone, mag-wp), prec, rnd)
		case 3:
			return MpfNeg(MpfLn(MpfAbs(x, 0, f), prec, rnd), 0, f)
		}
	}
	if typ == 1 {
		return mpfGamma(MpfAdd(x, Fone, 0, f), prec, rnd, 0)
	}
	tenWp := I(int64(10 * wp))
	if exp >= -1 {
		if isInteger {
			if gammaSize.Cmp(tenWp) < 0 {
				fac := Ifac(int(n.Int64()) - 1)
				switch typ {
				case 0:
					return FromInt(fac, prec, rnd)
				case 2:
					return FromRational(I(1), fac, prec, rnd)
				case 3:
					return MpfLn(FromInt(fac, 0, f), prec, rnd)
				}
			}
		}
		if n.Cmp(I(100)) < 0 || gammaSize.Cmp(tenWp) < 0 {
			nn := int(n.Int64())
			if sign == 1 {
				w := SqrtpiFixed(wp)
				fv := Ifac2(2*nn + 1)
				if nn%2 == 0 {
					fv = neg(fv)
				}
				switch typ {
				case 0:
					return MpfShift(FromRational(w, fv, prec, rnd), -wp+nn+1)
				case 2:
					return MpfShift(FromRational(fv, w, prec, rnd), wp-nn-1)
				case 3:
					return MpfLn(MpfShift(FromRational(w, absI(fv), prec, rnd), -wp+nn+1), prec, rnd)
				}
			} else if nn == 0 {
				switch typ {
				case 0:
					return MpfSqrtpi(prec, rnd)
				case 2:
					return MpfDiv(Fone, MpfSqrtpi(wp, f), prec, rnd)
				case 3:
					return MpfLn(MpfSqrtpi(wp, f), prec, rnd)
				}
			} else {
				w := SqrtpiFixed(wp)
				wf := FromManExp(mul(w, Ifac2(2*nn-1)), -wp-nn, 0, f)
				switch typ {
				case 0:
					return MpfPos(wf, prec, rnd)
				case 2:
					return MpfDiv(Fone, wf, prec, rnd)
				case 3:
					return MpfLn(MpfAbs(wf, 0, f), prec, rnd)
				}
			}
		}
	}
	absxman := lsh(man, exp+wp)
	if typ == 3 && sign == 0 {
		one := lsh(bigOne, wp)
		oneDist := absI(sub(absxman, one))
		twoDist := absI(sub(absxman, muli(one, 2)))
		md := oneDist
		if twoDist.Cmp(md) < 0 {
			md = twoDist
		}
		cancellation := wp - md.BitLen()
		if cancellation > 10 {
			xsub1 := MpfSub(Fone, x, 0, f)
			xsub2 := MpfSub(Ftwo, x, 0, f)
			xsub1mag := xsub1.Exp + xsub1.BC
			xsub2mag := xsub2.Exp + xsub2.BC
			if xsub1mag < -wp {
				return MpfMul(MpfEuler(wp, f), MpfSub(Fone, x, 0, f), prec, rnd)
			}
			if xsub2mag < -wp {
				return MpfMul(MpfSub(Fone, MpfEuler(wp, f), 0, f), MpfSub(x, Ftwo, 0, f), prec, rnd)
			}
			wp += maxInt(-xsub1mag, -xsub2mag)
			absxman = lsh(man, exp+wp)
		}
	}
	nForStirling := int(gammaStirlingBeta * float64(wp))
	if n.Cmp(I(int64(maxInt(100, nForStirling)))) < 0 && wp < maxGammaTaylorPrec {
		if sign == 1 {
			absxman = neg(absxman)
		}
		return gammaFixedTaylor(x, absxman, wp, prec, rnd, typ)
	}
	xorig := x
	var r *big.Int
	var xabs MPF
	if n.Cmp(I(int64(nForStirling))) < 0 {
		one := lsh(bigOne, wp)
		r = one
		d := nForStirling - int(n.Int64())
		for k := 0; k < d; k++ {
			r = rsh(mul(r, absxman), wp)
			absxman = add(absxman, one)
		}
		xabs = FromManExp(absxman, -wp, 0, f)
	} else {
		xabs = MpfAbs(x, 0, f)
	}
	y := realStirlingSeries(absxman, wp)
	u := ToFixed(MpfLn(xabs, wp, f), wp)
	u = rsh(mul(sub(absxman, lsh(bigOne, wp-1)), u), wp)
	y = add(y, u)
	w := FromManExp(y, -wp, 0, f)
	hasR := r != nil && r.Sign() != 0
	if sign == 1 {
		A := MpfMul(MpfSinPi(xorig, wp, f), xorig, wp, f)
		B := MpfNeg(MpfPi(wp, f), 0, f)
		if typ == 0 || typ == 2 {
			A = MpfMul(A, MpfExp(w, wp, f), 0, f)
			if hasR {
				B = MpfMul(B, FromManExp(r, -wp, 0, f), wp, f)
			}
			if typ == 0 {
				return MpfDiv(B, A, prec, rnd)
			}
			return MpfDiv(A, B, prec, rnd)
		}
		if hasR {
			B = MpfMul(B, FromManExp(r, -wp, 0, f), wp, f)
		}
		A = MpfAdd(MpfLn(MpfAbs(A, 0, f), wp, f), w, wp, f)
		return MpfSub(MpfLn(MpfAbs(B, 0, f), wp, f), A, prec, rnd)
	}
	switch typ {
	case 0:
		if hasR {
			return MpfDiv(MpfExp(w, wp, f), FromManExp(r, -wp, 0, f), prec, rnd)
		}
		return MpfExp(w, prec, rnd)
	case 2:
		if hasR {
			return MpfDiv(FromManExp(r, -wp, 0, f), MpfExp(w, wp, f), prec, rnd)
		}
		return MpfExp(MpfNeg(w, 0, f), prec, rnd)
	default:
		if hasR {
			return MpfSub(w, MpfLn(FromManExp(r, -wp, 0, f), wp, f), prec, rnd)
		}
		return MpfPos(w, prec, rnd)
	}
}

// MpcGamma is mpc_gamma.
func MpcGamma(z MPC, prec int, rnd Rnd) MPC { return mpcGamma(z, prec, rnd, 0) }

func mpcGamma(z MPC, prec int, rnd Rnd, typ int) MPC {
	a, b := z.Re, z.Im
	f := RoundFast
	if b.Eq(Fzero) {
		if typ == 3 && a.Sign == 1 {
			re := mpfGamma(a, prec, rnd, 3)
			n := rsh(neg(a.Man), -a.Exp)
			im := MpfMulInt(MpfPi(prec+10, f), n, prec, rnd)
			return MPC{re, im}
		}
		return MPC{mpfGamma(a, prec, rnd, typ), Fzero}
	}
	if (a.Man.Sign() == 0 && a.Exp != 0) || (b.Man.Sign() == 0 && b.Exp != 0) {
		return MPC{Fnan, Fnan}
	}
	wp := prec + 20
	amag := a.Exp + a.BC
	bmag := b.Exp + b.BC
	mag := bmag
	if a.Man.Sign() != 0 {
		mag = maxInt(amag, bmag)
	}
	if mag < -8 {
		if mag < -wp {
			v := MpcAdd(z, MpcMulMpf(MpcMul(z, z, wp, f), MpfEuler(wp, f), wp, f), wp, f)
			switch typ {
			case 0:
				return MpcReciprocal(v, prec, rnd)
			case 1:
				return MpcDiv(z, v, prec, rnd)
			case 2:
				return MpcPos(v, prec, rnd)
			case 3:
				return MpcLn(MpcReciprocal(v, prec, f), prec, rnd)
			}
		} else if typ != 1 {
			wp += -mag
		}
	}
	if typ == 3 && mag > wp && (a.Sign == 0 || bmag >= amag) {
		return MpcSub(MpcMul(z, MpcLn(z, wp, f), wp, f), z, prec, rnd)
	}
	if typ == 1 {
		return mpcGamma(MPC{MpfAdd(a, Fone, 0, f), b}, prec, rnd, 0)
	}
	an := absI(ToInt(a, f))
	bn := absI(ToInt(b, f))
	absn := an
	if bn.Cmp(absn) > 0 {
		absn = bn
	}
	gammaSize := muli(absn, int64(mag))
	if typ != 3 {
		wp += gammaSize.BitLen()
	}
	needReflection := a.Sign == 1
	zorig := z
	if needReflection {
		z = MpcNeg(z, 0, f)
		a, b = z.Re, z.Im
	}
	var yfinal *MPC
	balancePrec := 0
	if bmag < -10 {
		if typ == 3 {
			zsub1 := MpcSubMpf(z, Fone, 0, f)
			var cancel1 int
			if zsub1.Re.Eq(Fzero) {
				cancel1 = -bmag
			} else {
				cancel1 = -maxInt(zsub1.Re.Exp+zsub1.Re.BC, bmag)
			}
			if cancel1 > wp {
				pi := MpfPi(wp, f)
				x := MpcMulMpf(zsub1, pi, wp, f)
				x = MpcMul(x, x, wp, f)
				x = MpcDivMpf(x, FromInt64(12), wp, f)
				y := MpcMulMpf(zsub1, MpfNeg(MpfEuler(wp, f), 0, f), wp, f)
				yf := MpcAdd(x, y, wp, f)
				yfinal = &yf
				if !needReflection {
					return MpcPos(yf, prec, rnd)
				}
			} else if cancel1 > 0 {
				wp += cancel1
			}
			zsub2 := MpcSubMpf(z, Ftwo, 0, f)
			var cancel2 int
			if zsub2.Re.Eq(Fzero) {
				cancel2 = -bmag
			} else {
				cancel2 = -maxInt(zsub2.Re.Exp+zsub2.Re.BC, bmag)
			}
			if cancel2 > wp {
				pi := MpfPi(wp, f)
				t := MpfSub(MpfMul(pi, pi, 0, f), FromInt64(6), 0, f)
				x := MpcMulMpf(MpcMul(zsub2, zsub2, wp, f), t, wp, f)
				x = MpcDivMpf(x, FromInt64(12), wp, f)
				y := MpcMulMpf(zsub2, MpfSub(Fone, MpfEuler(wp, f), 0, f), wp, f)
				yf := MpcAdd(x, y, wp, f)
				yfinal = &yf
				if !needReflection {
					return MpcPos(yf, prec, rnd)
				}
			} else if cancel2 > 0 {
				wp += cancel2
			}
		}
		if bmag < -wp {
			pp := 2 * (wp + 10)
			aabs := MpfAbs(a, 0, f)
			eps := MpfShift(Fone, amag-wp)
			x1 := mpfGamma(aabs, pp, f, typ)
			x2 := mpfGamma(MpfAdd(aabs, eps, 0, f), pp, f, typ)
			xprime := MpfDiv(MpfSub(x2, x1, pp, f), eps, pp, f)
			y := MpfMul(b, xprime, prec, rnd)
			yf := MPC{x1, y}
			yfinal = &yf
			if !needReflection {
				return MpcPos(yf, prec, rnd)
			}
		} else {
			balancePrec += -bmag
		}
	}
	wp += balancePrec
	nForStirling := int(gammaStirlingBeta * float64(wp))
	afix := ToFixed(a, wp)
	bfix := ToFixed(b, wp)
	var r *MPC
	var y MPC
	if yfinal == nil {
		zprered := z
		if absn.Cmp(I(int64(nForStirling))) < 0 {
			anf, _ := new(big.Float).SetInt(an).Float64()
			nfs := I(int64(nForStirling))
			base := sub(addi(mul(nfs, nfs), 1), mul(bn, bn))
			basef, _ := new(big.Float).SetInt(base).Float64()
			d := int(float64(PyPow(basef, 0.5)) - anf)
			one := lsh(bigOne, wp)
			rre, rim := one, I(0)
			for k := 0; k < d; k++ {
				rre, rim = rsh(sub(mul(afix, rre), mul(bfix, rim)), wp), rsh(add(mul(afix, rim), mul(bfix, rre)), wp)
				afix = add(afix, one)
			}
			rr := MPC{FromManExp(rre, -wp, 0, f), FromManExp(rim, -wp, 0, f)}
			r = &rr
			a = FromManExp(afix, -wp, 0, f)
			z = MPC{a, b}
		}
		yre, yim := complexStirlingSeries(afix, bfix, wp)
		l := MpcLn(z, wp, f)
		lre := ToFixed(l.Re, wp)
		lim := ToFixed(l.Im, wp)
		yre = add(sub(rsh(sub(mul(lre, afix), mul(lim, bfix)), wp), rsh(lre, 1)), yre)
		yim = add(sub(rsh(add(mul(lre, bfix), mul(lim, afix)), wp), rsh(lim, 1)), yim)
		y = MPC{FromManExp(yre, -wp, 0, f), FromManExp(yim, -wp, 0, f)}
		if r != nil && typ == 3 {
			y = MpcSub(y, MpcLn(*r, wp, f), wp, f)
			zfa := ToFloat(zprered.Re, false, f)
			zfb := ToFloat(zprered.Im, false, f)
			zfabs := PyHypot(zfa, zfb)
			yfb := ToFloat(y.Im, false, f)
			u := PyAtan2(zfb, zfa)
			var gi float64
			if zfabs <= 0.5 {
				gi = float64(0.577216*zfb) - u
			} else {
				gi = -zfb - float64(0.5*u)
				gi = gi + float64(zfa*u)
				gi = gi + float64(zfb*PyLog(zfabs))
			}
			n := int(math.Floor(float64((gi-yfb)/float64(2*math.Pi)) + 0.5))
			y = MPC{y.Re, MpfAdd(y.Im, MpfMulInt(MpfPi(wp, f), I(int64(2*n)), wp, f), wp, f)}
		}
	}
	if needReflection {
		if typ == 0 || typ == 2 {
			A := MpcMul(MpcSinPi(zorig, wp, f), zorig, wp, f)
			B := MPC{MpfNeg(MpfPi(wp, f), 0, f), Fzero}
			if yfinal != nil {
				if typ == 2 {
					A = MpcDiv(A, *yfinal, wp, f)
				} else {
					A = MpcMul(A, *yfinal, wp, f)
				}
			} else {
				A = MpcMul(A, MpcExp(y, wp, f), wp, f)
			}
			if r != nil {
				B = MpcMul(B, *r, wp, f)
			}
			if typ == 0 {
				return MpcDiv(B, A, prec, rnd)
			}
			return MpcDiv(A, B, prec, rnd)
		}
		var s1 MPC
		if yfinal != nil {
			s1 = MpcNeg(*yfinal, 0, f)
		} else {
			s1 = MpcNeg(y, 0, f)
		}
		s1 = MpcSub(s1, MpcLn(MpcNeg(zorig, 0, f), wp, f), wp, f)
		rezfloor := MpfFloor(zorig.Re, 0, f)
		imzsign := MpfSign(zorig.Im)
		pi := MpfPi(wp, f)
		t := MpfMul(pi, rezfloor, 0, f)
		t = MpfMulInt(t, I(int64(imzsign)), wp, f)
		s1 = MPC{s1.Re, MpfAdd(s1.Im, t, wp, f)}
		s1 = MpcAddMpf(s1, MpfLn(pi, wp, f), wp, f)
		tc := MpcSinPi(MpcSubMpf(zorig, rezfloor, 0, f), wp, f)
		tc = MpcLn(tc, wp, f)
		s1 = MpcSub(s1, tc, wp, f)
		if imzsign == 0 {
			t := MpfMul(pi, MpfFloor(rezfloor, 0, f), wp, f)
			s1 = MPC{s1.Re, MpfSub(s1.Im, t, wp, f)}
		}
		return MpcPos(s1, prec, rnd)
	}
	switch typ {
	case 0:
		if r != nil {
			return MpcDiv(MpcExp(y, wp, f), *r, prec, rnd)
		}
		return MpcExp(y, prec, rnd)
	case 2:
		if r != nil {
			return MpcDiv(*r, MpcExp(y, wp, f), prec, rnd)
		}
		return MpcExp(MpcNeg(y, 0, f), prec, rnd)
	default:
		return MpcPos(y, prec, rnd)
	}
}

func MpfFactorial(x MPF, prec int, rnd Rnd) MPF { return mpfGamma(x, prec, rnd, 1) }
func MpcFactorial(z MPC, prec int, rnd Rnd) MPC { return mpcGamma(z, prec, rnd, 1) }
func MpfRgamma(x MPF, prec int, rnd Rnd) MPF    { return mpfGamma(x, prec, rnd, 2) }
func MpcRgamma(z MPC, prec int, rnd Rnd) MPC    { return mpcGamma(z, prec, rnd, 2) }

func MpfLoggamma(x MPF, prec int, rnd Rnd) MPF {
	if x.Sign == 1 {
		Raise("ComplexResult", "")
	}
	return mpfGamma(x, prec, rnd, 3)
}

func MpcLoggamma(z MPC, prec int, rnd Rnd) MPC {
	a, b := z.Re, z.Im
	if b.Eq(Fzero) && a.Sign == 1 {
		re := mpfGamma(a, prec, rnd, 3)
		n := rsh(neg(a.Man), -a.Exp)
		im := MpfMulInt(MpfPi(prec+10, RoundFast), n, prec, rnd)
		return MPC{re, im}
	}
	return mpcGamma(z, prec, rnd, 3)
}

// MpfGammaInt is gamma(n) for a positive integer n.
func MpfGammaInt(n *big.Int, prec int, rnd Rnd) MPF {
	if n.Cmp(I(smallFactorialCacheSize)) < 0 {
		return MpfPos(smallFactorialCache[n.Int64()-1], prec, rnd)
	}
	return MpfGamma(FromInt(n, 0, RoundFast), prec, rnd)
}
