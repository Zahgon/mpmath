// Copyright (c) 2005-2026 Fredrik Johansson and mpmath contributors.
// Go port distributed under the same BSD 3-clause licence (see LICENSE).

package libmp

import (
	"math"
	"math/big"
	"strings"
	"sync"
)

// Port of mpmath/libmp/libintmath.py with the pure-Python backend on
// CPython 3.12, where isqrt is math.isqrt (exact floor square root).

// GiantSteps returns the precision steps of a quadratically (n = 2) or more
// quickly convergent iteration from start to target.
func GiantSteps(start, target, n int) []int {
	L := []int{target}
	for L[len(L)-1] > start*n {
		L = append(L, floorDivInt(L[len(L)-1], n)+2)
	}
	for i, j := 0, len(L)-1; i < j; i, j = i+1, j-1 {
		L[i], L[j] = L[j], L[i]
	}
	return L
}

// Trailing counts the trailing zero bits of abs(n).
func Trailing(n *big.Int) int {
	if n.Sign() == 0 {
		return 0
	}
	a := absI(n)
	return int(a.TrailingZeroBits())
}

// BinToRadix converts x * 2**xbits to floor(x * base**bdigits).
func BinToRadix(x *big.Int, xbits, base, bdigits int) *big.Int {
	return rsh(mul(x, powi(I(int64(base)), bdigits)), xbits)
}

const stddigits = "0123456789abcdefghijklmnopqrstuvwxyz"

// Numeral is the string numeral of n in the given base.
func Numeral(n *big.Int, base int) string {
	return n.Text(base)
}

// Isqrt is the floor square root of a nonnegative integer.
func Isqrt(x *big.Int) *big.Int {
	return new(big.Int).Sqrt(x)
}

// Sqrtrem returns floor(sqrt(x)) and the remainder x - y*y.
func Sqrtrem(x *big.Int) (*big.Int, *big.Int) {
	y := new(big.Int).Sqrt(x)
	return y, sub(x, mul(y, y))
}

// SqrtFixed is isqrt(x << prec).
func SqrtFixed(x *big.Int, prec int) *big.Int {
	return Isqrt(lsh(x, prec))
}

// Gcd is math.gcd.
func Gcd(a, b *big.Int) *big.Int {
	return new(big.Int).GCD(nil, nil, absI(a), absI(b))
}

var ifibCache sync.Map

// Ifib is the nth Fibonacci number for integer n.
func Ifib(n int) *big.Int {
	if n < 0 {
		v := Ifib(-n)
		if (-n+1)%2 != 0 {
			return neg(v)
		}
		return v
	}
	if v, ok := ifibCache.Load(n); ok {
		return v.(*big.Int)
	}
	m := n
	a, b, p, q := I(1), I(0), I(0), I(1)
	for n != 0 {
		if n&1 == 1 {
			aq := mul(a, q)
			a, b = add(add(mul(b, q), aq), mul(a, p)), add(mul(b, p), aq)
			n--
		} else {
			qq := mul(q, q)
			p, q = add(mul(p, p), qq), add(qq, mul(muli(p, 2), q))
			n >>= 1
		}
	}
	ifibCache.Store(m, b)
	return b
}

var (
	ifacMu    sync.Mutex
	ifacCache = map[int]*big.Int{}
)

// Ifac is n! for an integer n >= 0.
func Ifac(n int) *big.Int {
	if n < 0 {
		Raise("ValueError", "factorial() not defined for negative values")
	}
	ifacMu.Lock()
	v, ok := ifacCache[n]
	ifacMu.Unlock()
	if ok {
		return v
	}
	v = new(big.Int).MulRange(1, int64(n))
	if n <= 1 {
		v = I(1)
	}
	ifacMu.Lock()
	ifacCache[n] = v
	ifacMu.Unlock()
	return v
}

// Ifac2 is n!! for an integer n >= 0.
func Ifac2(n int) *big.Int {
	p := I(1)
	for k := 2 - n&1; k <= n; k += 2 {
		p = muli(p, int64(k))
	}
	return p
}

// ListPrimes returns the primes up to n.
func ListPrimes(n int) []int {
	n++
	sieve := make([]int, n)
	for i := range sieve {
		sieve[i] = i
	}
	if n > 0 {
		sieve[0] = 0
	}
	if n > 1 {
		sieve[1] = 0
	}
	lim := int(math.Sqrt(float64(n))) + 1
	for i := 2; i < lim; i++ {
		if sieve[i] != 0 {
			for j := i * i; j < n; j += i {
				sieve[j] = 0
			}
		}
	}
	var out []int
	for _, p := range sieve {
		if p != 0 {
			out = append(out, p)
		}
	}
	return out
}

var smallOddPrimes = []int64{3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47}

// Isprime is the original's primality test.
func Isprime(n *big.Int) bool {
	if !isOdd(n) {
		return n.Cmp(bigTwo) == 0
	}
	if n.Cmp(I(50)) < 0 {
		for _, p := range smallOddPrimes {
			if n.Cmp(I(p)) == 0 {
				return true
			}
		}
		return false
	}
	for _, p := range smallOddPrimes {
		if pyMod(n, I(p)).Sign() == 0 {
			return false
		}
	}
	m := sub(n, bigOne)
	s := Trailing(m)
	d := rsh(m, s)
	test := func(a int64) bool {
		x := new(big.Int).Exp(I(a), d, n)
		if x.Cmp(bigOne) == 0 || x.Cmp(m) == 0 {
			return true
		}
		for r := 1; r < s; r++ {
			x = pyMod(mul(x, x), n)
			if x.Cmp(m) == 0 {
				return true
			}
		}
		return false
	}
	var witnesses []int64
	switch {
	case n.Cmp(I(1373653)) < 0:
		witnesses = []int64{2, 3}
	case n.Cmp(I(341550071728321)) < 0:
		witnesses = []int64{2, 3, 5, 7, 11, 13, 17}
	default:
		witnesses = smallOddPrimes
	}
	for _, a := range witnesses {
		if !test(a) {
			return false
		}
	}
	return true
}

// Moebius is the Moebius function of n.
func Moebius(n int) int {
	n = absInt(n)
	if n < 2 {
		return n
	}
	var factors []int
	for p := 2; p <= n; p++ {
		if n%p == 0 {
			if n%(p*p) == 0 {
				return 0
			}
			s := 0
			for _, f := range factors {
				s += p % f
			}
			if s == 0 {
				factors = append(factors, p)
			}
		}
	}
	if len(factors)%2 == 1 {
		return -1
	}
	return 1
}

// Eulernum is the Euler number E(m).
func Eulernum(m int) *big.Int {
	if m&1 == 1 {
		return I(0)
	}
	a := []*big.Int{I(0), I(0), I(1), I(0), I(0), I(0)}
	suma := I(1)
	n := m
	for i := 1; i <= m; i++ {
		for j := i + 1; j > -1; j -= 2 {
			a[j+1] = add(muli(a[j], int64(j-1)), muli(a[j+2], int64(j+1)))
		}
		a = append(a, I(0))
		suma = I(0)
		for k := i + 1; k > -1; k -= 2 {
			suma = add(suma, a[k+1])
		}
	}
	if (n/2)%2 == 1 {
		suma = neg(suma)
	}
	return floorDiv(suma, lsh(bigOne, n))
}

// Stirling1 is the Stirling number of the first kind.
func Stirling1(n, k int) *big.Int {
	if n < 0 || k < 0 {
		Raise("ValueError", "")
	}
	if k >= n {
		if n == k {
			return I(1)
		}
		return I(0)
	}
	if k < 1 {
		return I(0)
	}
	L := make([]*big.Int, k+1)
	for i := range L {
		L[i] = I(0)
	}
	L[1] = I(1)
	for m := 2; m <= n; m++ {
		for j := minInt(k, m); j > 0; j-- {
			L[j] = add(muli(L[j], int64(m-1)), L[j-1])
		}
	}
	if (n+k)%2 == 1 {
		return neg(L[k])
	}
	return L[k]
}

// Stirling2 is the Stirling number of the second kind.
func Stirling2(n, k int) *big.Int {
	if n < 0 || k < 0 {
		Raise("ValueError", "")
	}
	if k >= n {
		if n == k {
			return I(1)
		}
		return I(0)
	}
	if k <= 1 {
		if k == 1 {
			return I(1)
		}
		return I(0)
	}
	s := I(0)
	t := I(1)
	for j := 0; j <= k; j++ {
		term := mul(t, powi(I(int64(j)), n))
		if (k+j)&1 == 1 {
			s = sub(s, term)
		} else {
			s = add(s, term)
		}
		t = floorDiv(muli(t, int64(k-j)), I(int64(j+1)))
	}
	return floorDiv(s, Ifac(k))
}

// padLeft is str.rjust(n, c).
func padLeft(s string, n int, c byte) string {
	if len(s) >= n {
		return s
	}
	return strings.Repeat(string(c), n-len(s)) + s
}
