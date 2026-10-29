// Copyright (c) 2005-2026 Fredrik Johansson and mpmath contributors.
// Go port distributed under the same BSD 3-clause licence (see LICENSE).

// Package libmp is a port of mpmath's libmp: raw binary floating-point values
// (sign, mantissa, exponent, bit count) and the correctly rounded arithmetic,
// elementary and special functions built on them.
package libmp

import "math/big"

// Integer helpers with Python semantics. Every helper returns a new value and
// never modifies its arguments, so a *big.Int held by an MPF is immutable.

var (
	bigZero  = big.NewInt(0)
	bigOne   = big.NewInt(1)
	bigTwo   = big.NewInt(2)
	bigThree = big.NewInt(3)
	bigFive  = big.NewInt(5)
	bigTen   = big.NewInt(10)
)

// I returns n as a new big integer.
func I(n int64) *big.Int { return big.NewInt(n) }

func add(a, b *big.Int) *big.Int { return new(big.Int).Add(a, b) }
func sub(a, b *big.Int) *big.Int { return new(big.Int).Sub(a, b) }
func mul(a, b *big.Int) *big.Int { return new(big.Int).Mul(a, b) }
func neg(a *big.Int) *big.Int    { return new(big.Int).Neg(a) }
func absI(a *big.Int) *big.Int   { return new(big.Int).Abs(a) }

func addi(a *big.Int, n int64) *big.Int { return new(big.Int).Add(a, big.NewInt(n)) }
func muli(a *big.Int, n int64) *big.Int { return new(big.Int).Mul(a, big.NewInt(n)) }

// lsh is a << n; a negative n shifts right with floor rounding.
func lsh(a *big.Int, n int) *big.Int {
	if n < 0 {
		return new(big.Int).Rsh(a, uint(-n))
	}
	return new(big.Int).Lsh(a, uint(n))
}

// rsh is a >> n (floor rounding); a negative n shifts left.
func rsh(a *big.Int, n int) *big.Int {
	if n < 0 {
		return new(big.Int).Lsh(a, uint(-n))
	}
	return new(big.Int).Rsh(a, uint(n))
}

// floorDiv is Python's a // b.
func floorDiv(a, b *big.Int) *big.Int {
	if b.Sign() == 0 {
		Raise("ZeroDivisionError", "integer division or modulo by zero")
	}
	q, m := new(big.Int).QuoRem(a, b, new(big.Int))
	if m.Sign() != 0 && (m.Sign() < 0) != (b.Sign() < 0) {
		q.Sub(q, bigOne)
	}
	return q
}

// pyMod is Python's a % b.
func pyMod(a, b *big.Int) *big.Int {
	if b.Sign() == 0 {
		Raise("ZeroDivisionError", "integer modulo by zero")
	}
	m := new(big.Int).Rem(a, b)
	if m.Sign() != 0 && (m.Sign() < 0) != (b.Sign() < 0) {
		m.Add(m, b)
	}
	return m
}

// divmod is Python's divmod(a, b).
func divmod(a, b *big.Int) (*big.Int, *big.Int) {
	if b.Sign() == 0 {
		Raise("ZeroDivisionError", "integer division or modulo by zero")
	}
	q, m := new(big.Int).QuoRem(a, b, new(big.Int))
	if m.Sign() != 0 && (m.Sign() < 0) != (b.Sign() < 0) {
		q.Sub(q, bigOne)
		m.Add(m, b)
	}
	return q, m
}

func cmpI(a, b *big.Int) int       { return a.Cmp(b) }
func cmpi(a *big.Int, n int64) int { return a.Cmp(big.NewInt(n)) }
func isZero(a *big.Int) bool       { return a.Sign() == 0 }
func isOdd(a *big.Int) bool        { return a.Bit(0) == 1 }
func bitLen(a *big.Int) int        { return a.BitLen() }
func powi(a *big.Int, n int) *big.Int {
	return new(big.Int).Exp(a, big.NewInt(int64(n)), nil)
}
func pow10(n int) *big.Int { return powi(bigTen, n) }

// floorDivInt is Python's a // b for machine integers.
func floorDivInt(a, b int) int {
	q := a / b
	if (a%b != 0) && ((a < 0) != (b < 0)) {
		q--
	}
	return q
}

// pyModInt is Python's a % b for machine integers.
func pyModInt(a, b int) int {
	m := a % b
	if m != 0 && ((m < 0) != (b < 0)) {
		m += b
	}
	return m
}

func absInt(a int) int {
	if a < 0 {
		return -a
	}
	return a
}

// intBitLen is int.bit_length() for a machine integer.
func intBitLen(n int) int {
	return big.NewInt(int64(n)).BitLen()
}

func maxInt(a int, rest ...int) int {
	for _, b := range rest {
		if b > a {
			a = b
		}
	}
	return a
}

func minInt(a int, rest ...int) int {
	for _, b := range rest {
		if b < a {
			a = b
		}
	}
	return a
}
