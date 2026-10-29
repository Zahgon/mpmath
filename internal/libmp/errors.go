// Copyright (c) 2005-2026 Fredrik Johansson and mpmath contributors.
// Go port distributed under the same BSD 3-clause licence (see LICENSE).

package libmp

import "fmt"

// PyError is a Python exception raised by the ported code: its class name and
// the text str() gives for it.
type PyError struct {
	Class string
	Msg   string
}

// Error renders the exception as "Class: message".
func (e *PyError) Error() string { return e.Class + ": " + e.Msg }

// Raise panics with a Python exception. The public package recovers it at its
// API boundary and returns it as an error.
func Raise(class, format string, args ...any) {
	msg := format
	if len(args) > 0 {
		msg = fmt.Sprintf(format, args...)
	}
	panic(&PyError{Class: class, Msg: msg})
}

// Is reports whether err is a raised exception of one of the given classes,
// following the original's class hierarchy.
func Is(err *PyError, classes ...string) bool {
	if err == nil {
		return false
	}
	for _, c := range classes {
		for k := err.Class; k != ""; k = parentClass[k] {
			if k == c {
				return true
			}
		}
	}
	return false
}

var parentClass = map[string]string{
	"ComplexResult":       "ValueError",
	"ValueError":          "Exception",
	"ZeroDivisionError":   "ArithmeticError",
	"OverflowError":       "ArithmeticError",
	"ArithmeticError":     "Exception",
	"NoConvergence":       "Exception",
	"TypeError":           "Exception",
	"NotImplementedError": "RuntimeError",
	"RuntimeError":        "Exception",
	"IndexError":          "LookupError",
	"KeyError":            "LookupError",
	"LookupError":         "Exception",
	"UnboundLocalError":   "NameError",
	"NameError":           "Exception",
}

// Try runs f and returns the Python exception it raised, if any. Any other
// panic is propagated unchanged.
func Try(f func()) (err *PyError) {
	defer func() {
		if r := recover(); r != nil {
			if pe, ok := r.(*PyError); ok {
				err = pe
				return
			}
			panic(r)
		}
	}()
	f()
	return nil
}
