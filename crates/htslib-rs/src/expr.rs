//! HTSlib-style filter expression parsing and evaluation.

use regex::Regex;
use std::{error, fmt};

/// An expression value.
#[derive(Clone, Debug, PartialEq)]
pub struct Value {
    is_str: bool,
    is_true: bool,
    string: Option<String>,
    number: f64,
}

impl Value {
    /// Creates a numeric value.
    pub fn number(number: f64) -> Self {
        Self {
            is_str: false,
            is_true: false,
            string: None,
            number,
        }
    }

    /// Creates a string value.
    pub fn string(value: impl Into<String>) -> Self {
        Self {
            is_str: true,
            is_true: false,
            string: Some(value.into()),
            number: 0.0,
        }
    }

    /// Creates an undefined value.
    pub fn undefined() -> Self {
        Self {
            is_str: false,
            is_true: false,
            string: None,
            number: f64::NAN,
        }
    }

    /// Forces the value to evaluate as true.
    pub fn with_true(mut self) -> Self {
        self.is_true = true;
        self
    }

    /// Returns whether this value is a string.
    pub fn is_string(&self) -> bool {
        self.is_str
    }

    /// Returns whether this value is explicitly true.
    pub fn is_true(&self) -> bool {
        self.is_true
    }

    /// Returns the string value, if this is a defined string.
    pub fn as_str(&self) -> Option<&str> {
        self.string.as_deref()
    }

    /// Returns the numeric value.
    pub fn number_value(&self) -> f64 {
        self.number
    }

    /// Returns whether this value is defined.
    pub fn exists(&self) -> bool {
        if self.is_str {
            self.string.is_some()
        } else {
            !self.number.is_nan()
        }
    }

    /// Returns whether this value is defined or explicitly true.
    pub fn exists_true(&self) -> bool {
        self.is_true || self.exists()
    }

    /// Returns the HTSlib expression truth value.
    pub fn truth(&self) -> bool {
        self.is_true
            || (self.is_str && self.string.is_some())
            || (!self.is_str && self.exists() && self.number != 0.0)
    }

    fn boolean(value: bool) -> Self {
        Self {
            is_str: false,
            is_true: value,
            string: None,
            number: f64::from(value),
        }
    }

    fn set_undefined(&mut self) {
        *self = Self::undefined();
    }

    fn finalize(&mut self) {
        if self.is_str {
            if self.string.is_some() {
                self.is_true = true;
            }
            self.number = f64::from(self.is_true);
        } else if self.exists() && self.number != 0.0 {
            self.is_true = true;
        }
    }
}

/// A parsed HTSlib filter expression.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Filter {
    src: String,
}

impl Filter {
    /// Creates a filter expression.
    pub fn new(src: impl Into<String>) -> Self {
        Self { src: src.into() }
    }

    /// Evaluates the expression using a symbol lookup callback.
    pub fn eval_with<F>(&self, lookup: F) -> Result<Value, Error>
    where
        F: Fn(&str) -> Option<(Value, usize)>,
    {
        let mut parser = Parser {
            src: &self.src,
            pos: 0,
            lookup: &lookup,
        };

        let mut value = parser.parse_expression()?;
        parser.skip_ws();

        if !parser.is_eof() {
            return Err(Error::new(parser.pos, "trailing input"));
        }

        value.finalize();
        Ok(value)
    }
}

/// Expression parser error.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Error {
    position: usize,
    message: &'static str,
}

impl Error {
    fn new(position: usize, message: &'static str) -> Self {
        Self { position, message }
    }

    /// Returns the byte position where parsing failed.
    pub fn position(&self) -> usize {
        self.position
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} at byte {}", self.message, self.position)
    }
}

impl error::Error for Error {}

struct Parser<'a, F>
where
    F: Fn(&str) -> Option<(Value, usize)>,
{
    src: &'a str,
    pos: usize,
    lookup: &'a F,
}

impl<F> Parser<'_, F>
where
    F: Fn(&str) -> Option<(Value, usize)>,
{
    fn parse_expression(&mut self) -> Result<Value, Error> {
        self.parse_and_or()
    }

    fn parse_and_or(&mut self) -> Result<Value, Error> {
        let mut lhs = self.parse_eq()?;

        loop {
            self.skip_ws();

            if self.consume("&&") {
                let rhs = self.parse_eq()?;

                if !lhs.exists_true() || !rhs.exists_true() {
                    lhs = Value::number(0.0);
                    lhs.set_undefined();
                    lhs.number = 0.0;
                } else {
                    lhs = Value::boolean(lhs.truth() && rhs.truth());
                }
            } else if self.consume("||") {
                let rhs = self.parse_eq()?;

                if (!lhs.exists_true() && !rhs.exists_true())
                    || (!lhs.exists_true() && !rhs.truth())
                    || (!rhs.exists_true() && !lhs.truth())
                {
                    lhs = Value::number(0.0);
                    lhs.set_undefined();
                    lhs.number = 0.0;
                } else {
                    lhs = Value::boolean(lhs.truth() || rhs.truth());
                }
            } else {
                break;
            }
        }

        Ok(lhs)
    }

    fn parse_eq(&mut self) -> Result<Value, Error> {
        let mut lhs = self.parse_cmp()?;
        self.skip_ws();

        let Some(op) = self.consume_any(["==", "!=", "=~", "!~"]) else {
            return Ok(lhs);
        };

        let rhs = self.parse_eq()?;

        match op {
            "==" => {
                if !lhs.exists() || !rhs.exists() {
                    lhs.set_undefined();
                } else {
                    let eq = if lhs.is_str {
                        lhs.string.is_some() && rhs.string.is_some() && lhs.string == rhs.string
                    } else {
                        !rhs.is_str && lhs.number == rhs.number
                    };
                    lhs = Value::boolean(eq);
                }
            }
            "!=" => {
                if !lhs.exists() || !rhs.exists() {
                    lhs.set_undefined();
                } else {
                    let ne = if lhs.is_str {
                        if lhs.string.is_some() && rhs.string.is_some() {
                            lhs.string != rhs.string
                        } else {
                            true
                        }
                    } else {
                        rhs.is_str || lhs.number != rhs.number
                    };
                    lhs = Value::boolean(ne);
                }
            }
            "=~" | "!~" => {
                if !lhs.is_str || !rhs.is_str {
                    return Err(Error::new(self.pos, "regex operands must be strings"));
                }

                let matched = match (lhs.string.as_deref(), rhs.string.as_deref()) {
                    (Some(haystack), Some(pattern)) => Regex::new(pattern)
                        .map_err(|_| Error::new(self.pos, "invalid regex"))?
                        .is_match(haystack),
                    _ => false,
                };

                lhs = Value::boolean(if op == "=~" { matched } else { !matched });
            }
            _ => unreachable!(),
        }

        Ok(lhs)
    }

    fn parse_cmp(&mut self) -> Result<Value, Error> {
        let mut lhs = self.parse_bitor()?;
        self.skip_ws();

        let Some(op) = self.consume_any([">=", ">", "<=", "<"]) else {
            return Ok(lhs);
        };

        let rhs = self.parse_cmp()?;

        if !lhs.exists() || !rhs.exists() {
            lhs.set_undefined();
            return Ok(lhs);
        }

        let result = match (lhs.string.as_deref(), rhs.string.as_deref()) {
            (Some(a), Some(b)) if lhs.is_str && rhs.is_str => match op {
                ">=" => a >= b,
                ">" => a > b,
                "<=" => a <= b,
                "<" => a < b,
                _ => unreachable!(),
            },
            _ if !lhs.is_str && !rhs.is_str => match op {
                ">=" => lhs.number >= rhs.number,
                ">" => lhs.number > rhs.number,
                "<=" => lhs.number <= rhs.number,
                "<" => lhs.number < rhs.number,
                _ => unreachable!(),
            },
            _ => false,
        };

        Ok(Value::boolean(result))
    }

    fn parse_bitor(&mut self) -> Result<Value, Error> {
        self.parse_bitwise(Self::parse_bitxor, "|", |a, b| a | b)
    }

    fn parse_bitxor(&mut self) -> Result<Value, Error> {
        self.parse_bitwise(Self::parse_bitand, "^", |a, b| a ^ b)
    }

    fn parse_bitand(&mut self) -> Result<Value, Error> {
        self.parse_bitwise(Self::parse_add, "&", |a, b| a & b)
    }

    fn parse_bitwise(
        &mut self,
        next: fn(&mut Self) -> Result<Value, Error>,
        op: &str,
        f: fn(i64, i64) -> i64,
    ) -> Result<Value, Error> {
        let mut lhs = next(self)?;
        let mut undef = false;

        loop {
            self.skip_ws();

            if (op == "&" && self.starts_with("&&")) || (op == "|" && self.starts_with("||")) {
                break;
            }

            if !self.consume(op) {
                break;
            }

            let rhs = next(self)?;

            if !lhs.exists() || !rhs.exists() {
                undef = true;
            } else if lhs.is_str || rhs.is_str {
                return Err(Error::new(self.pos, "bitwise operands must be numeric"));
            } else {
                lhs = Value::number(f(lhs.number as i64, rhs.number as i64) as f64);
                lhs.is_true = lhs.number != 0.0;
            }
        }

        if undef {
            lhs.set_undefined();
        }

        Ok(lhs)
    }

    fn parse_add(&mut self) -> Result<Value, Error> {
        let mut lhs = self.parse_mul()?;

        loop {
            self.skip_ws();
            let Some(op) = self.consume_any(["+", "-"]) else {
                break;
            };

            let rhs = self.parse_mul()?;

            if !lhs.exists() || !rhs.exists() {
                lhs.set_undefined();
            } else if lhs.is_str || rhs.is_str {
                return Err(Error::new(self.pos, "arithmetic operands must be numeric"));
            } else if op == "+" {
                lhs.number += rhs.number;
                lhs.is_true = lhs.number != 0.0;
            } else {
                lhs.number -= rhs.number;
                lhs.is_true = lhs.number != 0.0;
            }
        }

        Ok(lhs)
    }

    fn parse_mul(&mut self) -> Result<Value, Error> {
        let mut lhs = self.parse_unary()?;

        loop {
            self.skip_ws();
            let Some(op) = self.consume_any(["*", "/", "%"]) else {
                break;
            };

            let rhs = self.parse_unary()?;

            if !lhs.exists() || !rhs.exists() {
                lhs.set_undefined();
            } else if lhs.is_str || rhs.is_str {
                return Err(Error::new(self.pos, "arithmetic operands must be numeric"));
            } else {
                match op {
                    "*" => lhs.number *= rhs.number,
                    "/" => lhs.number /= rhs.number,
                    "%" if rhs.number != 0.0 => {
                        lhs.number = (lhs.number as i64 % rhs.number as i64) as f64;
                    }
                    "%" => lhs.set_undefined(),
                    _ => unreachable!(),
                }

                if lhs.number.is_nan() {
                    lhs.set_undefined();
                } else {
                    lhs.is_true = lhs.number != 0.0;
                }
            }
        }

        Ok(lhs)
    }

    fn parse_unary(&mut self) -> Result<Value, Error> {
        self.skip_ws();

        if self.consume("+") {
            let mut value = self.parse_simple()?;
            if !value.exists() {
                value.set_undefined();
            } else if value.is_str {
                return Err(Error::new(self.pos, "unary operand must be numeric"));
            } else {
                value.is_true = value.number != 0.0;
            }
            Ok(value)
        } else if self.consume("-") {
            let mut value = self.parse_simple()?;
            if !value.exists() {
                value.set_undefined();
            } else if value.is_str {
                return Err(Error::new(self.pos, "unary operand must be numeric"));
            } else {
                value.number = -value.number;
                value.is_true = value.number != 0.0;
            }
            Ok(value)
        } else if self.consume("!") {
            let mut value = self.parse_unary()?;

            if value.is_true {
                value = Value::number(0.0);
            } else if !value.exists() {
                value = Value::boolean(true);
            } else if value.is_str {
                value = Value::boolean(value.string.is_none());
            } else {
                value = Value::boolean((value.number as i64) == 0);
            }

            Ok(value)
        } else if self.consume("~") {
            let mut value = self.parse_unary()?;

            if !value.exists() {
                value.set_undefined();
            } else if value.is_str {
                return Err(Error::new(self.pos, "unary operand must be numeric"));
            } else {
                value.number = !(value.number as i64) as f64;
                value.is_true = value.number != 0.0;
            }

            Ok(value)
        } else {
            self.parse_simple()
        }
    }

    fn parse_simple(&mut self) -> Result<Value, Error> {
        self.skip_ws();

        if self.consume("(") {
            let value = self.parse_expression()?;
            self.skip_ws();
            if !self.consume(")") {
                return Err(Error::new(self.pos, "missing ')'"));
            }
            return Ok(value);
        }

        if self.starts_with("\"") {
            return self.parse_string();
        }

        if let Some(value) = self.parse_number()? {
            return Ok(value);
        }

        let rest = &self.src[self.pos..];
        if let Some((value, len)) = (self.lookup)(rest) {
            self.pos += len;
            return Ok(value);
        }

        self.parse_function()
    }

    fn parse_function(&mut self) -> Result<Value, Error> {
        let start = self.pos;

        if self.consume("avg(") {
            let mut value = self.parse_expression()?;
            self.expect_function_close()?;

            if !value.is_str {
                return Err(Error::new(start, "avg expects a string"));
            }

            let Some(s) = value.string.take() else {
                value.set_undefined();
                return Ok(value);
            };

            let len = s.len();
            value = Value::number(if len == 0 {
                0.0
            } else {
                s.bytes().map(f64::from).sum::<f64>() / len as f64
            });
            return Ok(value);
        }

        if self.consume("default(") {
            let value = self.parse_expression()?;
            self.skip_ws();
            if !self.consume(",") {
                return Err(Error::new(self.pos, "missing ','"));
            }
            let default = self.parse_expression()?;
            self.expect_function_close()?;
            return Ok(if value.exists_true() { value } else { default });
        }

        if self.consume("exists(") {
            let value = self.parse_expression()?;
            self.expect_function_close()?;
            return Ok(Value::boolean(value.exists_true()));
        }

        if self.consume("exp(") {
            return self.parse_numeric_function(f64::exp);
        }

        if self.consume("length(") {
            let mut value = self.parse_expression()?;
            self.expect_function_close()?;

            if !value.is_str {
                return Err(Error::new(start, "length expects a string"));
            }

            value = match value.string.as_ref() {
                Some(s) => Value::number(s.len() as f64),
                None => Value::undefined(),
            };

            return Ok(value);
        }

        if self.consume("log(") {
            return self.parse_numeric_function(f64::ln);
        }

        if self.consume("min(") {
            return self.parse_string_byte_function(|bytes| bytes.iter().copied().min());
        }

        if self.consume("max(") {
            return self.parse_string_byte_function(|bytes| bytes.iter().copied().max());
        }

        if self.consume("pow(") {
            let lhs = self.parse_expression()?;
            self.skip_ws();
            if !self.consume(",") {
                return Err(Error::new(self.pos, "missing ','"));
            }
            let rhs = self.parse_expression()?;
            self.expect_function_close()?;

            if !lhs.exists() || !rhs.exists() {
                return Ok(Value::undefined());
            }

            if lhs.is_str || rhs.is_str {
                return Err(Error::new(start, "pow expects numbers"));
            }

            let value = lhs.number.powf(rhs.number);
            return Ok(if value.is_nan() {
                Value::undefined()
            } else {
                Value::number(value)
            });
        }

        if self.consume("sqrt(") {
            return self.parse_numeric_function(f64::sqrt);
        }

        Err(Error::new(start, "expected expression"))
    }

    fn parse_numeric_function(&mut self, f: fn(f64) -> f64) -> Result<Value, Error> {
        let value = self.parse_expression()?;
        self.expect_function_close()?;

        if !value.exists() {
            return Ok(Value::undefined());
        }

        if value.is_str {
            return Err(Error::new(self.pos, "function expects a number"));
        }

        let number = f(value.number);
        Ok(if number.is_nan() {
            Value::undefined()
        } else {
            Value::number(number)
        })
    }

    fn parse_string_byte_function(&mut self, f: fn(&[u8]) -> Option<u8>) -> Result<Value, Error> {
        let value = self.parse_expression()?;
        self.expect_function_close()?;

        if !value.is_str {
            return Err(Error::new(self.pos, "function expects a string"));
        }

        let Some(s) = value.string else {
            return Ok(Value::undefined());
        };

        Ok(f(s.as_bytes()).map_or_else(Value::undefined, |b| Value::number(f64::from(b))))
    }

    fn expect_function_close(&mut self) -> Result<(), Error> {
        self.skip_ws();
        if self.consume(")") {
            Ok(())
        } else {
            Err(Error::new(self.pos, "missing ')'"))
        }
    }

    fn parse_string(&mut self) -> Result<Value, Error> {
        self.consume("\"");
        let mut value = String::new();

        while !self.is_eof() {
            let Some(ch) = self.bump_char() else {
                break;
            };

            match ch {
                '"' => return Ok(Value::string(value)),
                '\\' => {
                    let Some(escaped) = self.bump_char() else {
                        value.push('\\');
                        break;
                    };

                    match escaped {
                        '"' => value.push('"'),
                        '\\' => value.push('\\'),
                        't' => value.push('\t'),
                        'n' => value.push('\n'),
                        'r' => value.push('\r'),
                        other => {
                            value.push('\\');
                            value.push(other);
                        }
                    }
                }
                other => value.push(other),
            }
        }

        Err(Error::new(self.pos, "unterminated string"))
    }

    fn parse_number(&mut self) -> Result<Option<Value>, Error> {
        let rest = &self.src[self.pos..];

        if let Some(hex) = rest
            .strip_prefix("0x")
            .or_else(|| rest.strip_prefix("0X"))
            .and_then(hex_prefix_len)
        {
            let token = &rest[2..2 + hex];
            let number =
                i64::from_str_radix(token, 16).map_err(|_| Error::new(self.pos, "bad hex"))?;
            self.pos += 2 + hex;
            return Ok(Some(Value::number(number as f64)));
        }

        let len = number_prefix_len(rest);
        if len == 0 {
            return Ok(None);
        }

        let token = &rest[..len];
        let number = token
            .parse::<f64>()
            .map_err(|_| Error::new(self.pos, "bad number"))?;
        self.pos += len;
        Ok(Some(Value::number(number)))
    }

    fn skip_ws(&mut self) {
        while let Some(ch) = self.src[self.pos..].chars().next() {
            if ch == ' ' || ch == '\t' {
                self.pos += ch.len_utf8();
            } else {
                break;
            }
        }
    }

    fn consume(&mut self, s: &str) -> bool {
        if self.starts_with(s) {
            self.pos += s.len();
            true
        } else {
            false
        }
    }

    fn consume_any<const N: usize>(&mut self, values: [&'static str; N]) -> Option<&'static str> {
        values.into_iter().find(|&value| self.consume(value))
    }

    fn starts_with(&self, s: &str) -> bool {
        self.src[self.pos..].starts_with(s)
    }

    fn bump_char(&mut self) -> Option<char> {
        let ch = self.src[self.pos..].chars().next()?;
        self.pos += ch.len_utf8();
        Some(ch)
    }

    fn is_eof(&self) -> bool {
        self.pos >= self.src.len()
    }
}

fn number_prefix_len(s: &str) -> usize {
    let bytes = s.as_bytes();
    let mut i = 0;
    let mut saw_digit = false;
    let mut saw_dot = false;

    while i < bytes.len() {
        match bytes[i] {
            b'0'..=b'9' => {
                saw_digit = true;
                i += 1;
            }
            b'.' if !saw_dot => {
                saw_dot = true;
                i += 1;
            }
            b'e' | b'E' if saw_digit => {
                let exp = i;
                i += 1;
                if i < bytes.len() && (bytes[i] == b'+' || bytes[i] == b'-') {
                    i += 1;
                }
                let digit_start = i;
                while i < bytes.len() && bytes[i].is_ascii_digit() {
                    i += 1;
                }
                if i == digit_start {
                    return exp;
                }
                return i;
            }
            _ => break,
        }
    }

    if saw_digit { i } else { 0 }
}

fn hex_prefix_len(s: &str) -> Option<usize> {
    let len = s.bytes().take_while(u8::is_ascii_hexdigit).count();
    (len > 0).then_some(len)
}

#[cfg(test)]
mod tests {
    use super::{Filter, Value};

    #[test]
    fn evaluates_arithmetic_and_null() {
        assert!(
            Filter::new("(1+2)*3==9")
                .eval_with(|_| None)
                .expect("eval")
                .truth()
        );

        let value = Filter::new("1+null").eval_with(test_lookup).expect("eval");
        assert!(!value.exists());
    }

    fn test_lookup(src: &str) -> Option<(Value, usize)> {
        match src {
            s if s.starts_with("null") => Some((Value::undefined(), 4)),
            _ => None,
        }
    }
}
