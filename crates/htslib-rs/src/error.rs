//! Error helpers for preserving operation and path context.

use std::{
    error, fmt, io,
    path::{Path, PathBuf},
};

/// An error annotated with the operation that failed and, optionally, its path.
#[derive(Debug)]
pub struct ContextError {
    operation: &'static str,
    path: Option<PathBuf>,
    source: Box<dyn error::Error + Send + Sync + 'static>,
}

impl ContextError {
    /// Creates a contextual error for an operation.
    pub fn new<E>(operation: &'static str, source: E) -> Self
    where
        E: error::Error + Send + Sync + 'static,
    {
        Self {
            operation,
            path: None,
            source: Box::new(source),
        }
    }

    /// Creates a contextual error for an operation on a path.
    pub fn with_path<E, P>(operation: &'static str, path: P, source: E) -> Self
    where
        E: error::Error + Send + Sync + 'static,
        P: AsRef<Path>,
    {
        Self {
            operation,
            path: Some(path.as_ref().to_path_buf()),
            source: Box::new(source),
        }
    }

    /// Returns the operation that failed.
    pub fn operation(&self) -> &'static str {
        self.operation
    }

    /// Returns the path associated with the failed operation.
    pub fn path(&self) -> Option<&Path> {
        self.path.as_deref()
    }
}

impl fmt::Display for ContextError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.path {
            Some(path) => write!(
                f,
                "failed to {} {}: {}",
                self.operation,
                path.display(),
                self.source
            ),
            None => write!(f, "failed to {}: {}", self.operation, self.source),
        }
    }
}

impl error::Error for ContextError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        Some(self.source.as_ref())
    }
}

/// Adds context to ordinary result errors.
pub trait ResultExt<T> {
    /// Adds operation context to an error.
    fn with_context(self, operation: &'static str) -> Result<T, ContextError>;

    /// Adds operation and path context to an error.
    fn with_path_context<P>(self, operation: &'static str, path: P) -> Result<T, ContextError>
    where
        P: AsRef<Path>;
}

impl<T, E> ResultExt<T> for Result<T, E>
where
    E: error::Error + Send + Sync + 'static,
{
    fn with_context(self, operation: &'static str) -> Result<T, ContextError> {
        self.map_err(|e| ContextError::new(operation, e))
    }

    fn with_path_context<P>(self, operation: &'static str, path: P) -> Result<T, ContextError>
    where
        P: AsRef<Path>,
    {
        self.map_err(|e| ContextError::with_path(operation, path, e))
    }
}

/// Adds context while preserving an `io::Result` API.
pub trait IoResultExt<T> {
    /// Adds operation context to an I/O result.
    fn with_io_context(self, operation: &'static str) -> io::Result<T>;

    /// Adds operation and path context to an I/O result.
    fn with_io_path_context<P>(self, operation: &'static str, path: P) -> io::Result<T>
    where
        P: AsRef<Path>;
}

impl<T> IoResultExt<T> for io::Result<T> {
    fn with_io_context(self, operation: &'static str) -> io::Result<T> {
        self.map_err(|e| {
            let kind = e.kind();
            io::Error::new(kind, ContextError::new(operation, e))
        })
    }

    fn with_io_path_context<P>(self, operation: &'static str, path: P) -> io::Result<T>
    where
        P: AsRef<Path>,
    {
        self.map_err(|e| {
            let kind = e.kind();
            io::Error::new(kind, ContextError::with_path(operation, path, e))
        })
    }
}

/// Creates an `io::Error` with contextual source information.
pub fn io_error_with_path<E, P>(
    kind: io::ErrorKind,
    operation: &'static str,
    path: P,
    source: E,
) -> io::Error
where
    E: error::Error + Send + Sync + 'static,
    P: AsRef<Path>,
{
    io::Error::new(kind, ContextError::with_path(operation, path, source))
}
