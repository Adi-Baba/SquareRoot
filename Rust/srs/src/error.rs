use thiserror::Error;

#[derive(Error, Debug)]
pub enum SrsError {
    #[error("Invalid input: {0}")]
    InvalidInput(String),

    #[error("Factorization failed: {0}")]
    FactorizationFailed(String),

    #[error("Precision too low: need at least {0} bits")]
    PrecisionTooLow(usize),

    #[error("Arithmetic error: {0}")]
    ArithmeticError(String),

    #[error("Certification failed: {0}")]
    CertificationFailed(String),
}

impl From<&str> for SrsError {
    fn from(s: &str) -> Self {
        SrsError::InvalidInput(s.to_string())
    }
}

impl From<String> for SrsError {
    fn from(s: String) -> Self {
        SrsError::InvalidInput(s)
    }
}

pub type Result<T> = std::result::Result<T, SrsError>;