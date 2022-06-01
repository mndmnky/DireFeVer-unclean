use std::fmt;
use std::error::Error;

#[derive(Debug)]
pub enum ImportError {
    IoError(std::io::Error),
    InputMalformedError,
    BadIntError(std::num::ParseIntError),
}

impl From<std::io::Error> for ImportError {
    fn from(e: std::io::Error) -> ImportError {
        ImportError::IoError(e)
    }
}

impl From<std::num::ParseIntError> for ImportError {
    fn from(e: std::num::ParseIntError) -> ImportError {
        ImportError::BadIntError(e)
    }
}

impl fmt::Display for ImportError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::IoError(_) => write!(f, "Import: IoError"),
            Self::InputMalformedError => write!(f, "Import: Input is malformed."),
            Self::BadIntError(_) => write!(f, "Import: Integer is malformed."),
        }
    }
}

impl Error for ImportError {}

#[derive(Debug)]
pub enum GraphError {
    RebuildError,
    DeletionError,
    BrokenGraphError,
    InsertionError,
    NothingToRebuildError,
    PopPaceholderOutOfOrder,
}

impl fmt::Display for GraphError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::RebuildError => write!(f, "Graph: Graph could not be recovered"),
            Self::NothingToRebuildError => write!(f, "Graph: There is nothing to recover"),
            Self::DeletionError => write!(f, "Graph: Element could not be deleted."),
            Self::BrokenGraphError => write!(f, "Graph: The graph is broken."),
            Self::InsertionError => write!(f, "Graph: Element could not be inserted."),
            Self::PopPaceholderOutOfOrder => write!(f, "DFVSI: Popped placeholder is not the next in `merge_nodes`."),
        }
    }
}

impl Error for GraphError {}

#[derive(Debug)]
pub enum InterruptError {
    TimeOut,
    SigTerm,
    Unclear,
}

impl fmt::Display for InterruptError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TimeOut => write!(f, "Time ran out."),
            Self::SigTerm => write!(f, "SigTerm was send."),
            Self::Unclear => write!(f, "Interrupter was set."),
        }
    }
}

impl Error for InterruptError {}
