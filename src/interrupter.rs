use std::time::Instant;
use crate::cust_errors::InterruptError;

#[derive(Debug, Default, Clone, Eq, PartialEq)]
pub struct Interrupter {
    start_time: Option<Instant>,
    time_limit: Option<u128>,
}

impl Interrupter {

    /// Creates a new Interrupter. 
    /// If `time_limit` is given, `self.check_interrupt()` becomes true after the `time_limit` has
    /// passed. 
    /// Otherwise only sigterm (int) events are considered (TODO: not yet implemented).
    pub fn new(time_limit: Option<u128>) -> Self {
        if time_limit.is_some() {
            let start = Instant::now();
            return Interrupter {
                start_time: Some(start),
                time_limit,
            }
        }
        Interrupter {
            start_time: None,
            time_limit: None,
        }
    }

    /// Checks if an sigterm (int?) was send, or the allowed time has expired. 
    /// TODO: Sigterm not yet implemented.
    ///
    /// On default this should always return false.
    pub fn check_interrupt(&self) -> bool {
        if let Some(dur) = self.time_limit {
            let duration = self.start_time.expect("since `time_limit` is some").elapsed();
            return duration.as_millis() >= dur
        }
        false
    }

    /// Sends an `InterruptError` if any interrupt was set.
    pub fn send_interrupt(&self) -> Result<(),InterruptError> {
        if let Some(dur) = self.time_limit {
            let duration = self.start_time.expect("since `time_limit` is some").elapsed();
            if duration.as_millis() >= dur {
                return Err(InterruptError::TimeOut);
            } 
            // TODO: check for sigterm
        }
        Ok(())
    }

}

