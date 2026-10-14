//! Crossover studies: the subject's course (from the fixed-duration engine) split into
//! treatment periods and washouts.
//!
//! Period `k` starts on its dosing day and ends the day after its in-house stay (when the
//! washout starts) or, for the last period attended, at the end of treatment. A subject who
//! stops during a washout ends in that washout.

use trialforge_core::Date;
use trialforge_model::Epoch;

use crate::timeline::{ElementSpan, Timeline};

/// Replaces the treatment element of `t` with period and washout elements. `period_days` are
/// the planned study days of each period's dose; `in_house_days` the days after a dose before
/// discharge.
pub fn apply(t: &mut Timeline, period_days: &[i32], in_house_days: i32) {
    let first = t.first_dose();
    let n = period_days.len();
    let starts: Vec<Date> = period_days
        .iter()
        .map(|d| first.add_days(d - 1))
        .filter(|d| *d <= t.treatment_end)
        .collect();
    let washout = |k: usize| Epoch::Washout(if n == 2 { 0 } else { k as u8 + 1 });
    let mut elements: Vec<ElementSpan> = t
        .elements
        .iter()
        .filter(|e| e.epoch == Epoch::Screening)
        .cloned()
        .collect();
    for (k, &start) in starts.iter().enumerate() {
        let discharge = start.add_days(in_house_days + 1);
        let epoch = Epoch::Period(k as u8 + 1);
        match starts.get(k + 1) {
            Some(&next) => {
                let end = discharge.min(next);
                elements.push(ElementSpan { epoch, start, end });
                if end < next {
                    elements.push(ElementSpan {
                        epoch: washout(k),
                        start: end,
                        end: next,
                    });
                }
            }
            None if k + 1 < n && t.treatment_end > discharge => {
                elements.push(ElementSpan {
                    epoch,
                    start,
                    end: discharge,
                });
                elements.push(ElementSpan {
                    epoch: washout(k),
                    start: discharge,
                    end: t.treatment_end,
                });
            }
            None => elements.push(ElementSpan {
                epoch,
                start,
                end: t.treatment_end,
            }),
        }
    }
    if t.last_contact > t.treatment_end {
        elements.push(ElementSpan {
            epoch: Epoch::FollowUp,
            start: t.treatment_end,
            end: t.last_contact,
        });
    }
    t.elements = elements;
}
