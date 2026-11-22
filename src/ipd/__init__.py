"""Recovering what the published series hide: patients inside groups, and curves inside summaries.

    patient_level   the series that print their recipients one by one (Botha 2010, Chan 2011,
                    Ogura 2010): transcription checks and what each table supports
    ecological      the slope for a patient rather than for a group mean, by integrating
                    over the within-group spread of the pressure load
    km_reconstruct  individual survival data solved back out of a published Kaplan-Meier
                    curve and its numbers at risk (Guyot et al. 2012); validated on simulated
                    cohorts, not yet applied to a series
    run_ipd         command line

No module creates information: each uses what the original authors published. See
THEORY.md, section 8.
"""
