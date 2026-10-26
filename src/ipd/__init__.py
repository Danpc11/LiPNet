"""Recovering what the published series hide: patients inside groups, and curves inside summaries.

    km_reconstruct  individual survival data solved back out of a published Kaplan-Meier
                    curve and its numbers at risk (Guyot et al. 2012)
    ecological      the slope for a patient rather than for a group mean, by integrating
                    over the within-group spread of the pressure load

Neither module creates information. The first inverts a function the original authors
published; the second uses a dispersion that six of the groups already report. See
docs/IPD.md.
"""
