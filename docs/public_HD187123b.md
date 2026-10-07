# The First Keck/HIRES Exoplanet

*J. Xavier Prochaska and Claude*

In the autumn of 1998, a short note appeared in the *Publications of the
Astronomical Society of the Pacific* announcing that a star in Cygnus was moving
towards and away from us, over and over, every three days. The star is
[HD 187123](https://simbad.cds.unistra.fr/simbad/sim-id?Ident=HD+187123), a G2V
dwarf about 46 parsecs — 150 light years — away. Its surface temperature is
5790 K; the Sun's is 5772 K. Its composition, mass, age and rotation rate are
all close to solar. If you wanted a stand-in for the Sun, you could do far
worse.

The thing making it move was a planet about half the mass of Jupiter, whipping
around it once every 3.097 days. The measurement was made with HIRES, the High
Resolution Echelle Spectrometer on the 10-metre Keck I telescope on Maunakea,
and it was the first planet the Keck telescope ever found.

The whole effect amounted to a speed of 72 metres per second — roughly the
cruising speed of a car on a motorway, measured across 150 light years of empty
space. This is the story of how that was possible, and why one small glass cell
full of iodine vapour mattered more than almost anything else in the building.

---

## Why a planet makes its star move

A planet does not orbit its star. The two of them orbit each other, about a
shared balance point called the centre of mass. The star is enormously heavier,
so it sits almost on top of that point and traces only a small circle — but it
does trace one.

![Left: a star and planet orbiting their common centre of mass, with the star's
motion greatly exaggerated. Right: the velocity an observer measures along the
sightline as the star goes around.](figs/fig1_wobble.png)

*Figure 1. Left: the star and the planet both circle their common centre of
mass. The planet itself is never seen — beside its star it is hopelessly faint —
so the star's motion is the entire measurement. Right: a spectrograph cannot see
the star's circle either. It records only the part of that motion aimed at us or
away from us, which is the circle's shadow cast along our sightline. When the
star is travelling across our view the shadow is momentarily still and there is
nothing to measure; when the star swings towards or away from Earth the shadow
moves fastest. Follow that shadow over time and you get the curve on the right.
The star's circle is drawn far larger than reality.*

The size of the star's circle follows from a single piece of bookkeeping. If the
planet has mass *m* and the star has mass *M*, then

> *m* × *v*<sub>planet</sub> = *M* × *v*<sub>star</sub>

Jupiter is about a thousandth of the Sun's mass, so the Sun moves about a
thousandth as fast as Jupiter does. For HD 187123, the planet is heavy and
extremely close in, which makes it move fast, which drags the star along at a
respectable 72 m/s.

We cannot see that circle directly. What we can measure is the part of the
motion pointed at us or away from us, through the Doppler effect. Light from a
receding star is stretched to longer wavelengths; light from an approaching star
is squeezed to shorter ones. The relation is simple algebra:

> Δλ / λ = *v* / *c*

where *c* is the speed of light. And it is at this point that the difficulty
becomes clear.

---

## The signal is one hundredth of a line

Take a typical absorption line in the star's spectrum — a narrow dip at a
wavelength of 550 nm — green light — where atoms in the star's atmosphere have
absorbed light. Plug 72 m/s into the equation above, and the line shifts by
0.00013 nm.

A stellar absorption line is about 0.015 nm wide. So the entire signal is roughly
one hundredth of the width of the feature you are trying to measure it against.

![Left: an absorption line at rest and shifted by 72 m/s, the two curves lying
on top of each other. Right: the same line magnified 60 times, where the shift
finally becomes visible.](figs/fig2_doppler.png)

*Figure 2. The shift caused by a 72 m/s motion, drawn to scale. On the left the
shifted and unshifted lines are indistinguishable. Only at 60× magnification, on
the steep flank of the line, does the entire signal become visible.*

Two things make this measurable at all. The first is that a stellar spectrum
contains not one line but thousands, and each carries the same shift, so the
information adds up. The second is that you do not need to know the absolute
wavelength — only whether it changed since last month.  Relative experiments
are the most sensitive in science.

That turns out to be the hard part, and it is where most attempts had failed.
Spectrographs are physical objects. They warm up and cool down, flex under their
own weight as the telescope moves, and drift as air pressure changes. A
temperature change of a fraction of a degree moves the spectrum across the
detector by far more than 0.00013 nm. Worse, the *shape* of the instrument's
response — how a single sharp wavelength gets smeared across several pixels,
known as the point spread function — changes too, and an asymmetric smear
imitates a Doppler shift almost perfectly.

Calibrating against a separate lamp does not solve this. Lamp light travels a
slightly different path through the instrument than starlight, and it is
recorded at a different moment. You end up calibrating something adjacent to
your measurement, rather than the measurement itself.

---

## The iodine cell

The fix is disarmingly physical: put a sealed glass cell of iodine vapour in the
beam, directly in front of the spectrograph, so the starlight passes through it
on its way in.

Molecular iodine has a forest of thousands of sharp absorption lines between
about 500 and 620 nm. Their wavelengths are fixed by molecular physics. They do
not drift, they do not care about the weather, and they will be in exactly the
same places in fifty years.

![Three stacked panels: the stellar spectrum, the dense iodine absorption
forest, and the product of the two, which is what the detector
records.](figs/fig3_iodine.png)

*Figure 3. The iodine cell imposes a reference spectrum -- a veritable ruler -- directly onto the
starlight. Every observation records the ruler and the thing being measured
together, through the same optics, at the same instant. This panel is an
illustration: the iodine forest shown here is synthetic, drawn to have the
character of the real thing rather than to reproduce it.*

This solves both problems at once, and that is the point worth dwelling on.

The iodine lines are a wavelength ruler that shares the starlight's exact
optical path, so any drift in the instrument moves the ruler and the spectrum
together and cancels out. And because the true iodine spectrum is already known
to extraordinary precision from laboratory measurements, any blurring or
asymmetry seen in the recorded iodine lines must have been introduced by the
instrument. The iodine spectrum therefore *measures the point spread function*
of that particular exposure — and once you know how the instrument smeared the
light, you can account for the smear instead of mistaking it for a Doppler
shift.

Analysing such a spectrum means modelling it: take a known iodine spectrum, take
a spectrum of the star recorded without the cell, shift the star by a trial
velocity, multiply them together, blur the result by a trial instrumental
profile, and compare to what was actually recorded. Adjust, and repeat. The
velocity that best reproduces the observation is the answer.

Geoff Marcy and Paul Butler had developed this approach at Lick Observatory and
in 1996 published the method paper that gave it its name —
[*Attaining Doppler Precision of 3 m/s*](https://doi.org/10.1086/133755). When
Steve Vogt, who designed and built HIRES, brought the technique to Keck, he
commissioned Marcy and Butler to build a custom iodine cell for it.

---

## Five nights in July

The Keck Planet Search began on 10 July 1996, with Marcy, Vogt and Butler
observing on time initially lent to them by Ben Zuckerman at UCLA. A new
velocity survey does not produce planets immediately; you need a few years to
build up enough observations of enough stars.

HD 187123 was not on the original list. It was added on the recommendation of
Kevin Apps, then an undergraduate at the University of Sussex, who had gone
through catalogues of nearby stars picking out metal-rich ones — a sensible bet,
since metal-rich stars turn out to host giant planets more often. He was right,
he was named as an author on the discovery paper, and the team took to calling
the result the "Planet of the Apps".

The first velocities of the star were taken in December 1997. Then, over five
consecutive nights in July 1998, ten observations caught the star swinging
through a complete cycle and back again — twice.

![Ten Keck/HIRES velocities of HD 187123 taken on five consecutive nights, 15 to
19 July 1998, tracking a 3.097-day cycle.](figs/fig4_july1998.png)

*Figure 4. Five nights in July 1998. Each point is one measurement of the star's
velocity towards or away from Earth: above the line it is approaching, below it
is receding. The star swings through about 125 m/s and back in roughly three
days, and then does it again. Once a run like this exists, the period is no
longer in doubt. Velocities are modern re-reductions of the original
discovery-epoch HIRES spectra, from the catalogue of Teklu et al. (2025); the
error bars, about 1.2 m/s, are smaller than the points.*

By the time the paper was submitted, 30 velocities had been collected between
December 1997 and September 1998. Folded on the 3.097-day period — that is, with
every observation plotted according to where it falls in the orbital cycle,
regardless of which month it was taken — they trace a clean curve.

![Thirty Keck/HIRES velocities of HD 187123 from December 1997 to September
1998, folded on the 3.097-day orbital period, following a smooth 69 m/s
curve.](figs/fig5_phasefold.png)

*Figure 5. All thirty discovery-epoch velocities, folded on the orbital period.
Observations taken nine months apart fall on the same curve, which is what makes
the interpretation a planet rather than a mood. Fitting these modern
re-reductions gives a semiamplitude of 69 m/s, against the 72 m/s published by
Butler and colleagues in 1998.*

The paper — [Butler, Marcy, Vogt & Apps
(1998)](https://doi.org/10.1086/316287) — appeared on 22 October 1998. It
reported a planet of at least 0.52 Jupiter masses on a near-circular orbit
0.042 AU from the star. The authors were careful to note that starspots and
pulsations could not be formally excluded, though the evidence pointed firmly at
a planet.

---

## A year that lasts three days

An orbit at 0.042 AU is about a ninth of Mercury's distance from the Sun, and
only about nine stellar radii from the surface of the star.

![The orbit of HD 187123 b drawn inside Mercury's orbit, showing how much
smaller it is.](figs/fig6_orbit_scale.png)

*Figure 6. HD 187123 b's orbit compared with Mercury's, to scale. The star is
drawn three times its true size to be visible at all.*

In 1998 this was still deeply strange. A gas giant cannot form that close to its
star, so it must have formed further out and migrated inwards — a conclusion the
first few hot Jupiters forced on a planet-formation community that had built its
theories around a single, orderly example.

## A note on "first"

There is a wrinkle worth being honest about. Seven weeks before the HD 187123
note, [Marcy and colleagues published the discovery of a planet around Gliese
876](https://doi.org/10.1086/311623), an M dwarf — and that paper does contain
Keck/HIRES velocities. But the detection rested on four years of data from Lick
Observatory; the authors state plainly that their single year of Keck data was
"not yet adequate" on its own.

So: Gliese 876 b is the first published planet paper containing HIRES data.
HD 187123 b is the first planet HIRES actually found. This document is about the
second.

## What happened next

HIRES kept watching. The catalogue now holds 137 velocities of HD 187123
spanning 24 years, and they revealed something invisible in 1998: a *second*
planet, HD 187123 c, roughly twice Jupiter's mass on a ten-year orbit. The same
instrument that caught a three-day wobble also held its calibration steadily
enough, across decades, to find a planet with a decade-long year. That is the
iodine cell's real legacy.

This repository is an attempt to go back to the beginning. The original 1998
spectra are archived, and [PypeIt](https://pypeit.readthedocs.io/) is a modern,
open-source reduction pipeline. We re-reduced those frames from scratch to see
how close we could come to the result that Butler, Marcy, Vogt and Apps
extracted from them a quarter of a century ago. The rest of this document is
what we found.

---

## Going back to the 1998 spectra

Every exposure Keck has ever taken is kept in the
[Keck Observatory Archive](https://koa.ipac.caltech.edu/). For HD 187123 in the
discovery era that is 36 exposures of the star, spread over 20 nights between
December 1997 and September 1998. There are also the lamp and flat-field
calibration frames taken on those nights. We did the work in three stages,
each built on the one before.

**Stage one: turning pictures into spectra.** An echelle spectrograph spreads
starlight into 37 short strips, stacked on the detector like the lines of a
page. Getting from the raw image to a spectrum means finding each strip,
correcting the detector's quirks and assigning a wavelength to every pixel.
PypeIt had never been run on HIRES data this old; the detector was replaced in
2004. Three bugs had to be fixed before the spectra came out right:
- the image was read 21 columns out of place;
- the map of bad detector columns landed 20 pixels from the real defects;
- the pixel binning was read backwards.

After that, all 20 nights reduced without any hand-tuning. The fixes are
public, and a report is going back to the PypeIt developers.

**Stage two: using a lamp as the ruler.** The obvious first attempt at
velocities uses the wavelength scale from the thorium–argon lamp exposures
taken each night. This is exactly the approach the iodine cell was invented
to replace, and the result shows why. Between exposures the velocities
scattered by **about 700 m/s**, ten times the planet's signal (Figure 7,
left). Fitting the known orbit to them gives an amplitude of 408 ± 188 m/s.
That is not a detection, and nothing in it would have looked different if the
planet did not exist.

The failure was informative, though. Within a single exposure, the different
strips of the spectrum agreed with each other to about 50 m/s. It was the
whole spectrum that moved, by up to 2 km/s, from one exposure to the next:
- **Drift through the night.** On one night the lamp-based wavelengths drifted
  by 1.6 km/s over seven hours, as the instrument flexed.
- **A borrowed lamp.** One night took no lamp exposures at all, so it borrowed
  the previous night's, and sat 2 km/s out.

This is the problem described in "The iodine cell" above, measured on real
data: a lamp calibrates something *adjacent* to the measurement, not the
measurement itself.

**Stage three: the iodine cell.** So we did what Butler and colleagues did,
and modelled every exposure. The model combines:
- a template spectrum of the star, taken without the cell;
- a laboratory spectrum of iodine;
- a trial velocity, a trial wavelength scale and a trial instrumental blur.

Each exposure was split into 700 short stretches, each fitted separately. We
used [`pyodine`](https://github.com/pepeheeren/pyodine), an open
implementation of the method by Heeren and colleagues (2023), adapting it to HIRES
and fixing what we found along the way. Two things we did not have:
- **the laboratory spectrum of the Keck cell itself.** We used one scanned
  from a different iodine cell in 2022;
- **the decades of refinement** that went into the pipeline built for this
  instrument.

![The same 31 exposures of HD 187123 folded on the orbital period, measured
three ways: with a lamp for the wavelength scale, with the iodine cell and
open software, and with the iodine cell and the modern professional
pipeline.](figs/fig7_three_ways.png)

*Figure 7. The same photons, measured three ways, folded on the 3.097-day
period. The grey curve is the orbit from the modern catalogue.*
- *Left: velocities from a lamp-based wavelength scale (stage two). Note the
  vertical scale, which is sixteen times larger than the other two panels;
  the shaded band is their entire range. The planet is lost.*
- *Middle: our iodine-cell measurements (stage three). The planet is plainly
  there.*
- *Right: the modern re-reduction of the same observations by the team that
  maintains the Keck velocity catalogue (Teklu et al. 2025), with error bars
  smaller than the points.*

**The planet came back.** Fitting a circular orbit at the known period to our
iodine velocities gives an amplitude of **72 ± 7 m/s**. Butler and colleagues
published 72 m/s in 1998, and the modern catalogue gives 69 m/s. That is a
detection at ten times its own uncertainty. It does not depend on knowing the
answer in advance:
- **The period, found blind.** Asked only "is there any repeating signal
  between half a day and a month?", the data single out **3.095 days**. The
  chance of a peak that strong arising from noise is about one in 40 billion.
- **The timing.** The planet's timing in its orbit agrees with the modern
  catalogue to within a few degrees.

**The iodine also repaired the lamp's failures.**
- **The borrowed-lamp night** sits on the orbit with the others.
- **The drift through the night** is gone: the seven-hour, 1.6 km/s drift
  disappears into the noise.

The iodine cell did the job it was built for. It carried its own ruler through
the instrument with the starlight, and the ruler made the lamp's problems
irrelevant.

## Where we fell short, and why

The honest comparison is Figure 8. Each of our iodine measurements scatters
by **26 m/s** about the orbit.
- That is well inside the 72 m/s signal, which is why the planet is detected.
- It is 28 times better than the lamp-based attempt.
- It is roughly **ten times worse** than the 3 m/s the method was designed to
  reach.
- It is roughly ten times worse than the **2.2 m/s** the modern pipeline gets
  from these same exposures.

![Horizontal bars comparing how much a single velocity measurement scatters
about the orbit for three approaches, against the 72 m/s signal of the
planet.](figs/fig8_precision.png)

*Figure 8. How much one velocity measurement scatters about the orbit, on a
logarithmic scale, for the three approaches in Figure 7. The dashed line is
the planet's signal. A method detects the planet when its bar ends well to the
left of the line.*

The shortfall is not the light itself. The number of photons in each exposure
would support about 1.5 m/s. What limits us is the modelling:
- **The wrong iodine cell.** Every iodine cell is slightly different, and the
  laboratory spectrum we used is of a different cell from the one in HIRES.
  Its lines are in the right places but not exactly the right shapes and
  depths, so the model never fits the data perfectly. Much of that mismatch is
  the same in every exposure and cancels when exposures are compared. The
  rest does not. A laboratory scan of the actual Keck cell, if one survives,
  would be the single biggest improvement.
- **A simpler model of the instrument.** Our description of how HIRES blurs a
  sharp line is deliberately simple. The Keck team refined theirs over many
  years.
- **The red end of the band.** At the long-wavelength end of the iodine
  region the iodine lines are weak. Those parts of the spectrum carry little
  information, and they add noise.
- **A few poor exposures.** In several frames, taken at the end of a night,
  only a small part of the spectrum was usable.

Leaving out the red end and the poorest exposures brings the scatter down to
about 12 m/s. But we decided to try that *after* seeing the results, which is
exactly the kind of choice that can flatter a measurement. We report 26 m/s as
our result.

So we did not reproduce what Butler, Marcy, Vogt and Apps achieved in 1998.
What we did show is that the discovery can be made again, independently, from
the archived photons, with open software and without any of the original
team's code. And the reason it works is the one that matters: a small glass
cell of iodine, put in the beam a quarter of a century ago.

Along the way we found problems in both of the open codes we relied on. They
go back to their authors: fifteen items for PypeIt and twenty-one for `pyodine`.
The full technical record, with every number in this section and the scripts
that produce them, is in this repository.

---

## References

- Butler, R. P., Marcy, G. W., Vogt, S. S., & Apps, K. 1998, *A Planet with a
  3.1 Day Period around a Solar Twin*, PASP, 110, 1389.
  [doi:10.1086/316287](https://doi.org/10.1086/316287)
- Butler, R. P., Marcy, G. W., Williams, E., et al. 1996, *Attaining Doppler
  Precision of 3 m/s*, PASP, 108, 500.
  [doi:10.1086/133755](https://doi.org/10.1086/133755)
- Marcy, G. W., Butler, R. P., Vogt, S. S., Fischer, D., & Lissauer, J. J. 1998,
  *A Planetary Companion to a Nearby M4 Dwarf, Gliese 876*, ApJ, 505, L147.
  [doi:10.1086/311623](https://doi.org/10.1086/311623)
- Wright, J. T., Marcy, G. W., Fischer, D. A., et al. 2007, *Four New Exoplanets
  and Hints of Additional Substellar Companions to Exoplanet Host Stars*, ApJ,
  657, 533. [doi:10.1086/510553](https://doi.org/10.1086/510553)
- Heeren, P., Tronsgaard, R., Grundahl, F., et al. 2023, A&A, 674, A164.
  [arXiv:2306.13615](https://arxiv.org/abs/2306.13615)
- Teklu, J. T., Perdelwitz, V., Butler, R. P., et al. 2025, *An updated catalog
  of HIRES/Keck radial velocity measurements*, A&A, 702, A68.
  [doi:10.1051/0004-6361/202555034](https://doi.org/10.1051/0004-6361/202555034)

Figures 4 and 5 use velocities from the Teklu et al. (2025) catalogue
([VizieR J/A+A/702/A68](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A%2BA/702/A68)),
retrieved 19 September 2026. These are modern, systematics-corrected
re-reductions of the original Keck/HIRES observations, not the velocity values
printed in the 1998 paper. Figures 7 and 8 add this project's own velocities
from the same exposures:
[`data/xcorr_velocities.csv`](../first_hires_exoplanet/data/xcorr_velocities.csv)
(lamp-calibrated) and
[`data/iodine_velocities.csv`](../first_hires_exoplanet/data/iodine_velocities.csv)
(iodine cell). All figures are generated by
[`first_hires_exoplanet/figs.py`](../first_hires_exoplanet/figs.py).
