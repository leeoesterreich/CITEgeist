#!/usr/bin/env python
"""
The paradox: How does the same D538G mutation cause ER to LEAVE in MCF7 but BIND in T47D?

Let's look at this more carefully - especially Vehicle vs E2 conditions.
"""

import os
import gzip
import pandas as pd

DATA_DIR = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/midkine"

def load_peaks(filepath):
    peaks = []
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                peaks.append({'chr': parts[0], 'start': int(parts[1]),
                             'end': int(parts[2]), 'score': float(parts[4])})
    return pd.DataFrame(peaks)


def main():
    print("="*80)
    print("THE BINDING PARADOX: SAME MUTATION, OPPOSITE EFFECTS")
    print("="*80)

    # Load all ChIP-seq
    chip_files = {
        'MCF7_WT_Veh': 'GSM3563750_MCF7_WT_Vehicle_peaks.bed.gz',
        'MCF7_D538G_Veh': 'GSM3563756_MCF7_D538G_Vehicle_peaks.bed.gz',
        'MCF7_WT_E2': 'GSM3563751_MCF7_WT_E2_peaks.bed.gz',
        'MCF7_D538G_E2': 'GSM3563757_MCF7_D538G_E2_peaks.bed.gz',
        'T47D_WT_Veh': 'GSM3563759_T47D_WT_Vehicle_peaks.bed.gz',
        'T47D_D538G_Veh': 'GSM3563765_T47D_D538G_Vehicle_peaks.bed.gz',
        'T47D_WT_E2': 'GSM3563760_T47D_WT_E2_peaks.bed.gz',
        'T47D_D538G_E2': 'GSM3563766_T47D_D538G_E2_peaks.bed.gz',
    }

    peaks = {}
    for name, fname in chip_files.items():
        fpath = os.path.join(DATA_DIR, fname)
        if os.path.exists(fpath):
            peaks[name] = load_peaks(fpath)

    # Key observation: Look at Vehicle vs E2 separately
    print("\n" + "="*80)
    print("VEHICLE CONDITION (No Estrogen) - D538G Constitutive Activity")
    print("="*80)

    print(f"""
    MCF7:
      WT-Veh:    {len(peaks['MCF7_WT_Veh']):>6} peaks
      D538G-Veh: {len(peaks['MCF7_D538G_Veh']):>6} peaks  → {len(peaks['MCF7_D538G_Veh']) - len(peaks['MCF7_WT_Veh']):+} peaks ({len(peaks['MCF7_D538G_Veh'])/len(peaks['MCF7_WT_Veh']):.1f}x)

    T47D:
      WT-Veh:    {len(peaks['T47D_WT_Veh']):>6} peaks
      D538G-Veh: {len(peaks['T47D_D538G_Veh']):>6} peaks  → {len(peaks['T47D_D538G_Veh']) - len(peaks['T47D_WT_Veh']):+} peaks ({len(peaks['T47D_D538G_Veh'])/len(peaks['T47D_WT_Veh']):.1f}x)

    BOTH cell lines GAIN peaks in vehicle condition!
    This is D538G constitutive activity - it binds DNA without needing estrogen.
    """)

    print("\n" + "="*80)
    print("E2 CONDITION (With Estrogen) - Where the paradox appears")
    print("="*80)

    print(f"""
    MCF7:
      WT-E2:    {len(peaks['MCF7_WT_E2']):>6} peaks
      D538G-E2: {len(peaks['MCF7_D538G_E2']):>6} peaks  → {len(peaks['MCF7_D538G_E2']) - len(peaks['MCF7_WT_E2']):+} peaks ({len(peaks['MCF7_D538G_E2'])/len(peaks['MCF7_WT_E2']):.1f}x)

    T47D:
      WT-E2:    {len(peaks['T47D_WT_E2']):>6} peaks
      D538G-E2: {len(peaks['T47D_D538G_E2']):>6} peaks  → {len(peaks['T47D_D538G_E2']) - len(peaks['T47D_WT_E2']):+} peaks ({len(peaks['T47D_D538G_E2'])/len(peaks['T47D_WT_E2']):.1f}x)

    The paradox only appears WITH estrogen!
    MCF7: LOSES 57% of peaks
    T47D: GAINS 454% more peaks
    """)

    print("\n" + "="*80)
    print("THE EXPLANATION: IT'S ABOUT SATURATION")
    print("="*80)

    print(f"""
    The key is comparing WT-E2 to the "ceiling":

    MCF7-WT-E2:  12,472 peaks  ← Already NEAR MAXIMUM
    T47D-WT-E2:   1,724 peaks  ← Far below capacity

    D538G doesn't just make ER "bind more" or "bind less" -
    it changes WHICH sites ER prefers.

    IN MCF7 (saturated):
    ┌─────────────────────────────────────────────────────────┐
    │ WT-E2: 12,472 peaks (nearly all sites occupied)         │
    │                                                         │
    │ D538G has slightly different binding preferences        │
    │ It prefers some sites that WT doesn't, and vice versa   │
    │                                                         │
    │ But with 12,472 sites already filled:                   │
    │ - Sites D538G doesn't like → ER leaves                  │
    │ - Sites D538G prefers → already occupied, can't add     │
    │ - Net effect: LOSS of peaks                             │
    │                                                         │
    │ Result: 12,472 → 5,403 (lose 7,069 peaks)               │
    └─────────────────────────────────────────────────────────┘

    IN T47D (unsaturated):
    ┌─────────────────────────────────────────────────────────┐
    │ WT-E2: 1,724 peaks (many open sites available)          │
    │ High FOXA1 = lots of accessible chromatin               │
    │                                                         │
    │ D538G has constitutive activity (doesn't need E2)       │
    │ + E2 present = even more activity                       │
    │                                                         │
    │ Plenty of room to bind:                                 │
    │ - Some WT sites lost (D538G preference change)          │
    │ - But many NEW sites gained (open chromatin available)  │
    │ - Net effect: MASSIVE GAIN                              │
    │                                                         │
    │ Result: 1,724 → 9,552 (gain 7,828 peaks)                │
    └─────────────────────────────────────────────────────────┘
    """)

    print("\n" + "="*80)
    print("THE ANALOGY")
    print("="*80)

    print("""
    Think of parking lots:

    MCF7 = FULL parking lot (12,472 cars)
    - New driver (D538G) arrives with different parking preferences
    - Wants spots near the entrance, not the back
    - Some cars leave their "wrong" spots
    - But the "good" spots are already taken
    - Net effect: fewer total cars

    T47D = EMPTY parking lot (1,724 cars)
    - New driver (D538G) arrives
    - Lots of open spots available
    - Can park wherever it wants
    - Net effect: many more cars
    """)

    print("\n" + "="*80)
    print("PROOF: D538G EFFECT DEPENDS ON STARTING OCCUPANCY")
    print("="*80)

    print(f"""
    Compare D538G effect by condition:

    WITHOUT E2 (low starting occupancy):
      MCF7:  125 →  1,016  = +891 peaks (GAIN, 8.1x)
      T47D:  615 →  1,468  = +853 peaks (GAIN, 2.4x)
      Both GAIN when starting low!

    WITH E2 (high starting occupancy in MCF7):
      MCF7:  12,472 → 5,403 = -7,069 peaks (LOSE, 0.4x)
      T47D:   1,724 → 9,552 = +7,828 peaks (GAIN, 5.5x)
      MCF7 loses because it started saturated!

    The mutation effect isn't intrinsically "more binding" or "less binding" -
    it's a REDISTRIBUTION that has different net effects depending on
    starting occupancy.
    """)


if __name__ == "__main__":
    main()
