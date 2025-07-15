# Randomized Switchback Method for Measurement and Verification #
## Project overview ##
Conventional measurement and verification (M&V) methods for estimating energy savings rely on comparing pre- and post-retrofit performance. They are often time-consuming and unreliable, especially when non-routine events, such as step changes or more gradual changes in building operation, occur during the M&V process. When those events are unrelated to the retrofit intervention but may significantly affect building energy consumption, then when the analyst applies the conventional M&V method the results can be confounded. In this study, we demonstrated that switchable interventions, such as most HVAC control retrofits, can benefit from a new M&V method that randomly samples whether to implement the baseline or the intervention strategy at a fixed interval (e.g., daily). We tested this novel randomized M&V method on a large public dataset (639 buildings) covering various climate zones and commercial building types, using a virtual chilled water supply temperature reset based on outdoor weather as the intervention. The results show that, compared to the conventional method, the randomized method provides more accurate savings estimations with a median of 74% accuracy improvement and is even faster (typically 36 weeks instead of 104 weeks, ~65% reduction). Additionally, we found that when non-routine events are present (e.g., occupancy pattern change), the randomized method estimates savings that are much closer to the ground-truth savings than the conventional method, demonstrating much improved reliability. We also assessed the impact of normalizing for different weather, starting the M&V at different dates of the year, continuing randomization with a different sampling ratio after satisfying all criteria, and dropping samples affected by carryover effects when switching between strategies. For each scenario, we identified the optimal sampling interval using the large dataset.

## Highlights ##
* Randomized switchback method leads to time reduction and increased accuracy
* Addressed the influence of non-routine events and also more subtle long-term 'operational drift' in buildings
* Investigated the trade-off of selecting longer sampling intervals to account for carryover effect and weather dependency of the method
* Discussed the impact of selecting different sampling ratio on the M&V accuracy

## Graphical summary ##
![graphical_abs](https://github.com/CenterForTheBuiltEnvironment/genome_mnv/blob/main/figs/manuscript/graphical_abs.jpg?raw=true)

## Structure ##
This repository includes: 
* [code](code/): R Markdown files for manuscript and supplementary material
* [readfiles](readfiles/): all input files and cleaned data
* [functions](functions/): defined functions for the analysis
* [paper](paper/): formatted journal submission
* [figs](figs/): generated figures for visualization