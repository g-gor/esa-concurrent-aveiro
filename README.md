## Orbit Transfer Optimization for Comet 311P
This repository contains MATLAB scripts for optimizing interplanetary orbit transfers, specifically focusing on a mission from Earth to Comet 311P. The core functionality revolves around the Lambert-Battin algorithm to determine optimal transfer trajectories.

```bash
main.m
```
This script utilizes the Lambert-Battin algorithm to calculate the optimal transfer date from Earth to Comet 311P, aiming to minimize the required ΔV (change in velocity).

Functionality: Determines the single best transfer date for a minimum ΔV trajectory.

Customization: You can easily adapt this script for other celestial bodies by adding their position and velocity vectors, initialized at 2032/12/31.

```bash
main_porkchops.m
```
This script generates a "porkchop plot," a powerful visualization tool for mission planning. It systematically sweeps through a range of launch and arrival dates to map the ΔV cost across a broad timeframe.

Functionality: Displays a contour plot showing the ΔV required for transfers as a function of both launch and arrival dates. This helps identify optimal launch windows.

Timeframes: The script currently sweeps through a 6-year launch window and a 6-year arrival window.

Adjustments: The lr_offset variable can be modified to adjust the minimum duration between launch and arrival, which is crucial for long-duration interplanetary missions.
