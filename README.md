# Phobos-And-Deimos-Exploration

**LandingSitedV**
This file contains the calculations to "hop" between a series of landing sites on both of Mars's moons Phobos and Deimos. A series of points found through an overlayed latitude and longitude is established at the start of the file, with the points spanning most of the moons. From here I wrote a Lamberts solver to calculate the dV needed to fly from one landing spot to the next (from a start point on the surface to a midpoint inbetween and back down to the next landing site. To verify the results of the Lamberts solver I propagated the transfers with a two-body propagation with the velocities from the Lamberts solver. The propagated results were then plotted on a 3D map to visualize the transfers.
