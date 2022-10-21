# Primer Designer
## Evan Huang, Swinburne Lab, UC Berkeley

The purpose of this program is to simplify designing primers specifically for surrounding homology arms for gene editing. It utilizes Primer3 to attempt to find suitable primers surrounding a certain insert site in a sequence. Rather than needing to manually alter parameter values when no suitable primers are found, this program will automate that process by automatically relaxing certain variables until a set of suitable primers is found. 

These variables include melting temperature, GC%, number of GC clamps, number of poly-x regions, and product size. While not necessary, the user can specify search sites and variable orders (what order to relax variables). There are default values set for all of these variables, but the user is able to change any of them. 