This is the code for onezone freefall tests with grackle.

To run the code, you need several arguments:

The first argument chooses test problems
    0. Freefall test with various metalliciites.
       You need the second argument to choose metallicity
           0--8: 10^-{0--8} solar metallicities
           9   : primordial
       Ref. [Omukai et al. (2005)](https://ui.adsabs.harvard.edu/abs/2005ApJ...626..627O) 

    1. Freefall test for clouds exposed to strong dissociating photons
       that can eventually host supermassive stars.
       *** grackle has not been compatible with this test ***
       Ref. [Omukai et al. (2008)](https://ui.adsabs.harvard.edu/abs/2008ApJ...686..801O)

    2. Freefall test for clouds exposed to moderate dissociating photons.
       2nd argument: metallicity
       3rd argument: choose FUV intensity
           0--8: 10^-{0--8} G0
           9   : no FUV
       Ref. [Omukai (2012)](https://ui.adsabs.harvard.edu/abs/2012PASJ...64..114O)

    3. Freefall test for clouds enriched by a single Pop III SN.
       You can also include grain growth.
       2nd argument: metallicity
       3rd argument: choose metal/dust abundances and dust size distribution
           0: local ISM (Pollack et al. 1994)
           1--4: normal core-collapse SN with progenitor masses of
                 1: 13 Msun, 2: 20 Msun, 3: 25 Msun, 4: 30 Msun
           5--8: faint SN 
                 5: 13 Msun, 6: 15 Msun, 7: 50 Msun, 8: 80 Msun
           9--10: pair-instability SN
                 9: 170 Msun, 10: 200 Msun
           11: simple dust model (Yajima et al. 2019)
       Ref. [Marassi et al. (2014)](https://ui.adsabs.harvard.edu/abs/2014ApJ...794..100M)
       Ref. [Chiaki et al. (2015)](https://ui.adsabs.harvard.edu/abs/2015MNRAS.446.2659C)

    4. Freefall test for clouds exposed to interstellar radiation fields
       with a simple dust model used in [Park et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022ApJ...936..116P)
       2nd argument: metallicity

    5. Freefall test for clouds enriched by multiple Pop III SNe.
       In this code, you can give metals ejected by a 30 Msun CCSN (50%) and
       15 Msun Faint SN (50%).
       2nd argument: metallicity
