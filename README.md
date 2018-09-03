# PTA_GWBminmax

This reporsitory includes codes for the following paper "The minimum and maximum gravitational-wave background from supermassive binary black holes" (https://arxiv.org/abs/1806.02346) by Zhu, Cui and Thrane 2018. It deals with a single number Ayr, which is the signal amplitude of the gravitational wave background from supermassive binary black holes (GWB) at a reference frequency of 1 per year.

1. AyrNaivefun.m
This is the code to use when you have a new SMBBH candidate and want to know its implication for the GWB. rate and Anaive will tell you if it would imply a creazily low or high merger rate or Ayr.

2. AyrLowerfun.m
This does things more properly and outputs the 95% lower limit and median estimate of Ayr that includes Possionian uncertainties of the merger rate calculations.

3. TheMinimum_Table1.m
It includes some examples to use AyrLowerfun.m and AyrNaivefun.m and reproduces Table 1 and Figure 2 in the paper.

4. TheMaximum.m
This piece of code reproduces results in Section 3 of the paper. It takes a black hole mass function and tells you the maximum gravitational wave background that can be expected given the constraints on the number density of supermassive black holes in the local universe.

