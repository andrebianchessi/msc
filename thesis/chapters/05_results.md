# Results

## PIM

- Time discretized into 8 buckets

![Dynamic response of a mass-spring-damper oscillator: PIM vs ETI prediction. Source: Author](figs/pimVsEtiSimpleSystem.png){#fig:pimVsEti width=80% style="scale:1;"}

- Time discretized into 20 buckets

![Dynamic response of a system of 6 masses, 9 springs and 6 dampers: PIM vs ETI prediction. Source: Author](figs/pimVsEti5Mass15LinksSystem.png){#fig:pimVsEti2 width=80% style="scale:1;"}

## E-GA vs P-GA vs Random {#sec:results_table}

The following table shows the *Efficiency* and *Quality* scores of the tests performed.
The first column indicates the number of non-fixed masses of the system, and the second indicates
the "id" of the experiment (3 random [COPs](#sec:cop) were solved for each number of non-fixed masses). See @sec:experiments for context.

| Masses | Test | Method       | Efficiency score | Quality score | Mean score |
| ------ | ---- | ------------ | ---------------- | ------------- | ---------- |
| 1      | 0    | P-GA         | 34.3654          | 71.0724       | 52.7189    |
| 1      | 0    | E-GA         | 0                | 100           | 50         |
| 1      | 0    | Random Guess | 100              | 0             | 50         |
|        |      |              |                  |               |            |
| 1      | 1    | P-GA         | 75.6739          | 0             | 37.837     |
| 1      | 1    | E-GA         | 0                | 100           | 50         |
| 1      | 1    | Random Guess | 100              | 0             | 50         |
|        |      |              |                  |               |            |
| 1      | 2    | P-GA         | 28.472           | 90.0197       | 59.2458    |
| 1      | 2    | E-GA         | 0                | 100           | 50         |
| 1      | 2    | Random Guess | 100              | 0             | 50         |
|        |      |              |                  |               |            |
| 2      | 0    | P-GA         | 72.3782          | 67.4949       | 69.9366    |
| 2      | 0    | E-GA         | 0                | 100           | 50         |
| 2      | 0    | Random Guess | 100              | 0             | 50         |
|        |      |              |                  |               |            |
| 2      | 1    | P-GA         | 47.5716          | 99.0683       | 73.3199    |
| 2      | 1    | E-GA         | 0                | 100           | 50         |
| 2      | 1    | Random Guess | 100              | 0             | 50         |
|        |      |              |                  |               |            |
| 2      | 2    | P-GA         | 62.6595          | 100           | 81.3298    |
| 2      | 2    | E-GA         | 0                | 98.1363       | 49.0681    |
| 2      | 2    | Random Guess | 100              | 0             | 50         |
|        |      |              |                  |               |            |
| 3      | 0    | P-GA         | 78.3921          | 60.019        | 69.2055    |
| 3      | 0    | E-GA         | 0                | 100           | 50         |
| 3      | 0    | Random Guess | 100              | 0             | 50         |
|        |      |              |                  |               |            |
| 3      | 1    | P-GA         | 67.7193          | 42.6812       | 55.2003    |
| 3      | 1    | E-GA         | 0                | 100           | 50         |
| 3      | 1    | Random Guess | 100              | 0             | 50         |
|        |      |              |                  |               |            |
| 3      | 2    | P-GA         | 79.8498          | 16.795        | 48.3224    |
| 3      | 2    | E-GA         | 0                | 100           | 50         |
| 3      | 2    | Random Guess | 100              | 0             | 50         |
|        |      |              |                  |               |            |
| 4      | 0    | P-GA         | 82.5093          | 0             | 41.2546    |
| 4      | 0    | E-GA         | 0                | 100           | 50         |
| 4      | 0    | Random Guess | 100              | 35.6312       | 67.8156    |
|        |      |              |                  |               |            |
| 4      | 1    | P-GA         | 64.3856          | 80.9193       | 72.6525    |
| 4      | 1    | E-GA         | 0                | 100           | 50         |
| 4      | 1    | Random Guess | 100              | 0             | 50         |
|        |      |              |                  |               |            |
| 4      | 2    | P-GA         | 65.827           | 70.105        | 67.966     |
| 4      | 2    | E-GA         | 0                | 100           | 50         |
| 4      | 2    | Random Guess | 100              | 0             | 50         |
|        |      |              |                  |               |            |
| 5      | 0    | P-GA         | 85.4621          | 27.0249       | 56.2435    |
| 5      | 0    | E-GA         | 0                | 100           | 50         |
| 5      | 0    | Random Guess | 100              | 0             | 50         |
|        |      |              |                  |               |            |
| 5      | 1    | P-GA         | 84.2239          | 23.8868       | 54.0554    |
| 5      | 1    | E-GA         | 0                | 100           | 50         |
| 5      | 1    | Random Guess | 100              | 0             | 50         |
|        |      |              |                  |               |            |
| 5      | 2    | P-GA         | 44.2273          | 15.1486       | 29.6879    |
| 5      | 2    | E-GA         | 0                | 100           | 50         |
| 5      | 2    | Random Guess | 100              | 0             | 50         |

## Conclusions {#sec:conclusions}