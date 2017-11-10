.DEFAULT_GOAL := sampler-problem

grid-regression: grid-regression.md
inla-from-scratch: inla-from-scratch.md
global-temperature: global-temperature.md
storm-tracks: storm-tracks.md
sampler-problem: sampler-problem.md

%.md: %.Rmd
	Rscript -e 'knitr::knit("$<")'

