.DEFAULT_GOAL := storm-tracks

inla-from-scratch: inla-from-scratch.md
global-temperature: global-temperature.md
storm-tracks: storm-tracks.md

%.md: %.Rmd
	Rscript -e 'knitr::knit("$<")'

