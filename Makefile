.DEFAULT_GOAL := grid-regression

grid-regression: grid-regression.md
inla-from-scratch: inla-from-scratch.md
global-temperature: global-temperature.md
storm-tracks: storm-tracks.md
inla-lecture: inla-lecture.md

%.md: %.Rmd
	Rscript -e 'knitr::knit("$<")'

