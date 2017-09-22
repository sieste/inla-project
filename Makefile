inla-from-scratch: inla-from-scratch.md

inla-from-scratch.md: inla-from-scratch.Rmd
	Rscript -e 'knitr::knit("inla-from-scratch.Rmd")'


