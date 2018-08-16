R_OPTS=--vanilla --no-save

pdf_with_diff: paper/latex/revision.tex
	cd paper/latex/;\
	pdflatex revision.tex;\
	bibtex revision;\
	pdflatex revision.tex;\
	pdflatex revision.tex;\
	latexdiff submitted.tex revision.tex > diff.tex;\
	pdflatex diff.tex;\
	bibtex diff;\
	pdflatex diff.tex;\
	pdflatex diff.tex

pdf: paper/latex/ref.bib paper/latex/revision.tex
	cd paper/latex/;\
	pdflatex revision.tex;\
	bibtex revision;\
	pdflatex revision.tex;\
	pdflatex revision.tex

sim_normal: sim/normal/nbiom10/simulate.R
	cd sim/normal/nbiom10/; nohup R CMD BATCH $(R_OPTS) simulate.R;\
	cd ..;	cd nbiom20/;	  nohup R CMD BATCH $(R_OPTS) simulate.R;\
	cd ..;	cd nbiom50/;	  nohup R CMD BATCH $(R_OPTS) simulate.R;\
	cd ..;	cd nbiom100/; 	nohup R CMD BATCH $(R_OPTS) simulate.R;\
	cd ..;	cd nbiom1000/; 	nohup R CMD BATCH $(R_OPTS) simulate.R;\
	cd ../..;	cd normal_block/nbiom10/; 	    nohup R CMD BATCH $(R_OPTS) simulate.R;\
	cd ../..;	cd normal_corr/nbiom10/; 	      nohup R CMD BATCH $(R_OPTS) simulate.R;\
	cd ..;	cd nbiom20/; 	nohup R CMD BATCH $(R_OPTS) simulate.R;\
	cd ..;	cd nbiom50/; 	nohup R CMD BATCH $(R_OPTS) simulate.R;\
	cd ..;	cd nbiom100/;	nohup R CMD BATCH $(R_OPTS) simulate.R;\
	cd ../..;	cd normal_tol.beta0/nbiom10/; 	nohup R CMD BATCH $(R_OPTS) simulate.R;\
	cd ..;	cd nbiom20/;  	nohup R CMD BATCH $(R_OPTS) simulate.R;\
	cd ..;	cd nbiom50/;  	nohup R CMD BATCH $(R_OPTS) simulate.R;\
	cd ..;	cd nbiom100/; 	nohup R CMD BATCH $(R_OPTS) simulate.R;\
	cd ..;	cd nbiom1000/; 	nohup R CMD BATCH $(R_OPTS) simulate.R;\


sim_survival:
	cd sim/survival/0predictive/;\
	nohup R CMD BATCH $(R_OPTS) "--args 100"  simulate.R simulate100.Rout;\
	nohup R CMD BATCH $(R_OPTS) "--args 400"  simulate.R simulate400.Rout;\
	nohup R CMD BATCH $(R_OPTS) "--args 700"  simulate.R simulate700.Rout;\
	nohup R CMD BATCH $(R_OPTS) "--args 1000" simulate.R simulate1000.Rout;\
	nohup R CMD BATCH $(R_OPTS) summary.R;\
	cd ..; cd 2predictive/;\
	nohup R CMD BATCH $(R_OPTS) "--args 100"  simulate.R simulate100.Rout;\
	nohup R CMD BATCH $(R_OPTS) "--args 400"  simulate.R simulate400.Rout;\
	nohup R CMD BATCH $(R_OPTS) "--args 700"  simulate.R simulate700.Rout;\
	nohup R CMD BATCH $(R_OPTS) "--args 1000" simulate.R simulate1000.Rout;\
	nohup R CMD BATCH $(R_OPTS) summary.R

sim:
	cd sim/normal/nbiom10/;\nohup R CMD BATCH $(R_OPTS) simulate.R;\
	cd ..;	cd nbiom20/;\	  nohup R CMD BATCH $(R_OPTS) simulate.R;\
	cd ..;	cd nbiom50/;\	  nohup R CMD BATCH $(R_OPTS) simulate.R;\
	cd ..;	cd nbiom100/;\ 	nohup R CMD BATCH $(R_OPTS) simulate.R;\
	cd ..;	cd nbiom1000/;\	nohup R CMD BATCH $(R_OPTS) simulate.R;\
	cd ../..;	cd normal_block/nbiom10/;\	    nohup R CMD BATCH $(R_OPTS) simulate.R;\
	cd ../..;	cd normal_corr/nbiom10/;\	      nohup R CMD BATCH $(R_OPTS) simulate.R;\
	cd ..;	cd nbiom20/;\ 	nohup R CMD BATCH $(R_OPTS) simulate.R;\
	cd ..;	cd nbiom50/;\ 	nohup R CMD BATCH $(R_OPTS) simulate.R;\
	cd ..;	cd nbiom100/;\	nohup R CMD BATCH $(R_OPTS) simulate.R;\
	cd ../..;	cd normal_tol.beta0/nbiom10/;\	nohup R CMD BATCH $(R_OPTS) simulate.R;\
	cd ..;	cd nbiom20/;\ 	nohup R CMD BATCH $(R_OPTS) simulate.R;\
	cd ..;	cd nbiom50/;\ 	nohup R CMD BATCH $(R_OPTS) simulate.R;\
	cd ..;	cd nbiom100/;\	nohup R CMD BATCH $(R_OPTS) simulate.R;\
	cd ..;	cd nbiom1000/;\	nohup R CMD BATCH $(R_OPTS) simulate.R;\
	cd ../..;\	cd survival/0predictive/;\
	nohup R CMD BATCH $(R_OPTS) "--args 100"  simulate.R simulate100.Rout;\
	nohup R CMD BATCH $(R_OPTS) "--args 400"  simulate.R simulate400.Rout;\
	nohup R CMD BATCH $(R_OPTS) "--args 700"  simulate.R simulate700.Rout;\
	nohup R CMD BATCH $(R_OPTS) "--args 1000" simulate.R simulate1000.Rout;\
	nohup R CMD BATCH $(R_OPTS) summary.R;\
	cd ..; cd 2predictive/;\
	nohup R CMD BATCH $(R_OPTS) "--args 100"  simulate.R simulate100.Rout;\
	nohup R CMD BATCH $(R_OPTS) "--args 400"  simulate.R simulate400.Rout;\
	nohup R CMD BATCH $(R_OPTS) "--args 700"  simulate.R simulate700.Rout;\
	nohup R CMD BATCH $(R_OPTS) "--args 1000" simulate.R simulate1000.Rout;\
	nohup R CMD BATCH $(R_OPTS) summary.R

figs: paper/figures/Rcode/Fig1.R paper/figures/Rcode/Fig2and3.R paper/figures/Rcode/Fig4and5.R paper/figures/Rcode/Fig6and7.R
	cd paper/figures/Rcode/;\
	nohup R CMD BATCH $(R_OPTS) Fig1.R;\
	nohup R CMD BATCH $(R_OPTS) Fig2and3.R;\
	nohup R CMD BATCH $(R_OPTS) Fig4and5.R;\
	nohup R CMD BATCH $(R_OPTS) Fig6and7.R;\
	cd ..; zip figures_pdf.zip *pdf; zip figures_eps.zip *eps

tables:	paper/tables/Rcode/Table1.R paper/tables/Rcode/Table2.R
	cd paper/tables/Rcode/;\
	nohup R CMD BATCH $(R_OPTS) Table1.R;\
	nohup R CMD BATCH $(R_OPTS) Table2.R

supplementary_srv:
	cd paper/SupplementaryMaterial/;\
	Rscript -e "Sys.setenv(PATH = '/usr/lib/rstudio-server/bin/pandoc/:$PATH');rmarkdown::render('SuppMaterial_v2.Rmd')"

supplementary:
	cd paper/SupplementaryMaterial/;\
	Rscript -e "rmarkdown::render('SuppMaterial_v2.Rmd')"

submission: paper/latex/revision.tex
	mkdir paper/submission
	cp paper/latex/revision.pdf paper/submission/
	cp paper/latex/diff.pdf paper/submission/
	cp paper/figures/*eps paper/submission/
	cp paper/SupplementaryMaterial/SuppMaterial_v2.pdf paper/submission/
	cd paper/;\zip -r submission submission/;\rm -r submission/

clean:
	find . -name "nohup.out" -type f -delete
	find . -name "Rplots.pdf" -type f -delete
