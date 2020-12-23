.PHONY: all

# run 'make' or 'make all' to render all webpages
all: index.html subpages/people.html \
	subpages/assignments.html subpages/seminars.html \
	subpages/lectures.html subpages/syllabus.html subpages/group_project_rubrics.html

# run 'make shallow' to render all webpages EXCEPT the seminars page	
shallow: index.html subpages/people.html \
	subpages/assignments.html  \
	subpages/lectures.html subpages/syllabus.html subpages/group_project_rubrics.html

# run 'make seminar' to grab seminars from 
# https://github.com/STAT540-UBC/STAT540-instructors-only/tree/master/seminars/seminars_winter_2019
# and remake the seminars webpage
seminar: subpages/seminars.html

deps: include/nav.html include/nothing.html

%.html: %.Rmd deps
	R -e "rmarkdown::render('$<')"

%.html: %.md deps
	R -e "rmarkdown::render('$<')"