version = 1.9
os = alpine
osversion = 3.11

all:
	docker build -t overcookedfrog/samtools:$(version)-$(os)$(osversion) .
	docker tag overcookedfrog/samtools:$(version)-$(os)$(osversion) overcookedfrog/samtools:$(version)
	#docker tag overcookedfrog/samtools:$(version)-$(os)$(osversion) overcookedfrog/samtools:latest

.PHONY: push

push:
	docker push overcookedfrog/samtools:$(version)-$(os)$(osversion)
	docker push overcookedfrog/samtools:$(version)
	#docker push overcookedfrog/samtools:latest

