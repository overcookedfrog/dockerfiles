version = 2.29.2
os = alpine
osversion = 3.11

all: docker
	
docker: Dockerfile
	docker build -t overcookedfrog/bedtools:$(version)-$(os)$(osversion) .
	docker tag overcookedfrog/bedtools:$(version)-$(os)$(osversion) overcookedfrog/bedtools:$(version)
	docker tag overcookedfrog/bedtools:$(version)-$(os)$(osversion) overcookedfrog/bedtools:latest


.PHONY: clean push

push:
	docker push overcookedfrog/bedtools:$(version)-$(os)$(osversion)
	docker push overcookedfrog/bedtools:$(version)
	docker push overcookedfrog/bedtools:latest

clean:
	docker rmi overcookedfrog/bedtools:$(version)-$(os)$(osversion)
	docker rmi overcookedfrog/bedtools:$(version)
	docker rmi overcookedfrog/bedtools:latest
