name = gtctovcf
version = 1.2.0
os = debian
osversion = 10-slim

all: docker

docker: Dockerfile
	docker build -t overcookedfrog/$(name):$(version)-$(os)$(osversion) .
	docker tag overcookedfrog/$(name):$(version)-$(os)$(osversion) overcookedfrog/$(name):$(version)
	docker tag overcookedfrog/$(name):$(version)-$(os)$(osversion) overcookedfrog/$(name):latest

.PHONY: push clean

push:
	docker push overcookedfrog/$(name):$(version)-$(os)$(osversion)
	docker push overcookedfrog/$(name):$(version)
	docker push overcookedfrog/$(name):latest

clean:
	docker rmi overcookedfrog/$(name):$(version)-$(os)$(osversion)
	docker rmi overcookedfrog/$(name):$(version)
	docker rmi overcookedfrog/$(name):latest
