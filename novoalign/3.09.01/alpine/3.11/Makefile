name = novoalign
version = 3.09.01
os = alpine
osversion = 3.11
dockeruser = overcookedfrog

all: docker

docker: Dockerfile
	docker build -t $(dockeruser)/$(name):$(version)-$(os)$(osversion) .
	docker tag $(dockeruser)/$(name):$(version)-$(os)$(osversion) $(dockeruser)/$(name):$(version)

.PHONY: push clean

push:
	docker push $(dockeruser)/$(name):$(version)-$(os)$(osversion)
	docker push $(dockeruser)/$(name):$(version)

clean:
	docker rmi $(dockeruser)/$(name):$(version)-$(os)$(osversion)
	docker rmi $(dockeruser)/$(name):$(version)
