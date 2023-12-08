make:
	g++ main.cpp -lssl -lcrypto -o main.out

test: make
	./main.out