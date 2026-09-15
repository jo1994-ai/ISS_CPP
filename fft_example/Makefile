# Everything is done by build.sh (it builds the bundled FFTW when needed)

run:
	./build.sh

clean:
	rm -rf fft_example output

# also remove the compiled FFTW (it is rebuilt from third_party/ on the next run)
distclean: clean
	rm -rf third_party/fftw third_party/*.log

.PHONY: run clean distclean
