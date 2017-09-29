microSN: gerlumph_part.cpp
	nvcc -std=c++11 -Wno-deprecated-gpu-targets -o gerlumph_part -lcufft -lpng -lCCfits -lgerlumph -ljsoncpp gerlumph_part.cpp

clean:
	$(RM) gerlumph_part
