add_library(oxford oxford_overall.c detectionStage.c postProcessingStage.c preProcessingStage.c scoringStage.c ringbuffer.c)
target_include_directories(oxford PUBLIC include)