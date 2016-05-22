# Copyright (c) 2016 NFLlib

set(CEREAL_PREFIX "" CACHE PATH "The path to the prefix of an CEREAL installation")

find_path(CEREAL_INCLUDE_DIR cereal
PATHS ${CEREAL_PREFIX}/include /usr/include /usr/local/include)

if(CEREAL_INCLUDE_DIR)
set(CEREAL_FOUND TRUE)
endif()