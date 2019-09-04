MESSAGE("** CURRENT SOURCE INCLUDE DIRECTORY: ${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND SOURCE_FILES
        "${CMAKE_CURRENT_LIST_DIR}/SinglePhase.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/SinglePhaseSym.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/SinglePhaseHelical.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/SinglePhaseSymX.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/SinglePhaseDirect.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/SeededSinglePhase.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/CEMapOptimizer.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/experimental.cpp"
        )

MESSAGE("** SOURCE FILES : ${SOURCE_FILES}")