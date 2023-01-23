exec_program(
    "git"
    ARGS "rev-parse" "HEAD"
    OUTPUT_VARIABLE VERSION )

message(STATUS ${VERSION})

configure_file(${SRC} ${DST})
