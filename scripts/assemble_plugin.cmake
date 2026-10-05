# Assemble only the files used by the Baltamatica plugin.
foreach(required HYMOS_SOURCE_DIR HYMOS_BUNDLE_DIR HYMOS_CORE_FILE HYMOS_CORE_SONAME HYMOS_VERSION)
    if(NOT DEFINED ${required} OR "${${required}}" STREQUAL "")
        message(FATAL_ERROR "Missing ${required}")
    endif()
endforeach()
get_filename_component(bundle_name "${HYMOS_BUNDLE_DIR}" NAME)
if(NOT IS_ABSOLUTE "${HYMOS_BUNDLE_DIR}" OR
   NOT bundle_name STREQUAL "HyMoS" OR
   NOT EXISTS "${HYMOS_BUNDLE_DIR}/main.so")
    message(FATAL_ERROR "Expected the built HyMoS plugin directory")
endif()
get_filename_component(bundle_real "${HYMOS_BUNDLE_DIR}" REALPATH)
get_filename_component(source_real "${HYMOS_SOURCE_DIR}" REALPATH)
if(bundle_real STREQUAL source_real OR IS_SYMLINK "${HYMOS_BUNDLE_DIR}")
    message(FATAL_ERROR "The bundle must be a separate build directory")
endif()

# Refresh managed directories so removed source assets cannot survive a rebuild.
foreach(subdir lib ui assets examples docs licenses scripts)
    set(destination "${HYMOS_BUNDLE_DIR}/${subdir}")
    if(IS_SYMLINK "${destination}")
        message(FATAL_ERROR "Refusing a symlink in the bundle: ${destination}")
    endif()
    get_filename_component(destination_real "${destination}" REALPATH)
    string(FIND "${destination_real}/" "${bundle_real}/" within_bundle)
    if(NOT within_bundle EQUAL 0)
        message(FATAL_ERROR "Bundle path escaped its directory: ${destination}")
    endif()
    file(REMOVE_RECURSE "${destination}")
endforeach()
file(MAKE_DIRECTORY "${HYMOS_BUNDLE_DIR}/lib" "${HYMOS_BUNDLE_DIR}/scripts")
file(COPY "${HYMOS_CORE_FILE}" DESTINATION "${HYMOS_BUNDLE_DIR}/lib")
get_filename_component(core_name "${HYMOS_CORE_FILE}" NAME)
file(CREATE_LINK "${core_name}" "${HYMOS_BUNDLE_DIR}/lib/${HYMOS_CORE_SONAME}" SYMBOLIC)
file(CREATE_LINK "${HYMOS_CORE_SONAME}" "${HYMOS_BUNDLE_DIR}/lib/libhymos.so" SYMBOLIC)
file(COPY "${HYMOS_SOURCE_DIR}/plugins/baltamatica/ui"
          "${HYMOS_SOURCE_DIR}/plugins/baltamatica/assets"
     DESTINATION "${HYMOS_BUNDLE_DIR}")
file(COPY "${HYMOS_SOURCE_DIR}/examples/ParallelPlate1D"
     DESTINATION "${HYMOS_BUNDLE_DIR}/examples")
foreach(directory docs licenses)
    if(EXISTS "${HYMOS_SOURCE_DIR}/${directory}")
        file(COPY "${HYMOS_SOURCE_DIR}/${directory}" DESTINATION "${HYMOS_BUNDLE_DIR}"
            PATTERN "*.aux" EXCLUDE PATTERN "*.log" EXCLUDE
            PATTERN "*.out" EXCLUDE PATTERN "*.toc" EXCLUDE
            PATTERN "*.bbl" EXCLUDE PATTERN "*.blg" EXCLUDE
            PATTERN "*.fls" EXCLUDE PATTERN "*.fdb_latexmk" EXCLUDE
            PATTERN "*.synctex.gz" EXCLUDE PATTERN "*.xdv" EXCLUDE)
    endif()
endforeach()
file(COPY "${HYMOS_SOURCE_DIR}/scripts/install.sh" DESTINATION "${HYMOS_BUNDLE_DIR}/scripts")
foreach(name README.md LICENSE LICENSE.md LICENSE.txt COPYING NOTICE NOTICE.md
             AUTHORS AUTHORS.md CITATION.cff AI_USAGE.md THIRD_PARTY_NOTICES.md)
    file(REMOVE "${HYMOS_BUNDLE_DIR}/${name}")
    if(EXISTS "${HYMOS_SOURCE_DIR}/${name}")
        file(COPY "${HYMOS_SOURCE_DIR}/${name}" DESTINATION "${HYMOS_BUNDLE_DIR}")
    endif()
endforeach()
foreach(document README.md docs/BUILD.md)
    set(document_path "${HYMOS_BUNDLE_DIR}/${document}")
    if(EXISTS "${document_path}")
        file(READ "${document_path}" document_text)
        string(REPLACE "(submission/README.md)"
            "(https://github.com/Liuyong2355/HyMoS/blob/main/submission/README.md)"
            document_text "${document_text}")
        string(REPLACE "(tests/README.md)"
            "(https://github.com/Liuyong2355/HyMoS/blob/main/tests/README.md)"
            document_text "${document_text}")
        string(REPLACE "(../tests/README.md)"
            "(https://github.com/Liuyong2355/HyMoS/blob/main/tests/README.md)"
            document_text "${document_text}")
        file(WRITE "${document_path}" "${document_text}")
    endif()
endforeach()
file(WRITE "${HYMOS_BUNDLE_DIR}/.hymos-plugin" "name=HyMoS\nformat=1\nversion=${HYMOS_VERSION}\n")
file(WRITE "${HYMOS_BUNDLE_DIR}/RUNTIME_DEPENDENCIES.txt"
"HyMoS ${HYMOS_VERSION} - Baltamatica Linux plugin
Bundled solver: lib/libhymos.so.5 (with its versioned library).
Host dependency: libbex.so is supplied by the Baltamatica installation.
System dependencies: glibc, libstdc++, libgcc_s, libm, libgomp (OpenMP).
The host SDK and system libraries are not redistributed in this package.
")
file(GLOB_RECURSE bundle_files LIST_DIRECTORIES false
     RELATIVE "${HYMOS_BUNDLE_DIR}" "${HYMOS_BUNDLE_DIR}/*")
list(SORT bundle_files)
file(WRITE "${HYMOS_BUNDLE_DIR}/SHA256SUMS" "")
foreach(name IN LISTS bundle_files)
    if(NOT name STREQUAL "SHA256SUMS")
        file(SHA256 "${HYMOS_BUNDLE_DIR}/${name}" digest)
        file(APPEND "${HYMOS_BUNDLE_DIR}/SHA256SUMS" "${digest}  ${name}\n")
    endif()
endforeach()
