# Source this file before running Nextflow on the assignment VM.
export JAVA_HOME="$HOME/.local/opt/java21"
export PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/compat:$JAVA_HOME/bin:$HOME/.local/bin:$PATH"
export NXF_VER=25.10.4
export R_LIBS_USER="$HOME/R/library"
