DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

export DELPHES_HOME="$DIR"
export PYTHONPATH="$DIR/python:${PYTHONPATH}"
export LD_LIBRARY_PATH="$DIR:${LD_LIBRARY_PATH}"
export LIBRARY_PATH="$DIR:${LIBRARY_PATH}"
