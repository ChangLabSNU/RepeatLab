#!/bin/bash
CONDA_DIR=$HOME/conda

if [ $# -lt 2 ]; then
  echo "Usage: $0 {env} {command} ..."
  exit 1
fi

ENV_NAME=$1
shift

__conda_setup="$($CONDA_DIR/bin/conda shell.bash hook 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "$CONDA_DIR/etc/profile.d/conda.sh" ]; then
        . "$CONDA_DIR/etc/profile.d/conda.sh"
    else
        export PATH="$CONDA_DIR/bin:$PATH"
    fi  
fi
unset __conda_setup

conda activate $ENV_NAME
eval "$@"
exit $?